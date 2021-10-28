quickactivate("NeonateTriCorr")

include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics

using TripleCorrelations, ProgressMeter, Random, JLD2


exclude_channels_by_patient = Dict(
    75 => ["ECG","Resp","Cz","Fz"],
    50 => ["ECG","Resp","Cz"],
    47 => ["ECG","Resp","Cz","Fz"],
    44 => ["ECG","Resp","Cz"],
    31 => ["ECG","Resp","Cz"],
    9 => ["ECG","Resp","Cz"]
)

function eval_class_contributions(patient_number; exclude_channels_by_patient=exclude_channels_by_patient)
PAT = "eeg$(patient_number)"
edf = EDF.read(datadir("exp_raw", "helsinki", "$(PAT).edf"))
eeg = ProcessedEEG(edf; exclude=exclude_channels_by_patient[patient_number])

#eeg_fig = draw_eeg_traces(eeg)
#save(plotsdir("traces", "traces_$(PAT).png"), eeg_fig)

λ_max = (9,25)
n_motif_classes = 14

n_seconds = floor(Int, eeg.duration)
eeg_motif_class_contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Float64, n_motif_classes, n_seconds-1))
sig_lock = ReentrantLock()
class_lock = ReentrantLock()
p = ProgressMeter.Progress(n_seconds-1)
Threads.@threads for i_sec ∈ 1:n_seconds-1
    i_start = round(Int, (i_sec-1)*eeg.sample_rate+1)
    i_end = round(Int, i_sec*eeg.sample_rate)
    snippet = lock(sig_lock) do
        eeg.signals[:,i_start:i_end]
    end
    normalize_01!(snippet)
    @assert all(snippet .>= 0)
    actual_contributions = sequence_class_tricorr_zeropad(snippet, λ_max...)
    lock(class_lock) do
        eeg_motif_class_contributions[:,i_sec] .= actual_contributions
    end
    ProgressMeter.next!(p)
end
jldsave(datadir("eeg_class_actual_SHORT_$(λ_max)_$(PAT).jld2"); class_contributions=eeg_motif_class_contributions)

# Clear bad values
annotations = load_binary_annotations(patient_number)[:]

fg_all = nothing#plot_contributions(eeg_motif_class_contributions; annotations=annotations, title=PAT)

(fg_all, eeg_motif_class_contributions, eeg)

end


res = map([31,44]) do pat_num#,31,44]
    eval_class_contributions(pat_num)
end
