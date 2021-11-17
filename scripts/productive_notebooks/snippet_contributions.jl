quickactivate("NeonateTriCorr")

include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics, StatsBase

using TripleCorrelations, ProgressMeter, Random, JLD2


exclude_channels_by_patient = Dict(
    75 => ["ECG","Resp","Cz","Fz"],
    50 => ["ECG","Resp","Cz"],
    47 => ["ECG","Resp","Cz","Fz"],
    44 => ["ECG","Resp","Cz"],
    31 => ["ECG","Resp","Cz"],
    9 => ["ECG","Resp","Cz"]
)

function calc_seizure_snippet_starts(annotations, snippets_duration)
    onsets, offsets = calc_seizure_bounds(annotations)
    bounds = zip(onsets, offsets)
    return vcat([collect(on:snippets_duration:off) for (on, off) ∈ bounds]...)
end

function calc_control_snippet_starts(annotations, snippets_duration, min_dist_to_seizure)
    seizure_onsets, seizure_offsets = calc_seizure_bounds(annotations)
    control_onsets = [0, seizure_offsets...]
    control_offsets = [seizure_onsets..., length(annotations) + 1]
    control_bounds = zip(control_onsets, control_offsets)
    return vcat([min(on+min_dist_to_seizure,off):snippets_duration:max(off-min_dist_to_seizure,on) for (on, off) ∈ control_bounds]...)
end

function calc_class_contributions(patient_number;  
        exclude_channels_by_patient=exclude_channels_by_patient,
        kwargs...)
    PAT = "eeg$(patient_number)"
    edf = EDF.read(datadir("exp_raw", "helsinki", "$(PAT).edf"))
    eeg = ProcessedEEG(edf; exclude=exclude_channels_by_patient[patient_number])

    calc_class_contributions(eeg; kwargs...) 
end

function calc_class_contributions(eeg::ProcessedEEGv1; λ_max = (9,25),
        n_motif_classes = 14, 
        n_seconds = floor(Int, eeg.duration),
        snippets_start_sec=0:(n_seconds-1),
        snippets_duration=1
    )
    eeg_motif_class_contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Float64, n_motif_classes, length(snippets_start_sec)))
    # sig_lock = ReentrantLock()
    # class_lock = ReentrantLock()
    p = ProgressMeter.Progress(length(snippets_start_sec))
    for (i_sec, snippet_start_sec) ∈ enumerate(snippets_start_sec)
        i_start = round(Int, (snippet_start_sec*eeg.sample_rate)+1)
        i_end = round(Int, (snippet_start_sec+snippets_duration)*eeg.sample_rate)
        snippet = eeg.signals[:,i_start:i_end]
        normalize_01!(snippet)
        @assert all(snippet .>= 0)
        actual_contributions = sequence_class_tricorr_zeropad(snippet, λ_max...)
        eeg_motif_class_contributions[:,i_sec] .= actual_contributions
        ProgressMeter.next!(p)
    end
    # jldsave(datadir("eeg_class_actual_$(λ_max)_$(PAT).jld2"); class_contributions=eeg_motif_class_contributions)

    fg_all = nothing#plot_contributions(eeg_motif_class_contributions; annotations=annotations, title=PAT)

    eeg_motif_class_contributions

end

function control_vs_seizure_class_contributions(patient_number; n_snippets, snippets_duration=1, min_dist_to_seizure=300, kwargs...)

    annotations = load_binary_annotations(patient_number)
    seizure_snippet_starts = calc_seizure_snippet_starts(annotations, snippets_duration)
    control_snippet_starts = calc_control_snippet_starts(annotations, snippets_duration, min_dist_to_seizure)

    if !(length(seizure_snippet_starts) > n_snippets && length(control_snippet_starts) > n_snippets)
        @warn "Patient $patient_number lacks sufficiently many separated snippets."
        return missing
    end

    seizure_snippet_starts = sort(sample(seizure_snippet_starts, n_snippets, replace=false))
    control_snippet_starts = sort(sample(control_snippet_starts, n_snippets, replace=false))

    seizure_snippets_class_contribution = calc_class_contributions(patient_number; snippets_start_sec = seizure_snippet_starts, snippets_duration = snippets_duration, kwargs...)
    control_snippets_class_contribution = calc_class_contributions(patient_number; snippets_start_sec = control_snippet_starts, snippets_duration = snippets_duration, kwargs...)

    return (seizure=seizure_snippets_class_contribution, control=control_snippets_class_contribution)
end

Random.seed!(1234)
res = map([31,44]) do pat_num#,31,44]
    contributions_by_case = control_vs_seizure_class_contributions(pat_num; n_snippets=3)
    seizure_contributions = [contributions.seizure for contributions in contributions_by_case]
    control_contributions = [contributions.control for contributions in contributions_by_case]
end
