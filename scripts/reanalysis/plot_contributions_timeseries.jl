using Distributed
using DrWatson
@quickactivate "NeonateTriCorr"

using Base.Threads
using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV
include(scriptsdir("include_src.jl"))

noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")

rms(xs) = sqrt(mean(xs .^ 2))

moving_average(vs, n) = [mean(@view vs[(i-n+1):i]) for i in n:length(vs)]

fig_eeg, fig_conts = let patients = [1, 47, 50];

snippets_duration_s = 1
save_prefix = "tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(snippets_duration_s)_lagextents8x25_helsinkiEEG"
plots_subdir = plotsdir("$(save_prefix)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))")
mkpath(plots_subdir)

target_match_patient = patient -> "tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(snippets_duration_s)_lagextents8x25_helsinkiEEG$(patient)_"
contributions_by_patient = map(patients) do patient
    jld_dict = load_most_recent_jld2(target_match_patient(patient), datadir("exp_pro"))
    contributions = jld_dict["contributions"]
end



fig_eeg = Figure(resolution=(3000,1600), font=noto_sans)
fig_conts = Figure(resolution=(3000,1600), font=noto_sans)

for (i, contributions) ∈ enumerate(contributions_by_patient)
    patient = patients[i]
    eeg = load_helsinki_eeg(patient)
    fig_eeg[1,i] = draw_eeg_traces!(fig_eeg, eeg; title = "EEG (Patient $(patient))", resolution=(1000,1600))
    times = get_times(eeg, sample_rate=1/snippets_duration_s)
    fig_conts[1,i] = plot_contributions!(fig_conts, eeg, times, contributions,; title="Motif Contributions (Patient $patient)", resolution=(1000,1600))
end
@show fig_conts
save(joinpath(plots_subdir, "example_eeg_traces.png"), fig_eeg)
save(joinpath(plots_subdir, "example_patients_contributions.png"), fig_conts)

# rms_timeseries = [rms(contributions[:,i_sec]) for i_sec ∈ 1:size(contributions,2)]
# lowpass_rms = moving_average(rms_timeseries, window)
# lowpass_times = times[window:end]
# (rms_fig, ax, l) = plot_contribution(eeg, lowpass_times, lowpass_rms;  resolution=(1200, 500), title="Patient $PAT, backward-lowpassed RMS motif contributions (blue seizure; red artifact; lowpass window = $(window)s)")

# save(joinpath(plots_subdir, "pat$(PAT)_rms_lowpass_$(window)s.png"), rms_fig)

# (rms_fig, conts_fig)
(fig_eeg, fig_conts)

end;
