using Distributed
using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"

using Random, JLD2

using HypothesisTests, CSV

rms(xs) = sqrt(mean(xs .^ 2))

let eeg = load_helsinki_eeg(PAT);

contributions = calc_class_contributions(eeg, Periodic(), AN_01norm;
        λ_max = (8,25),
        n_motif_classes = 14,
        snippets_duration=1
    )

save(datadir("exp_pro", "timeseries_AN_pat$(PAT)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).jld2"), Dict("contributions" => contributions))

plots_subdir = plotsdir("timeseries_AN_pat$(PAT)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))")
mkpath(plots_subdir)

eeg_fig = plot_eeg_traces(eeg)

contributions_fig = plot_contributions(contributions; eeg=eeg)
save(joinpath(plots_subdir, "timeseries_AN_allmotifs_pat$(PAT)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).$(ext)"), contributions_fig)

rms_timeseries = [rms(contributions[:,i_sec]) for i_sec ∈ 1:size(contributions,2)]
rms_fig = plot_contributions(rms_timeseries; eeg=eeg, n_motif_classes=1)
save(joinpath(plots_subdir, "timeseries_AN_RMS_pat$(PAT)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).$(ext)"), rms_fig)

lowpass = digitalfilter(Lowpass(0.2; fs=eeg.sample_rate), Butterworth(2))
lowpass_rms = filt(lowpass, rms_timeseries) # want row
lowpass_rms_fig = plot_contributions(lowpass_rms; eeg=eeg, n_motif_classes=1)
save(joinpath(plots_subdir, "timeseries_AN_lowpassRMS_pat$(PAT)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).$(ext)"), lowpass_rms_fig)

end