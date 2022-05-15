using Base.Threads
using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV
include(scriptsdir("include_src.jl"))

rms(xs) = sqrt(mean(xs .^ 2))

moving_average(vs, n) = [mean(@view vs[(i-n+1):i]) for i in n:length(vs)]

let eeg = load_helsinki_eeg(PAT), eeg = snip(eeg, 0, 0+300), 
    window=30, contributions_desc = "A_znorm",
    snippets_duration=1;
    λ_max = (8,25)
contributions = calc_class_contributions(eeg, Periodic(), snippet_contributions_fns[contributions_desc];
        λ_max = λ_max,
        n_motif_classes = 14,
        snippets_duration=snippets_duration
    )

save(datadir("exp_pro", "timeseries_$(contributions_desc)_$(λ_max[1])_$(λ_max[2])_pat$(PAT)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).jld2"), Dict("contributions" => contributions))

plots_subdir = plotsdir("timeseries_$(contributions_desc)_pat$(PAT)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))")
mkpath(plots_subdir)

eeg_fig = draw_eeg_traces(eeg; title = "EEG (Patient $PAT)", resolution=(1000,1600))
save(joinpath(plots_subdir, "pat$(PAT)_eeg_traces.png"), eeg_fig)

times = get_times(eeg, sample_rate=snippets_duration)
conts_fig = plot_contributions(times, contributions, eeg; title="Motif Contributions (Patient $PAT)", resolution=(1000,1600))
save(joinpath(plots_subdir, "pat$(PAT)_contributions.png"), conts_fig)

rms_timeseries = [rms(contributions[:,i_sec]) for i_sec ∈ 1:size(contributions,2)]
lowpass_rms = moving_average(rms_timeseries, window)
lowpass_times = times[window:end]
(rms_fig, ax, l) = plot_contribution(lowpass_times, lowpass_rms; eeg=eeg, resolution=(1200, 500), title="Patient $PAT, backward-lowpassed RMS motif contributions (blue seizure; red artifact; lowpass window = $(window)s)")

save(joinpath(plots_subdir, "pat$(PAT)_rms_lowpass_$(window)s.png"), rms_fig)

# eeg_fig = draw_eeg_traces(eeg; title = "Patient $PAT")
# save(joinpath(plots_subdir, "pat$(PAT)_eeg_traces.png"), eeg_fig)

(rms_fig, conts_fig)

end

# single thread:   66.249 s (90203946 allocations: 8.51 GiB)
