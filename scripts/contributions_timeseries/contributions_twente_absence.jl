using Distributed
using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV
using Base.Threads
include(scriptsdir("include_src.jl"))

rms(xs) = sqrt(mean(xs .^ 2))

moving_average(vs, n) = [mean(@view vs[(i-n+1):i]) for i in n:length(vs)]

let eeg_name = EEG_NAME, eeg = load_twente_eeg(eeg_name), window=30,
    snippets_duration=1, 
    λ_max = (8,25);
@info "Calculating contributions..."
contributions = calc_class_contributions(eeg, Periodic(), AN_01norm;
        λ_max = λ_max,
        n_motif_classes = 14,
        snippets_duration=snippets_duration
    )
@info "Saving data..."
save(datadir("exp_pro", "timeseries_AN_$(λ_max[1])_$(λ_max[2])_twente_$(eeg_name)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).jld2"), Dict("contributions" => contributions))

plots_subdir = plotsdir("timeseries_AN_twente_$(eeg_name)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))")
mkpath(plots_subdir)

eeg_fig = draw_eeg_traces(eeg; title = "EEG ($(eeg_name))", resolution=(1000,1600))
save(joinpath(plots_subdir, "$(eeg_name)_eeg_traces.png"), eeg_fig)

times = get_times(eeg, sample_rate=snippets_duration)
conts_fig = plot_contributions(times, contributions, eeg; title="Motif Contributions ($(eeg_name))", resolution=(1000,1600))
save(joinpath(plots_subdir, "$(eeg_name)_contributions.png"), conts_fig)

rms_timeseries = [rms(contributions[:,i_sec]) for i_sec ∈ 1:size(contributions,2)]
lowpass_rms = moving_average(rms_timeseries, window)
lowpass_times = times[window:end]
(rms_fig, ax, l) = plot_contribution(lowpass_times, lowpass_rms; eeg=eeg, resolution=(1200, 500), title="$(eeg_name), backward-lowpassed RMS motif contributions (blue seizure; red artifact; lowpass window = $(window)s)")

save(joinpath(plots_subdir, "$(eeg_name)_rms_lowpass_$(window)s.png"), rms_fig)

# eeg_fig = draw_eeg_traces(eeg; title = "$(eeg_name)")
# save(joinpath(plots_subdir, "$(eeg_name)_eeg_traces.png"), eeg_fig)

(rms_fig, conts_fig)

end

# single thread:   66.249 s (90203946 allocations: 8.51 GiB)
