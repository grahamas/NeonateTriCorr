using Distributed
using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

using Base.Threads
using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV
include(scriptsdir("include_src.jl"))


rms(xs) = sqrt(mean(xs .^ 2))

moving_average(vs, n) = [mean(@view vs[(i-n+1):i]) for i in n:length(vs)]



plots_subdir = plotsdir("AN_twente_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))")
mkpath(plots_subdir)

variability = zeros(Float64, 50)
average_amp = zeros(Float64, 50)

let discard_first_len=0, window=30, contributions_sampling_rate=1;
    
eeg_name = "rolandic"# "EEG$(lpad(PAT,3,"0"))"

eeg = snip_start(load_twente_eeg(eeg_name), discard_first_len)
processed_data_filenames = readdir(datadir("exp_pro"))
timeseries_filenames = processed_data_filenames[startswith.(processed_data_filenames, "timeseries_AN_8_25_twente")]

eeg_fig = draw_eeg_traces(eeg; title = "$eeg_name", resolution=(1000,1600))
save(joinpath(plots_subdir, "$(eeg_name)_eeg_traces.png"), eeg_fig)

target_filenames = timeseries_filenames[occursin.("$(eeg_name)_", timeseries_filenames)]
DATA = load(datadir("exp_pro", target_filenames[end]))
contributions = set_artifacts_missing(DATA["contributions"][:, discard_first_len+1:end], eeg, sample_rate=contributions_sampling_rate)

times = get_times(eeg, sample_rate=contributions_sampling_rate)
conts_fig = plot_contributions(eeg, times, contributions,; title="Motif Contributions ($(eeg_name))", resolution=(1000,1600))
save(joinpath(plots_subdir, "$(eeg_name)_contributions.png"), conts_fig)

rms_timeseries = [rms(contributions[:,i_sec]) for i_sec ∈ 1:size(contributions,2)]
(rms_fig, ax, l) = plot_contribution(times, rms_timeseries; eeg=eeg, resolution=(1200, 500), title="$(eeg_name), RMS motif contributions (blue seizure; red artifact)")
save(joinpath(plots_subdir, "$(eeg_name)_rms.png"), rms_fig)


lowpass_rms = moving_average(rms_timeseries, window)
lowpass_times = times[window:end]
(rms_fig, ax, l) = plot_contribution(lowpass_times, lowpass_rms; eeg=eeg, resolution=(1200, 500), title="$(eeg_name), backward-lowpassed RMS motif contributions (blue seizure; red artifact; lowpass window = $(window)s)")

save(joinpath(plots_subdir, "$(eeg_name)_rms_lowpass_$(window)s.png"), rms_fig)

# (rms_fig, conts_fig)
end#let