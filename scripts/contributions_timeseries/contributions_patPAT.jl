using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

## Saves contributions timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV
include(scriptsdir("include_src.jl"))

rms(xs) = sqrt(mean(xs .^ 2))
moving_average(vs, n) = [mean(@view vs[(i-n+1):i]) for i in n:length(vs)]

PAT = 1

contributions_PAT = let eeg = load_helsinki_eeg(PAT),# eeg = snip(eeg, 600, 975),
    snippets_duration=1, preproc! = zscore!, postproc! = zscore!,
    assumption = IndStdNormal(), conditioned_on = None(),
    lag_extents = (8,25), plot_traces=true;

unique_id = if @isdefined(parent_session_id)
    parent_session_id
else
    Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")
end
session_name = "tricorr_ts_$(fn2str(preproc!))_$(fn2str(postproc!))_$(obj2str(assumption))_$(obj2str(conditioned_on))_snippets$(snippets_duration)_lagextents$(lag_extents[1])x$(lag_extents[2])_helsinkiEEG$(PAT)_$(unique_id)"
@show session_name

if preproc! == zscore!
    preproc! = (o,i) -> zscore!(o,i,mean(i),std(i))
end
contributions = calc_class_contributions(eeg, Periodic(), 
        preproc!, postproc!,
        assumption, conditioned_on
        ;
        lag_extents = lag_extents,
        n_motif_classes = 14,
        snippets_duration=snippets_duration
    )

save(datadir("exp_pro", "$(session_name).jld2"), Dict("contributions" => contributions))

# contributions = load(datadir("exp_pro", "$(session_name).jld2"))["contributions"]

if plot_traces
    plots_subdir = plotsdir(session_name)
    mkpath(plots_subdir)

    eeg_fig = draw_eeg_traces(eeg; title = "EEG (Patient $PAT)", resolution=(1000,1600))
    save(joinpath(plots_subdir, "pat$(PAT)_eeg_traces.png"), eeg_fig)

    times = get_times(eeg, sample_rate=snippets_duration)
    conts_fig = plot_contributions(eeg, times, contributions; title="Motif Contributions (Patient $PAT)", resolution=(1000,1600))
    save(joinpath(plots_subdir, "pat$(PAT)_contributions.png"), conts_fig)
end

# rms_timeseries = [rms(contributions[:,i_sec]) for i_sec ∈ 1:size(contributions,2)]
# lowpass_rms = moving_average(rms_timeseries, window)
# lowpass_times = times[window:end]
# (rms_fig, ax, l) = plot_contribution(lowpass_times, lowpass_rms; eeg=eeg, resolution=(1200, 500), title="Patient $PAT, backward-lowpassed RMS motif contributions (blue seizure; red artifact; lowpass window = $(window)s)")

# save(joinpath(plots_subdir, "pat$(PAT)_rms_lowpass_$(window)s.png"), rms_fig)

# eeg_fig = draw_eeg_traces(eeg; title = "Patient $PAT")
# save(joinpath(plots_subdir, "pat$(PAT)_eeg_traces.png"), eeg_fig)

# (rms_fig, conts_fig)

contributions

end

# single thread:   66.249 s (90203946 allocations: 8.51 GiB)
