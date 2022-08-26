force_recalculate_contributions = false

using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

## Saves contributions timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV

include(scriptsdir("include_src.jl"))

contributions_PAT = let excluded_artifact_grades = Int[],
    eeg = load_helsinki_eeg(PAT; 
        excluded_artifact_grades=excluded_artifact_grades
    ), #eeg = snip(eeg, 600, 975),
    snippets_duration=1, 
    preproc! = TripleCorrelations.zscore!, postproc! = TripleCorrelations.zscore!,
    assumption = IndStdNormal(), 
    conditioned_on = None(),
    lag_extents = (9,25), 
    plot_traces = true;

unique_id = if @isdefined(parent_session_id)
    parent_session_id
else
    Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")
end
target_match_str = "tricorr_ts_$(fn2str(preproc!))_$(fn2str(postproc!))_$(obj2str(assumption))_$(obj2str(conditioned_on))_snippets$(snippets_duration)_lagextents$(lag_extents[1])x$(lag_extents[2])_helsinkiEEG$(PAT)_artifacts_"
session_name = "$(target_match_str)$(unique_id)"
@show session_name

if preproc! == TripleCorrelations.zscore!
    preproc! = (o,i) -> TripleCorrelations.zscore!(o,i,mean(i),std(i))
end
maybe_jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
contributions = if isnothing(maybe_jld_dict) || force_recalculate_contributions
    calc_class_contributions(eeg, Periodic(), 
            preproc!, postproc!,
            assumption, conditioned_on
            ;
            lag_extents = lag_extents,
            n_motif_classes = 14,
            snippets_duration=snippets_duration
        )
else
    maybe_jld_dict["contributions"]
end

save(datadir("exp_pro", "$(session_name).jld2"), 
    Dict(
        "contributions" => contributions,
        "excluded_artifact_grades" => excluded_artifact_grades,
        "snippets_duration"=> snippets_duration,
        "lag_extents" => lag_extents,
        "assumption" => assumption,
        "conditioned_on" => conditioned_on,
        "preproc_str" => fn2str(preproc!),
        "postproc_str" => fn2str(postproc!)
    )
)

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

force_recalculate_contributions = false