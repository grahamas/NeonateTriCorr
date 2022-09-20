using DrWatson
@quickactivate "NeonateTriCorr"

## Saves contributions timeseries to datadir()/exp_pro/
## Excludes seconds marked as artifactual by

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

let min_reviewers_per_seizure = 2,
    excluded_artifact_grades = Int[1],
    epoch_s = 60,
    rolling_window_s = 60,
    snippets_duration_s = 15,
    patients_considered = [1:15..., 19,31,44,47,50,62];

params = Dict(
    :excluded_artifact_grades=>excluded_artifact_grades,
    :snippets_duration_s => snippets_duration_s,
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.zscore!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25)
)

rolling_window = Int(rolling_window_s / snippets_duration_s)
all_patient_results_df = load(datadir("aeeg_$(snippets_duration_s)_lower_margin_results_df_reviewers$(min_reviewers_per_seizure)_window$(rolling_window).jld2"))["channel_results_df"]

task_name = "detectaeeg"
save_root = make_filename_stem(task_name; params...)
save_dir = plotsdir("$(save_root)_reviewers$(min_reviewers_per_seizure)_grace$(epoch_s)_window$(rolling_window_s)_$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, patients_considered) do patient_num

    params[:patient_num] = patient_num

    results_df = filter(:patient => p -> p == patient_num, all_patient_results_df)

    @assert !isempty(excluded_artifact_grades)
    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades)
    aeeg = calculate_aeeg(eeg; window_len_s=snippets_duration_s)
    signals = aeeg_lower_margin(aeeg)'
    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    fig = plot_μ_and_σ_signals_and_roc(results_df, signals, signal_times, seizure_bounds; analysis_eeg=eeg, rolling_window_s=rolling_window_s, example_θ=3, n_signals_used=5, epoch_s=epoch_s, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)")

    save(joinpath(save_dir, "$(task_name)_roc_patient$(patient_num)_reviewers$(min_reviewers_per_seizure)$(artifacts_str(excluded_artifact_grades)).png"), fig)
end

end