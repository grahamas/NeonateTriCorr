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

function detect_patient_seizures(patient_num; save_dir,
        excluded_artifact_grades, min_reviewers_per_seizure, snippets_duration_s,
        task_name
    )
    @assert !isempty(excluded_artifact_grades)
    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades)

    target_match_str = make_signal_stem("tricorr"; 
        excluded_artifact_grades=excluded_artifact_grades,
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        snippets_duration_s=snippets_duration_s,
        remaining_params...
    )
    jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    signals = jld_dict["contributions"]
    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    fig = plot_μ_and_σ_signals_and_roc(results_df, signals, signal_times, seizure_bounds; analysis_eeg=eeg, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)")

    save(joinpath(save_dir, "$(task_name)_roc_patient$(patient_num)_reviewers$(min_reviewers_per_seizure).png"), fig)
end

let signal_type = "tricorr", reduction_type = "mostsignificant",
    patients_considered = [1:15..., 19,31,44,47,50,62];

params = Dict(
    :excluded_artifact_grades=>excluded_artifact_grades,
    :snippets_duration_s => snippets_duration_s,
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.zscore!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25), :patient_num => "",
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[1],
    :alert_grace_s => 60,
    :rolling_window_s => 60,
    :snippets_duration_s => 1
)

epochdiff_stem = make_epochdiff_stem(signal_type; params...)
maybe_dict = load_most_recent_jld2(epochdiff_stem, datadir())
results_df = if isnothing(maybe_dict)
    @error "Cannot find results_df like $(epochdiff_stem)"
else
    maybe_dict["results_df"]
end

detect_stem = make_detection_stem(signal_type, reduction_type; params...)
save_dir = plotsdir("$(detect_stem)$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, patients_considered) do patient_num
    patient_results = filter(:patient => p -> p == patient_num, results_df)
    detect_patient_seizures(patient_num; save_dir=save_dir, task_name="$(signal_type)$(reduction_type)", params...)
end

end