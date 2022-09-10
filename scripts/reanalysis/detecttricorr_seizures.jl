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
        task_name, remaining_params...
    )
    @assert !isempty(excluded_artifact_grades)
    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades)

    target_match_str = make_signal_stem("tricorr"; 
        excluded_artifact_grades=excluded_artifact_grades,
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        snippets_duration_s=snippets_duration_s,
        remaining_params...
    )
    maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    signals = if isnothing(maybe_dict)
        @error "No data like $target_match_str"
    else
        maybe_dict["contributions"]
    end
    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    fig = plot_μ_and_σ_signals_and_roc(signals, signal_times, seizure_bounds; analysis_eeg=eeg, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)", remaining_params...)

    save(joinpath(save_dir, "$(task_name)_roc_patient$(patient_num)_reviewers$(min_reviewers_per_seizure).png"), fig)
end

let signal_type = "tricorr", reduction_type = "meansignificant",
    patients_considered = [1:15..., 19,31,44,47,50,62];

params = Dict(
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.zscore!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25), 
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[1],
    :min_dist_to_seizure => 30,
    :alert_grace_s => 60,
    :rolling_window_s => 60,
    :snippets_duration_s => 1,
    :signals_reduction_params => Dict{Symbol,Any}(
        :n_signals_used => 5    
    )
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
    params[:signals_reduction_params][:results_df] = patient_results
    detect_patient_seizures(patient_num; save_dir=save_dir, task_name="$(signal_type)$(reduction_type)", params..., signals_reduction_name=reduction_type, patient_num=patient_num)
end

end