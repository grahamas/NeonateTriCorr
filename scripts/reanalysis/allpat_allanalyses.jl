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

let all_analyses_path = plotsdir("all_analyses_IntraStd_$(Dates.now())");
mkpath(all_analyses_path)


for signal_type ∈ ["tricorr", "aEEG"], 
    signals_reduction_name ∈ ["maxany"], 
    excluded_artifact_grades ∈ [Int[]], 
    min_reviewers_per_seizure ∈ 3:-1:1

    if signal_type == "tricorr"
        signals_reduction_name = "$(signals_reduction_name)abs"
    end

    patient_sets = if isempty(excluded_artifact_grades)
        [patients_all, patients_artifact_annotated]
    else
        [patients_artifact_annotated]
    end

    params = merge(common_params, analysis_particular_params[signal_type], 
        Dict(
            :excluded_artifact_grades => excluded_artifact_grades,
            :min_reviewers_per_seizure => min_reviewers_per_seizure
        )
    )

    detect_stem = make_detection_stem(signal_type, signals_reduction_name;
        params...
    )
    save_dir = joinpath(all_analyses_path, "$(detect_stem)$(Dates.now())")
    mkpath(save_dir)

    for patients_considered ∈ patient_sets
        detect_all_patients_seizures(patients_considered; 
            signal_type=signal_type, save_dir=save_dir, params...,
            signals_reduction_name=signals_reduction_name
        )
    end
end

end