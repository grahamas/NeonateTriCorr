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

let all_analyses_path = plotsdir("all_analyses_$(Dates.now())");
mkpath(all_analyses_path)


for signal_type ∈ ["tricorr", "aEEG"], 
    signals_reduction_name ∈ [raw"meanall", "maxany"], 
    artifacts_excluded ∈ [Int[], [1]], 
    reviewer_consensus ∈ 1:3

    if signal_type == "tricorr"
        signals_reduction_name = "$(signals_reduction_name)abs"
    end

    patient_sets = if isempty(artifacts_excluded)
        [patients_all, patients_artifact_annotated]
    else
        [patients_artifact_annotated]
    end

    detect_stem = make_detection_stem(signal_type, signals_reduction_name;
        params...
    )
    save_dir = joinpath(all_analyses_path, "$(detect_stem)$(Dates.now())")
    mkpath(save_dir)

    detect_all_patients_seizures(patients_considered; 
        signal_type=signal_type, save_dir=save_dir, params...,
        signals_reduction_name=signals_reduction_name
    )
end

