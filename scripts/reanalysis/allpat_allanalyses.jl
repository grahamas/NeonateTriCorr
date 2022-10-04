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

nts = []

for signal_type ∈ ["tricorr", "aEEG"], 
    signals_reduction_name ∈ ["maxany"],
    standardization ∈ ["within", "across"],
    (pos_name, pos_fn) ∈ [("seizure", posseizure_negepoch), ("patient", pospatient_negepoch)],
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
        r = detect_all_patients_seizures(patients_considered; 
            signal_type=signal_type, save_dir=save_dir, params...,
            signals_reduction_name=signals_reduction_name,
            standardization=standardization,
            calculate_targets_fn=pos_fn
        )
        push!(nts, (target=pos_name, standardization=standardization, signal=signal_type, reduction=signals_reduction_name, min_reviewers_per_seizure=min_reviewers_per_seizure, patient_set=patients_considered, excluded_artifact_grades=excluded_artifact_grades, auc=(1-auc(r)), tpr57=tpr_at_fphr(r,5.7), fphr80=fphr_at_tpr(r,0.8)))
    end
end

df = DataFrame(nts)

save(joinpath(all_analyses_path, "sensitivity_df.jld2"), Dict("sensitivity_df" => df))

CSV.write(joinpath(all_analyses_path, "sensitivity_df.csv"), df)

reduced_df = filter(:min_reviewers_per_seizure => ==(3), filter(:patient_set => ==("1:79"), df))
sort!(reduced_df, [:target, :standardization, :signal])

open(joinpath(all_analyses_path, "sensitivity_3rev_all_patients.tex"), "w") do io
    write(io, latexify_sensitivity(reduced_df))
end


en