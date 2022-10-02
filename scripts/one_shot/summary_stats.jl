
using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, linecolor=:black, linewidth=4)
set_theme!(font_theme)


using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

function summary_stats(patients::AbstractArray{<:Number}; kwargs...)
    eegs = load_helsinki_eeg.(patients; kwargs...);
    (summary_stats(eegs), eegs)
end

function summary_stats(eegs::AbstractArray{<:AbstractEEG})
    
    (
        n_patients = length(eegs),
        n_artifact_patients = count(.!isempty.(getproperty.(eegs, :artifact_annotations))),
        n_artifacts = sum(length.(getproperty.(eegs, :artifact_annotations))),
        std_n_artifacts = std(length.(getproperty.(eegs, :artifact_annotations))),
        artifacts_duration = sum(sum.(map.(((x,y),) -> y - x, getproperty.(eegs, :artifact_annotations)))),
        std_artifacts_duration = std(sum.(map.(((x,y),) -> y - x, getproperty.(eegs, :artifact_annotations)))),
        n_seizing_patients = count(.!isempty.(getproperty.(eegs, :seizure_annotations))),
        n_seizures = sum(length.(getproperty.(eegs, :seizure_annotations))),
        std_n_seizures = std(length.(filter(!isempty, getproperty.(eegs, :seizure_annotations)))),
        seizures_duration = sum(sum.(map.(((x,y),) -> y - x, getproperty.(eegs, :seizure_annotations)))),
        std_seizures_duration = std(sum.(map.(((x,y),) -> y - x, filter(!isempty, getproperty.(eegs, :seizure_annotations))))),
        total_duration = sum(getproperty.(eegs, :duration)),
        std_total_duration = std(getproperty.(eegs, :duration))
    )
end

summary_all, all_eegs = summary_stats(patients_all; excluded_artifact_grades=Int[], min_reviewers_per_seizure=3, discretization_s=1)

summary_no_artifacts, no_artifacts_eegs = summary_stats(patients_artifact_annotated; excluded_artifact_grades=Int[1], min_reviewers_per_seizure=3, discretization_s=1)

summary_annotated, annotated_eegs = summary_stats([patients_artifact_annotated..., 21]; excluded_artifact_grades=Int[], min_reviewers_per_seizure=3, discretization_s=1)

