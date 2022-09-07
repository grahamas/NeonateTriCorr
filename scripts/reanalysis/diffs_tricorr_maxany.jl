using DrWatson
@quickactivate "NeonateTriCorr"

# Assumes you have previously run contributions_patPAT.jl
#   and loads most recently calc'd contributions from datadir("exp_pro")

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using HypothesisTests, CSV, Distances, LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

patients_all = 1:75
patients_artifact_annotated = [1:15..., 19,31,44,47,50,62]
patients_unannotated = setdiff(patients_all, patients_artifact_annotated) 

let patients = patients_artifact_annotated,
    signal_type = "tricorr";

params = Dict(
    :excluded_artifact_grades=>Int[1],
    :snippets_duration_s => 1,
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.zscore!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25), :patient_num => "",
    :min_reviewers_per_seizure => 3,
    :min_dist_to_seizure => 30,
    :min_snippets_for_comparison => 80
)

(results_df, fig_sig, fig_Δμ, fig_Δσ) = calculate_epoch_differences_across_patients(signal_type, patients; params...)

end