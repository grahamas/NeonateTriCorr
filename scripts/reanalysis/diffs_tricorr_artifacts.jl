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

patients_all = 1:79
patients_artifact_annotated = [1:15..., 19,31,44,47,50,62]
patients_unannotated = setdiff(patients_all, patients_artifact_annotated) 

let patients = 76:79,
    signal_type = "tricorr";

params = Dict(
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.identity!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25),
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[],
    :min_dist_to_seizure => 30,
    :epoch_s => 60,
    :rolling_window_s => 60,
    :snippets_duration_s => 1,
    :min_snippets_for_comparison => 15
)


(results_df, fig_sig, fig_Δμ, fig_Δσ) = calculate_epoch_differences_across_patients(signal_type, patients; get_channel_label=offset_motif_numeral, params...)

end