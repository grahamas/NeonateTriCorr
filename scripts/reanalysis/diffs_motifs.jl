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
include(scriptsdir("meats", "epoch_differences.jl"))

motif_results = mapreduce(vcat, [1]) do PAT#[1:15..., 19,31,44,47,50,62]) do PAT

params = Dict(
    :excluded_artifact_grades=>Int[],
    :snippets_duration_s => 1,
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.zscore!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25), :patient_num => PAT
)
target_match_str = make_filename_stem("tricorr"; params...)

calculate_epoch_differences(target_match_str; params...)

end