using DrWatson
@quickactivate "NeonateTriCorr"

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using HypothesisTests, CSV, Distances, LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

let PATs = [1:15..., 19,31,44,47,50,62];
    # all pats are 1:75
    # artifact annotated pats are [1:15..., 19,31,44,47,50,62]

params = Dict(
    :excluded_artifact_grades=>Int[],
    :snippets_duration_s => 1,
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.zscore!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25), :patient_num => PAT
)

# download recordings if not present
download_helsinki_eegs(PATs)

# calculate triple correlations
for PAT in PATs
    eeg = load_helsinki_eeg(PAT; excluded_artifact_grades=params[:excluded_artifact_grades])
    contributions = calculate_patient_tricorr(PAT; eeg=eeg, params...)
    
