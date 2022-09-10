using DrWatson
@quickactivate "NeonateTriCorr"

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using HypothesisTests, CSV, Distances, LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

let signal_type = "tricorr", reduction_type = "distance",
    patients_considered = [1:15..., 19,31,44,47,50,62];

params = Dict(
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.zscore!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25), :patient_num => "",
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

# download recordings if not present
download_helsinki_eegs(PATs)

# calculate triple correlations
for PAT in PATs
    eeg = load_helsinki_eeg(PAT; excluded_artifact_grades=params[:excluded_artifact_grades])
    contributions = calculate_patient_tricorr(PAT; eeg=eeg, params...)
    
