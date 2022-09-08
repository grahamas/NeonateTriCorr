force_recalculate_contributions = false

using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

## Saves contributions timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV

include(scriptsdir("include_src.jl"))

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
    :snippets_duration_s => 15,
    :signals_reduction_params => Dict{Symbol,Any}(
        :n_signals_used => 5    
    )
)
contributions_PAT = calculate_patient_tricorr(PAT; 
    params...
)

force_recalculate_contributions = false