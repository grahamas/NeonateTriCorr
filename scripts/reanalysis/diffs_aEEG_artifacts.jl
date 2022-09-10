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

let patients = patients_all,
    signal_type = "aEEG";


    params = Dict(
        
        :min_reviewers_per_seizure => 1,
        :excluded_artifact_grades => Int[],
        :min_dist_to_seizure => 30,
        :alert_grace_s => 60,
        :rolling_window_s => 60,
        :signals_reduction_params => Dict{Symbol,Any}(
            :n_signals_used => 5    
        ),
        :lowpass_freq => 0.31,
        :snippets_duration_s => 15,
        :lower_margin_perc => 0.09,
        :upper_margin_perc => 0.93,
        :min_snippets_for_comparison => 15
    )


(results_df, fig_sig, fig_Δμ, fig_Δσ) = calculate_epoch_differences_across_patients(signal_type, patients; params...)

end