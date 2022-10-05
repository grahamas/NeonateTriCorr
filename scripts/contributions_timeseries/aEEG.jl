force_recalculate_aEEG = false

using DrWatson
@quickactivate "NeonateTriCorr"

## Saves aEEG timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV
using ROCAnalysis

include(scriptsdir("include_src.jl"))

params = merge(common_params, analysis_particular_params["aEEG"])
for PAT ∈ 1:79
aEEG_PAT = calculate_patient_aEEG(PAT; 
    params..., force_recalculate_aEEG=force_recalculate_aEEG,
    plot_traces = true
)
end

force_recalculate_aEEG = false