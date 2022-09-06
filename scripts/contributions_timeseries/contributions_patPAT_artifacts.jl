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

contributions_PAT = calculate_patient_tricorr(PAT; excluded_artifact_grades=Int[])

force_recalculate_contributions = false