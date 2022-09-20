force_recalculate_aEEG = false

using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

## Saves aEEG timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"
using Random, JLD2
using HypothesisTests, CSV

include(scriptsdir("include_src.jl"))

params = Dict(
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[1],
    :min_dist_to_seizure => 30,
    :alert_grace_s => 60,
    :rolling_window_s => 60,
    :rolling_window_s => 60,
    :window_fn => mean,
    :lowpass_freq => 0.31,
    :snippets_duration_s => 15,
    :lower_margin_perc => 0.09,
    :upper_margin_perc => 0.93,
    :min_snippets_for_comparison => 15,
    :n_θs => 100
)
for PAT ∈ 1:79
aEEG_PAT = calculate_patient_aEEG(PAT; 
    params..., force_recalculate_aEEG=force_recalculate_aEEG,
    plot_traces = true
)
end

force_recalculate_aEEG = false