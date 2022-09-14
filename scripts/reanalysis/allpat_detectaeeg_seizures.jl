using DrWatson
@quickactivate "NeonateTriCorr"

## Saves contributions timeseries to datadir()/exp_pro/
## Excludes seconds marked as artifactual by

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

patients_all = 1:79
patients_artifact_annotated = [1:15..., 19,31,44,47,50,62,75]
patients_unannotated = setdiff(patients_all, patients_artifact_annotated) 

let signal_type = "aEEG", signals_reduction_name = "meanall",
    patients_considered = [1:15..., 19,31,44,47,50,62,75];

params = Dict(
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[1],
    :min_dist_to_seizure => 30,
    :alert_grace_s => 60,
    :rolling_window_s => 60,
    :signals_reduction_params => Dict{Symbol,Any}(
        :n_signals_used => 5,
        :rolling_window_s => 60,
        :signal_sym => :Δμ,
        :window_fn => mean
    ),
    :lowpass_freq => 0.31,
    :snippets_duration_s => 15,
    :lower_margin_perc => 0.09,
    :upper_margin_perc => 0.93,
    :min_snippets_for_comparison => 15,
    :n_θs => 100
)

# epochdiff_stem = make_epochdiff_stem(signal_type; params...)
# maybe_dict = load_most_recent_jld2(epochdiff_stem, datadir())
# results_df = if isnothing(maybe_dict)
#     @error "Cannot find results_df like $(epochdiff_stem)"
# else
#     maybe_dict["results_df"]
# end

detect_stem = make_detection_stem(signal_type, signals_reduction_name; params...)
session_id = Dates.now()
save_dir = plotsdir("$(detect_stem)$(session_id)")
mkpath(save_dir)

detect_all_patients_seizures(patients_considered; signal_type=signal_type, save_dir=save_dir, session_id=session_id, params..., signals_reduction_name=signals_reduction_name)
end
