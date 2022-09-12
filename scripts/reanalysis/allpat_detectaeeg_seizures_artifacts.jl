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

let signal_type = "tricorr", signals_reduction_name = "meanall",
    patients_considered = 1:79;#[1:15..., 19,31,44,47,50,62];

params = Dict(
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.identity!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :lag_extents => (8,25),
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[],
    :min_dist_to_seizure => 30,
    :alert_grace_s => 60,
    :snippets_duration_s => 1,
    :signals_reduction_params => Dict{Symbol,Any}(
        :n_signals_used => 5,
        :rolling_window_s => 60,
        :signal_sym => :Δμ,
        :window_fn => mean
    ),
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
save_dir = plotsdir("$(detect_stem)$(Dates.now())")
mkpath(save_dir)

detect_all_patients_seizures(patients_considered; save_dir=save_dir, task_name="$(signal_type)$(signals_reduction_name)", params..., signals_reduction_name=signals_reduction_name, signal_from_dct_fn = (dct -> dct["contributions"]))
end
