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

let signal_type = "aEEG", reduction_type = "meanall",
    patients_considered = 1:75;#[1:15..., 19,31,44,47,50,62];

    params = Dict(
        
        :min_reviewers_per_seizure => 3,
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

epochdiff_stem = make_epochdiff_stem(signal_type; params...)
maybe_dict = load_most_recent_jld2(epochdiff_stem, datadir())
results_df = if isnothing(maybe_dict)
    @error "Cannot find results_df like $(epochdiff_stem)"
else
    maybe_dict["results_df"]
end

detect_stem = make_detection_stem(signal_type, reduction_type; params...)
save_dir = plotsdir("$(detect_stem)$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, patients_considered) do patient_num
    patient_results = filter(:patient => p -> p == patient_num, results_df)
    params[:signals_reduction_params][:results_df] = patient_results
    detect_patient_seizures(patient_num; save_dir=save_dir, task_name="$(signal_type)$(reduction_type)", params..., signals_reduction_name=reduction_type, patient_num=patient_num, signals_from_dct_fn = (dct -> aEEG_lower_margin(dct["aEEG"])))
end

end