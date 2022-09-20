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

let signals_reduction_name = "meanall",
    patients_considered = patients_all,
    resolution=(800,600);

base_params = Dict(
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[],
    :min_dist_to_seizure => 30,
    :alert_grace_s => 60,
    :rolling_window_s => 60,
    :rolling_window_s => 60,
    :window_fn => mean,
    :n_θs => 100
)

tricorr_params = Dict(base_params..., 
    :preproc! => TripleCorrelations.zscore!, 
    :postproc! => TripleCorrelations.identity!,
    :assumption => IndStdNormal(), :conditioned_on => None(),
    :snippets_duration_s => 1,
    :lag_extents => (8,25)
)
aEEG_params = Dict(base_params..., 
    :lowpass_freq => 0.31,
    :snippets_duration_s => 15,
    :lower_margin_perc => 0.09,
    :upper_margin_perc => 0.93,
    :min_snippets_for_comparison => 15
)

tricorr_ROC_df, loaded_tricorr_params = load_or_calculate_multipatient_ROC_df(patients_considered, "tricorr", signals_reduction_name; tricorr_params...)
aEEG_ROC_df, loaded_aEEG_params = load_or_calculate_multipatient_ROC_df(patients_considered, "aEEG", signals_reduction_name; aEEG_params...)

fig = Figure(resolution=resolution)

tricorr_color = :blue
aEEG_color = :red
combined_plt = plot_NODRAW_seizure_detection_ROC_standard!(tricorr_ROC_df, NamedTuple(loaded_tricorr_params).non_seizure_hours, color=tricorr_color) + plot_NODRAW_seizure_detection_ROC_standard!(aEEG_ROC_df, NamedTuple(loaded_aEEG_params).non_seizure_hours, color=aEEG_color)

plot_seizure_detection_ROC!(fig[1,1], DataFrame(); non_seizure_hours=NaN,title="Combined tricorr ($(string(tricorr_color))) and aEEG ($(string(aEEG_color))) ROC", roc_plt = combined_plt)

save(plotsdir("combined_aEEG_tricorr_patients$(length(patients_considered))_$(signals_reduction_name)_ROCs_STANDARD_$(Dates.now()).png"), fig)

return fig

end