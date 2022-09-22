using DrWatson
@quickactivate "NeonateTriCorr"

## Saves contributions timeseries to datadir()/exp_pro/
## Excludes seconds marked as artifactual by

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "svg"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

patients_all = 1:79
patients_artifact_annotated = [1:15..., 19,31,44,47,50,62,75]
patients_unannotated = setdiff(patients_all, patients_artifact_annotated) 

let aEEG_signals_reduction_name = "maxany",
    tricorr_signals_reduction_name = "$(aEEG_signals_reduction_name)abs",
    patients_considered = patients_all,
    standardization="within",
    resolution=(800,600);

aEEG_params = merge(common_params, analysis_particular_params["aEEG"], Dict(:standardization => standardization))
tricorr_params = merge(common_params, analysis_particular_params["tricorr"], Dict(:standardization => standardization))

tricorr_ROC_df, loaded_tricorr_params = load_or_calculate_multipatient_ROC_df(patients_considered, "tricorr", tricorr_signals_reduction_name; tricorr_params...)
aEEG_ROC_df, loaded_aEEG_params = load_or_calculate_multipatient_ROC_df(patients_considered, "aEEG", aEEG_signals_reduction_name; aEEG_params...)

fig = Figure(resolution=resolution)

tricorr_color = :blue
aEEG_color = :red
combined_plt = plot_NODRAW_seizure_detection_ROC!(tricorr_ROC_df, color=tricorr_color, epoch_s=aEEG_params[:epoch_s]) + plot_NODRAW_seizure_detection_ROC!(aEEG_ROC_df, color=aEEG_color, epoch_s=tricorr_params[:epoch_s])

plot_seizure_detection_ROC!(fig[1,1], DataFrame(); title="Tricorr ($(string(tricorr_color))) vs aEEG ($(string(aEEG_color))) ROC $(standardization)-patient standardized", roc_plt = combined_plt)
save(plotsdir("combined_aEEG_tricorr_$(standardization)_patients$(length(patients_considered))_$(aEEG_signals_reduction_name)_ROCs_$(Dates.now()).$(ext)"), fig)

return fig

end