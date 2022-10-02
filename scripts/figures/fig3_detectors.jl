using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, linecolor=:black, linewidth=4)
set_theme!(font_theme)


using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

eval_dct = Dict(
    "pospatient_negepoch" => (evaluate_detection_pospatient_negepoch, "patient"),
    "posseizure_negepoch" => (evaluate_detection_posseizure_negepoch, "seizure")
)


let aEEG_signals_reduction_name = "maxany",
    tricorr_signals_reduction_name = "$(aEEG_signals_reduction_name)abs",
    patients_considered = patients_all,
    resolution=(1200,800);



eval_str = "posseizure_negepoch"
eval_fn, eval_tpr = eval_dct[eval_str]

save_dir = plotsdir("fig3_detectors_$(eval_str)_$(Dates.now())")
mkpath(save_dir)

modified_common_params = merge(common_params, Dict(
    :evaluation_fn => eval_fn
))

aEEG_params = merge(modified_common_params, analysis_particular_params["aEEG"])
tricorr_params = merge(modified_common_params, analysis_particular_params["tricorr"])

epoch_s = aEEG_params[:epoch_s]
@assert tricorr_params[:epoch_s] == epoch_s

tricorr_within_ROC_df, loaded_tricorr_params = load_or_calculate_multipatient_ROC_df(patients_considered, "tricorr", tricorr_signals_reduction_name; tricorr_params..., standardization="within")
aEEG_within_ROC_df, loaded_aEEG_params = load_or_calculate_multipatient_ROC_df(patients_considered, "aEEG", aEEG_signals_reduction_name; aEEG_params..., standardization="within")

tricorr_across_ROC_df, loaded_tricorr_params = load_or_calculate_multipatient_ROC_df(patients_considered, "tricorr", tricorr_signals_reduction_name; tricorr_params..., standardization="across")
aEEG_across_ROC_df, loaded_aEEG_params = load_or_calculate_multipatient_ROC_df(patients_considered, "aEEG", aEEG_signals_reduction_name; aEEG_params..., standardization="across")

tricorr_color = :blue
aEEG_color = :red

tricorr_across_ROC_df[!, :signal] .= :tricorr
tricorr_across_ROC_df[!, :standardization] .= :across

tricorr_within_ROC_df[!, :signal] .= :tricorr
tricorr_within_ROC_df[!, :standardization] .= :within

aEEG_across_ROC_df[!, :signal] .= :aEEG
aEEG_across_ROC_df[!, :standardization] .= :across

aEEG_within_ROC_df[!, :signal] .= :aEEG
aEEG_within_ROC_df[!, :standardization] .= :within


all_ROC_df = vcat(tricorr_across_ROC_df, tricorr_within_ROC_df, aEEG_across_ROC_df, aEEG_within_ROC_df)

plt = data(all_ROC_df) * mapping(
    (:false_positives,:gt_negative) => 
            ((f, gt) -> f / (gt * epoch_s / (60 * 60))) => 
            "FP/hour", 
    (:true_positives,:gt_positive) => 
            ((t, gt) -> t / gt) => 
            "per-$(eval_tpr) TPR",
    color=:signal, col=:standardization
) * visual(Lines, linewidth=5)

fig = Figure(resolution=resolution)
drw = draw!(fig, plt; 
    axis=(limits=((0.,60.),(0.,1.)),)
)
leg = AlgebraOfGraphics.compute_legend(drw)
right_ax = content(fig[1,2])
hidedecorations!(right_ax; label = false, ticklabels = false, ticks = false, grid = true)
hidedecorations!(content(fig[1,1]); label = false, ticklabels = false, ticks = false, grid = true)
axislegend(right_ax, leg...; halign=:right, valign=:bottom)

save(joinpath(save_dir, "compare_aEEG_tricorr_standardizations_patients$(length(patients_considered))_$(aEEG_signals_reduction_name)_ROCs_$(Dates.now()).$(ext)"), fig)

return fig

end