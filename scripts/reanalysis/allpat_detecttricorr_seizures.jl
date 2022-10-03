using DrWatson
@quickactivate "NeonateTriCorr"

## Saves contributions timeseries to datadir()/exp_pro/
## Excludes seconds marked as artifactual by

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity
using ROCAnalysis

include(scriptsdir("include_src.jl"))

noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=24, Lines=Theme(color=:black), linewidth=4)
set_theme!(font_theme)

tricorr_targets_nontargets = let signal_type = "tricorr", signals_reduction_name = "maxany",
    patients_considered = patients_all;
 
params = merge(common_params, analysis_particular_params["tricorr"])

# epochdiff_stem = make_epochdiff_stem(signal_type; params...)
# maybe_dict = load_most_recent_jld2(epocwhdiff_stem, datadir())
# results_df = if isnothing(maybe_dict)
#     @error "Cannot find results_df like $(epochdiff_stem)"
# else
#     maybe_dict["results_df"]
# end

detect_stem = make_detection_stem(signal_type, signals_reduction_name; params...)
save_dir = plotsdir("$(detect_stem)$(Dates.now())")
mkpath(save_dir)

detect_all_patients_seizures(patients_considered; signal_type=signal_type, save_dir=save_dir, params..., signals_reduction_name=signals_reduction_name)
end

# tricorr_results = mapreduce(add_nts, zip(tricorr_sig_times_bounds...)) do (signal, times, bounds)
#     evaluate_detection_posseizure_negepoch(bounds, signal .> 1.0, times; snippets_duration_s=1, epoch_s=60)
# end