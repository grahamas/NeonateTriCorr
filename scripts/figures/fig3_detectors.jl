using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, xticklabelsize=20, yticklabelsize=20, linecolor=:black, linewidth=4)
set_theme!(font_theme)


using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity
using ROCAnalysis

include(scriptsdir("include_src.jl"))


let signals_reduction_name = "maxany",
    patients_considered = patients_all,
    resolution=(1600,1200);

sr_dict = Dict(
    "aEEG" => signals_reduction_name,
    "tricorr" => "$(signals_reduction_name)abs"
)
colors = Dict(
    "tricorr" => :blue,
    "aEEG" => :red
)

save_dir = plotsdir("fig3_detectors_$(Dates.now())")
mkpath(save_dir)

epoch_s = common_params[:epoch_s]

fig = Figure(resolution=resolution, figure_padding=(10,35,10,35))
axes = Array{Any}(undef, 2, 2)

nts = []

for (i_pos, (pos_name, pos_fn)) ∈ enumerate([("seizure", posseizure_negepoch), ("patient", pospatient_negepoch)])
    for (i_std, standardization) ∈ enumerate(["within", "across"])
        fig[i_pos, i_std] = axes[i_pos, i_std] = ax = Axis(fig)
        for signal_type ∈ ["tricorr", "aEEG"]
            sr = sr_dict[signal_type]
            r, p = load_or_calculate_multipatient_ROC(
                patients_considered, signal_type, sr
                ; 
                common_params..., 
                analysis_particular_params[signal_type]..., 
                standardization=standardization, 
                calculate_targets_fn=pos_fn
            )
            plot_roc_fphr!(ax, r; color=colors[signal_type], label=signal_type)
            push!(nts, (target=pos_name, standardization=fn2str(standardization), signal=signal_type, auc=auc(r), ))
        end
    end
end
axes[1,1].ylabel = "per-seizure TPR"
axes[2,1].ylabel = "per-patient TPR"
axes[2,1].xlabel = "FP/Hr"
axes[2,2].xlabel = "FP/Hr"

axislegend(axes[2,2]; halign=:right, valign=:bottom)

hidedecorations!.(axes; label = false, ticklabels = false, ticks = false, grid = false)
hidexdecorations!.(axes[1,1:2]; label=true, ticklabels=true, grid=false)
hideydecorations!.(axes[1:2,2]; grid=false)

colgap!(fig.layout,50)
rowgap!(fig.layout,50)

save(joinpath(save_dir, "compare_aEEG_tricorr_standardizations_patients$(length(patients_considered))_$(signals_reduction_name)_ROCs_$(Dates.now()).$(ext)"), fig)

return fig

end