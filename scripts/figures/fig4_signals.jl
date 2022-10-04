using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, xticklabelsize=20, yticklabelsize=20, linecolor=:black, linewidth=4, labelsize=36)
set_theme!(font_theme)


using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity
using ROCAnalysis

include(scriptsdir("include_src.jl"))

let signals_reduction_name = "maxany",
    patients_considered = [9,12,13,62],
    resolution=(1600,1200);

signal_types = ["tricorr", "aEEG"]

sr_dict = Dict(
    "aEEG" => signals_reduction_name,
    "tricorr" => "$(signals_reduction_name)abs"
)
colors = Dict(
    "tricorr" => :blue,
    "aEEG" => :red
)

save_dir = plotsdir("fig4_signals_$(Dates.now())")
mkpath(save_dir)

epoch_s = common_params[:epoch_s]

signals = Dict(
    pat_num => Dict(
        sig => load_signal(; common_params..., analysis_particular_params[sig]..., patient_num=pat_num) for sig in signal_types
    ) for pat_num in patients_considered
)

fig = Figure(resolution=resolution)

axes = Array{Any}(undef, 2, 3, 2)

map(eachindex(axes)) do idx
    axes[idx] = Axis(fig)
end

patients_layout = [
    9 13 12;
    9 62 12
]

for idx ∈ CartesianIndices(patients_layout)
    confusion_entry = GridLayout()
for (i_st, st) ∈ enumerate(signal_types)
    signal = signals[patients_layout[idx]][st]
    dt = analysis_particular_params[st][:snippets_duration_s]
    time = 0.:dt:((length(signal)-1)*dt)
    lines!(axes[idx, i_st], time, signal)
    confusion_entry[i_st] = axes[idx, i_st]
end
fig[idx] = confusion_entry
end

fig

end