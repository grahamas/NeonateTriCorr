using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, linecolor=:black, linewidth=4)
set_theme!(font_theme)

include(scriptsdir("plot_together/standardization_comparison.jl"))