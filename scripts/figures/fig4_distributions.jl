using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie, AlgebraOfGraphics
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, linecolor=:black, linewidth=4)
set_theme!(font_theme)

include(scriptsdir("include_src.jl"))

both_results = let signal_types = ["aEEG", "tricorr"],
    patients = patients_all;

map(signal_types) do st
params = merge(common_params, analysis_particular_params[st])

get_channel_label = if st == "tricorr"
    offset_motif_numeral
elseif st == "aEEG"
    nothing
end

(results, fig_sig, fig_mu, fig_sigma) = plot_epoch_differences_across_patients(st, patients; get_channel_label=get_channel_label, params...)

save_dir = plotsdir("$(st)_significance_results_$(Dates.time())")
mkpath(save_dir)
save(joinpath(save_dir, "$(st)_significance_facets.png"), fig_sig)
save(joinpath(save_dir, "$(st)_mean_change_facets.png"), fig_mu)
save(joinpath(save_dir, "$(st)_std_change_facets.png"), fig_sigma)

results

end
end