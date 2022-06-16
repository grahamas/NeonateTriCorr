using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

## Saves contributions timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using HypothesisTests, CSV, Distances, LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

motif_results = mapreduce(vcat, [1:15..., 19,31,44,47,50,62]) do PAT

eeg = load_helsinki_eeg(PAT)

min_snippets_for_comparison = 150

min_dist_to_seizure = 30
snippets_duration = 1
target_match_str = "tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(snippets_duration)_lagextents8x25_helsinkiEEG$(PAT)"
save_dir = plotsdir("$(target_match_str)_$(Dates.now())")
mkpath(save_dir)

seizure_bounds, consensus = load_helsinki_seizure_annotations(PAT)

if isempty(seizure_bounds)
    @warn "Patient $(PAT): no seizures."
    return []
end

jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
contributions = jld_dict["contributions"]

control_snippets = calc_control_snippet_starts_incl_artifacts(seizure_bounds, size(contributions, :time), snippets_duration, min_dist_to_seizure)
seizure_snippets = mapreduce(vcat, seizure_bounds) do (on, off)
    on:snippets_duration:off
end

if (length(control_snippets) < min_snippets_for_comparison)
    @warn "Not enough control snippets ($(length(control_snippets))) in Patient $PAT"
    return []
end
if (length(seizure_snippets) < min_snippets_for_comparison) 
    @warn "Not enough seizure snippets ($(length(seizure_snippets))) in Patient $PAT"
    return []
end


function estimate_JSDistance(control, seizure; show_plots=false)
    combined_kde = kde(vcat(control, seizure))
    common_grid = combined_kde.x
    control_kde = kde(control, common_grid)
    seizure_kde = kde(seizure, common_grid)

    if show_plots
        plt = plot(control_kde.x, control_kde.density)
        plot!(plt.axis, seizure_kde.x, seizure_kde.density, color=:red)
        display(plt)
    end

    control_dist = control_kde.density ./ sum(control_kde.density)
    seizure_dist = seizure_kde.density ./ sum(seizure_kde.density)

    sqrt(evaluate(JSDivergence(), control_dist, seizure_dist))
end

# THESE SHOULD NOT
control_contributions = contributions[:, control_snippets]
seizure_contributions = contributions[:, seizure_snippets]

motif_JSDistances_KSpvalues_effect = map(1:14) do motif_num
    motif_control_contributions = control_contributions[motif_num, :] |> skipmissing |> collect
    motif_seizure_contributions = seizure_contributions[motif_num, :] |> skipmissing |> collect

    (JS=estimate_JSDistance(motif_control_contributions, motif_seizure_contributions; show_plots=false),
     p=pvalue(ApproximateTwoSampleKSTest(motif_control_contributions, motif_seizure_contributions)),
     Δμ=mean(motif_seizure_contributions) - mean(motif_control_contributions),
     Δσ=std(motif_seizure_contributions) - std(motif_control_contributions),
     patient=PAT,
     motif=offset_motif_numeral(motif_num))
end

function plot_estimated_distributions(control_contributions, seizure_contributions; title)

    df = DataFrame(mapreduce(vcat, 1:14) do motif_num
        motif_control_contributions = control_contributions[motif_num, :] |> skipmissing |> collect
        motif_seizure_contributions = seizure_contributions[motif_num, :] |> skipmissing |> collect

        motif_contributions = vcat(motif_control_contributions, motif_seizure_contributions)
        motif_kde = kde(motif_contributions)
        common_grid = motif_kde.x
        typeof(motif_kde.x)

        motif_control_kde = kde(motif_control_contributions, common_grid)
        motif_seizure_kde = kde(motif_seizure_contributions, common_grid)

        seizure_nt = [(x=x, density=density, motif=offset_motif_numeral(motif_num), condition="seizure") for (x, density) in zip(motif_seizure_kde.x, motif_seizure_kde.density)]
        control_nt = [(x=x, density=density, motif=offset_motif_numeral(motif_num), condition="control") for (x, density) in zip(motif_control_kde.x, motif_control_kde.density)]
        vcat(seizure_nt, control_nt)
    end)

    plt = data(df) * mapping(:x, :density, color=:condition, layout=:motif => sorter(1:14 .|> offset_motif_numeral))
    fig = draw(plt; facet=(; linkxaxes=:none, linkyaxes=:none), figure=(;resolution = (2400, 2400), fontsize=38))
    Label(fig.figure[0,:], title, textsize=56)
    fig
end

fig_kdes = plot_estimated_distributions(control_contributions, seizure_contributions; title="Patient = $(PAT)")
save(joinpath(save_dir, "estimated_distributions_pat$(PAT).png"), fig_kdes)

rolling_window = 30
rolling_std = mapreduce(hcat, 1:14) do motif_num
    [std(contributions[motif_num,i:(i+rolling_window)]) for i ∈ 1:(size(contributions,:time)-rolling_window)]
end
rolling_mean = mapreduce(hcat, 1:14) do motif_num
    [mean(contributions[motif_num,i:(i+rolling_window)]) for i ∈ 1:(size(contributions,:time)-rolling_window)]
end

fig_std = plot_contributions(eeg, get_times(eeg, sample_rate=1)[1:end-rolling_window], rolling_std'; title="Rolling σ (window = $(rolling_window))", resolution=(800,1400))
save(joinpath(save_dir, "contributions_standard_deviations_pat$(PAT).png"), fig_std)

fig_mean = plot_contributions(eeg, get_times(eeg, sample_rate=1)[1:end-rolling_window], rolling_mean'; title="Rolling μ (window = $(rolling_window))", resolution=(800,1400))
save(joinpath(save_dir, "contributions_means_pat$(PAT).png"), fig_mean)

motif_JSDistances_KSpvalues_effect

end

function significance_category(x)
    if x > 0.05 / 14
        return "NS"
    elseif x > 0.01 / 14
        return "* 0.05/14"
    elseif x > 0.001 / 14
        return "** 0.01/14"
    elseif x > 0.0001 / 14
        return "*** 0.001/14"
    else
        return "**** 0.0001/14"
    end
end

complete_df = DataFrame(motif_results)

summary = combine(groupby(complete_df, :motif), :Δμ => std, :Δσ => std, :Δμ => mean, :Δσ => mean)
function normed_Δμ(xs, motifs)
    map(zip(xs, motifs)) do (x, motif)
        σ = only(filter(:motif => ==(motif), summary).Δμ_std)
        return x / σ
    end
end
function normed_Δσ(xs, motifs)
    map(zip(xs, motifs)) do (x, motif)
        σ = only(filter(:motif => ==(motif), summary).Δσ_std)
        return x / σ
    end
end

transform!(complete_df, [:Δμ, :motif] => normed_Δμ => :Δμ_z, [:Δσ, :motif] => normed_Δσ => :Δσ_z)

complete_df.significance = significance_category.(complete_df.p)
sig_categories = ([0.1, 0.04, 0.02, 0.002, 0.0002, 0.00002] ./ 14) .|> significance_category
colors = cgrad(:sunset, length(sig_categories), categorical=true)
plt = data(complete_df) * mapping(:patient => nonnumeric, :motif => sorter(1:14 .|> offset_motif_numeral), color=:significance => sorter(sig_categories) => "significance") * visual(markersize=30); fig_significance = draw(plt; palettes=(color=colors,))


max_Δμ = maximum(abs.(extrema(complete_df.Δμ_z)))
max_Δσ = maximum(abs.(extrema(complete_df.Δσ_z)))

plt = data(complete_df) * mapping(:patient => nonnumeric, :motif => sorter(1:14 .|> offset_motif_numeral),  color=:Δμ_z) * visual(markersize=30, colormap=:bluesreds, colorrange=(-max_Δμ,max_Δμ)); fig_Δμ = draw(plt)

plt = data(complete_df) * mapping(:patient => nonnumeric, :motif => sorter(1:14 .|> offset_motif_numeral),  color=:Δσ_z) * visual(markersize=30, colormap=:bluesreds, colorrange=(-max_Δσ,max_Δσ)); fig_Δσ = draw(plt)