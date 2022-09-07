
using TriCorrApplications
function TriCorrApplications.plot_reviewer_consensus!(ax, eeg::AbstractProcessedEEG)
    lines!(ax, get_times(eeg, sample_rate=1), eeg.seizure_reviewers_count)
    tightlimits!(ax); hidespines!(ax)
    ylims!(ax, (0, 3))
    hidedecorations!(ax, ticklabels=false)
end


function plot_estimated_distributions(control_contributions, seizure_contributions; title, get_label=offset_motif_numeral)

    df = DataFrame(mapreduce(vcat, 1:size(control_contributions,1)) do motif_num
        motif_control_contributions = control_contributions[motif_num, :] |> skipmissing |> collect
        motif_seizure_contributions = seizure_contributions[motif_num, :] |> skipmissing |> collect

        motif_contributions = vcat(motif_control_contributions, motif_seizure_contributions)
        motif_kde = kde(motif_contributions)
        common_grid = motif_kde.x
        typeof(motif_kde.x)

        motif_control_kde = kde(motif_control_contributions, common_grid)
        motif_seizure_kde = kde(motif_seizure_contributions, common_grid)

        seizure_nt = [(x=x, density=density, motif=get_label(motif_num), condition="seizure") for (x, density) in zip(motif_seizure_kde.x, motif_seizure_kde.density)]
        control_nt = [(x=x, density=density, motif=get_label(motif_num), condition="control") for (x, density) in zip(motif_control_kde.x, motif_control_kde.density)]
        vcat(seizure_nt, control_nt)
    end)

    plt = data(df) * mapping(:x, :density, color=:condition, layout=:motif => sorter(1:size(control_contributions,1) .|> get_label))
    fig = draw(plt; facet=(; linkxaxes=:none, linkyaxes=:none), figure=(;resolution = (2400, 2400), fontsize=38))
    Label(fig.figure[0,:], title, textsize=56)
    fig
end


function significance_category(x, n_comparisons)
    if x > 0.05 / n_comparisons
        return "NS"
    elseif x > 0.01 / n_comparisons
        return "* 0.05/$(n_comparisons)"
    elseif x > 0.001 / n_comparisons
        return "** 0.01/$(n_comparisons)"
    elseif x > 0.0001 / n_comparisons
        return "*** 0.001/$(n_comparisons)"
    else
        return "**** 0.0001/$(n_comparisons)"
    end
end

function draw_significances_plot!(df; all_motifs = 1:14 .|> offset_motif_numeral, draw_kwargs...)
    summary = combine(groupby(df, :motif), :Δμ => std, :Δσ => std, :Δμ => mean, :Δσ => mean)
    @show summary
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
    transform!(df, [:Δμ, :motif] => normed_Δμ => :Δμ_z, [:Δσ, :motif] => normed_Δσ => :Δσ_z)

    df.significance = significance_category.(df.p, length(all_motifs))
    sig_categories = ([0.1, 0.04, 0.02, 0.002, 0.0002, 0.00002] ./ length(all_motifs)) .|> Base.Fix2(significance_category, length(all_motifs))
    colors = cgrad(:sunset, length(sig_categories), categorical=true)
    plt = data(df) * mapping(:patient => nonnumeric, :motif => sorter(all_motifs), color=:significance => sorter(sig_categories) => "significance") * visual(markersize=30)
    return draw(plt; palettes=(color=colors,), draw_kwargs...)
end

function draw_Δμ_plot(df; all_motifs = 1:14 .|> offset_motif_numeral, draw_kwargs...)
    max_Δμ = maximum(abs.(extrema(df.Δμ_z)))
    plt = data(df) * mapping(:patient => nonnumeric, :motif => sorter(all_motifs),  color=:Δμ_z) * visual(markersize=30, colormap=:bluesreds, colorrange=(-max_Δμ,max_Δμ))
    return draw(plt; draw_kwargs...)
end

function draw_Δσ_plot(df; all_motifs = 1:14 .|> offset_motif_numeral, draw_kwargs...)
    max_Δσ = maximum(abs.(extrema(df.Δσ_z)))
    plt = data(df) * mapping(:patient => nonnumeric, :motif => sorter(all_motifs),  color=:Δσ_z) * visual(markersize=30, colormap=:bluesreds, colorrange=(-max_Δσ,max_Δσ))
    return draw(plt; draw_kwargs...)
end