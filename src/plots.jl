
using TriCorrApplications
function TriCorrApplications.plot_reviewer_consensus!(ax, eeg::AbstractProcessedEEG)
    lines!(ax, get_times(eeg, sample_rate=1), eeg.seizure_reviewers_count, color=:grey11)
    tightlimits!(ax); hidespines!(ax)
    ylims!(ax, (0, 3))
    hidedecorations!(ax, ticklabels=false)
end


function plot_estimated_distributions(control_contributions, seizure_contributions; title, all_channel_labels)

    df = DataFrame(mapreduce(vcat, 1:size(control_contributions,1)) do signal_num
        signal_control_contributions = control_contributions[signal_num, :] |> skipmissing |> collect
        signal_seizure_contributions = seizure_contributions[signal_num, :] |> skipmissing |> collect

        signal_contributions = vcat(signal_control_contributions, signal_seizure_contributions)
        signal_kde = kde(signal_contributions)
        common_grid = signal_kde.x
        typeof(signal_kde.x)

        signal_control_kde = kde(signal_control_contributions, common_grid)
        signal_seizure_kde = kde(signal_seizure_contributions, common_grid)

        seizure_nt = [(x=x, density=density, channel=all_channel_labels[signal_num], condition="seizure") for (x, density) in zip(signal_seizure_kde.x, signal_seizure_kde.density)]
        control_nt = [(x=x, density=density, channel=all_channel_labels[signal_num], condition="control") for (x, density) in zip(signal_control_kde.x, signal_control_kde.density)]
        vcat(seizure_nt, control_nt)
    end)

    plt = data(df) * mapping(:x, :density, color=:condition, layout=:channel => sorter(all_channel_labels))
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

function draw_significances_plot!(df; all_channel_labels, axis=(;), draw_kwargs...)
    n_channels = length(all_channel_labels)
    summary = combine(groupby(df, :channel), :Δμ => std, :Δσ => std, :Δμ => mean, :Δσ => mean)
    @show summary
    function normed_Δμ(xs, channels)
        map(zip(xs, channels)) do (x, channel)
            σ = only(filter(:channel => ==(channel), summary).Δμ_std)
            return x / σ
        end
    end
    function normed_Δσ(xs, channels)
        map(zip(xs, channels)) do (x, channel)
            σ = only(filter(:channel => ==(channel), summary).Δσ_std)
            return x / σ
        end
    end
    transform!(df, [:Δμ, :channel] => normed_Δμ => :Δμ_z, [:Δσ, :channel] => normed_Δσ => :Δσ_z)

    reverse_y_scale = y -> (10 ^ (-(10^y)))
    df.significance = significance_category.(df.p, n_channels)
    sig_categories = ([0.1, 0.04, 0.02, 0.002, 0.0002, 0.00002] ./ n_channels) .|> Base.Fix2(significance_category, n_channels)
    colors = cgrad(:sunset, length(sig_categories), categorical=true)
    plt = data(df) * mapping(:channel => sorter(all_channel_labels), :p => (log10 ∘ abs ∘ log10), layout=(:patient => nonnumeric)) * visual(markersize=15)
    return draw(plt; palettes=(color=colors,), draw_kwargs...)
end

# function draw_significances_plot!(df; all_channel_labels, draw_kwargs...)
#     n_channels = length(all_channel_labels)
#     summary = combine(groupby(df, :channel), :Δμ => std, :Δσ => std, :Δμ => mean, :Δσ => mean)
#     @show summary
#     function normed_Δμ(xs, channels)
#         map(zip(xs, channels)) do (x, channel)
#             σ = only(filter(:channel => ==(channel), summary).Δμ_std)
#             return x / σ
#         end
#     end
#     function normed_Δσ(xs, channels)
#         map(zip(xs, channels)) do (x, channel)
#             σ = only(filter(:channel => ==(channel), summary).Δσ_std)
#             return x / σ
#         end
#     end
#     transform!(df, [:Δμ, :channel] => normed_Δμ => :Δμ_z, [:Δσ, :channel] => normed_Δσ => :Δσ_z)

#     df.significance = significance_category.(df.p, n_channels)
#     sig_categories = ([0.1, 0.04, 0.02, 0.002, 0.0002, 0.00002] ./ n_channels) .|> Base.Fix2(significance_category, n_channels)
#     colors = cgrad(:sunset, length(sig_categories), categorical=true)
#     plt = data(df) * mapping(:patient => nonnumeric, :channel => sorter(all_channel_labels), color=:significance => sorter(sig_categories) => "significance") * visual(markersize=30)
#     return draw(plt; palettes=(color=colors,), draw_kwargs...)
# end

# function draw_Δμ_plot(df; all_channel_labels, draw_kwargs...)
#     max_Δμ = maximum(abs.(extrema(df.Δμ_z)))
#     plt = data(df) * mapping(:patient => nonnumeric, :channel => sorter(all_channel_labels),  color=:Δμ_z) * visual(markersize=30, colormap=:bluesreds, colorrange=(-max_Δμ,max_Δμ))
#     return draw(plt; draw_kwargs...)
# end

# function draw_Δσ_plot(df; all_channel_labels, draw_kwargs...)
#     max_Δσ = maximum(abs.(extrema(df.Δσ_z)))
#     plt = data(df) * mapping(:patient => nonnumeric, :channel => sorter(all_channel_labels),  color=:Δσ_z) * visual(markersize=30, colormap=:bluesreds, colorrange=(-max_Δσ,max_Δσ))
#     return draw(plt; draw_kwargs...)
# end

function draw_Δμ_plot(df; all_channel_labels, draw_kwargs...)
    n_channels = length(all_channel_labels)
    sig_categories = ([0.1, 0.04, 0.02, 0.002, 0.0002, 0.00002] ./ n_channels) .|> Base.Fix2(significance_category, n_channels)
    colors = cgrad(:sunset, length(sig_categories), categorical=true)
    max_Δμ = maximum(abs.(extrema(df.Δμ_z)))
    plt = data(df) * mapping(:channel => sorter(all_channel_labels), :Δμ_z 
=> "Δμ z-score", layout= :patient => nonnumeric, color=:significance => sorter(sig_categories) => "significance") * visual(Stem, markersize=15)
    return draw(plt; palettes=(color=colors,), legend=(markersize=15,), draw_kwargs...)
end

function draw_Δσ_plot(df; all_channel_labels, draw_kwargs...)
    n_channels = length(all_channel_labels)
    sig_categories = ([0.1, 0.04, 0.02, 0.002, 0.0002, 0.00002] ./ n_channels) .|> Base.Fix2(significance_category, n_channels)
    colors = cgrad(:sunset, length(sig_categories), categorical=true)
    max_Δσ = maximum(abs.(extrema(df.Δσ_z)))
    plt = data(df) * mapping(:channel => sorter(all_channel_labels), :Δσ_z => "Δσ z-score", layout= :patient => nonnumeric, color=:significance => sorter(sig_categories) => "significance") * visual(Stem, markersize=15)
    return draw(plt; palettes=(color=colors,), legend=(markersize=15,), draw_kwargs...)
end

using Makie: convert_arguments
Makie.convert_arguments(P::PointBased, r::Roc) = convert_arguments(P, r.pfa, 1 .- r.pmiss)

function plot_roc(r::Roc; xlabel="FPR", ylabel="TPR", title="AUC: $(1 - auc(r))")
    (fig, ax, plt) = lines(r)
    ax.title[] = title
    ax.xlabel[] = xlabel
    ax.ylabel[] = ylabel
    tightlimits!(ax)
    (fig, ax, plt)
end

function plot_roc_fphr(r::Roc; xlabel="FP/hr", ylabel="TPR", title="AUC: $(round(1 - auc(r), sigdigits=3))", resolution=(800,600))
    fig = Figure(resolution=resolution)
    fig[1,1] = ax = Axis(fig; xlabel=xlabel, ylabel=ylabel, title=title)
    plt = plot_roc_fphr!(ax, r)
    (fig, ax, plt)
end

function plot_roc_fphr!(ax::Axis, r::Roc; )
    plt = lines!(ax, r.pfa .* 60, 1 .- r.pmiss)
    tightlimits!(ax)
    return plt
end

function tpr_at_fphr(r::Roc, fphr)
    idx = findfirst(r.pfa .* 60. .<= fphr)
    return 1 - r.pmiss[idx]
end

function fphr_at_tpr(r::Roc, tpr)
    pmiss = 1 - tpr
    idx = findfirst(r.pmiss .>= pmiss)
    return r.pfa[idx] .* 60.
end