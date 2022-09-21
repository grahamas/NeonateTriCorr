function calc_seizure_snippet_starts(eeg::AbstractProcessedEEG, snippets_duration_s)
    calc_seizure_snippet_starts(eeg.seizure_annotations, snippets_duration_s)
end
function calc_seizure_snippet_starts(seizure_annotations, snippets_duration_s)
    return vcat([collect(on:snippets_duration_s:off) for (on, off) ∈ seizure_annotations]...)
end

function calc_control_snippet_starts(eeg::AbstractProcessedEEG, snippets_duration_s, min_dist_to_seizure)
    # should assert that seizure and artifacts do not overlap
    control_bounds = calc_control_bounds(eeg, snippets_duration_s, min_dist_to_seizure)
    return vcat([on:snippets_duration_s:off for (on, off) ∈ control_bounds]...)
end


function calc_control_bounds(eeg::AbstractProcessedEEG, snippets_duration_s, min_dist_to_seizure)
    non_control_bounds = merge_bounds(eeg.seizure_annotations, eeg.artifact_annotations)
    control_onsets = [0, first.(non_control_bounds)...]
    control_offsets = [last.(non_control_bounds)..., eeg.duration]
    control_bounds = zip(control_onsets, control_offsets)
    @show control_bounds
    return [(on+min_dist_to_seizure, off-min_dist_to_seizure) for (on, off) ∈ control_bounds if (off-min_dist_to_seizure)-(on+min_dist_to_seizure) > snippets_duration_s]
end

function calc_control_snippets_incl_artifacts(seizure_annotations, recording_duration, signal_times, min_dist_to_seizure)
    # should assert that seizure and artifacts do not overlap
    seizure_onsets = first.(seizure_annotations)
    seizure_offsets = last.(seizure_annotations)
    control_onsets = [0, seizure_offsets...]
    sort!(control_onsets)
    control_offsets = [seizure_onsets..., recording_duration]
    sort!(control_offsets)
    control_bounds = zip(control_onsets, control_offsets)
    return reduce((x,y) -> x .|| y, [(on+min_dist_to_seizure) .<= (signal_times .+ 1) .< (off-min_dist_to_seizure) for (on, off) ∈ control_bounds])
end

function control_vs_seizure_class_contributions(eeg::AbstractProcessedEEG; 
    boundary,
    contributions_desc,
    n_snippets, snippets_duration_s=1, min_dist_to_seizure=60, kwargs...)
    seizure_snippet_starts = calc_seizure_snippet_starts(eeg, snippets_duration_s)
    control_snippet_starts = calc_control_snippet_starts(eeg, snippets_duration_s, min_dist_to_seizure)

    if (length(control_snippet_starts) <= n_snippets)
        @warn "Patient lacks sufficiently many separated snippets."
        return missing
    end

    seizure_snippet_starts = sort(sample(seizure_snippet_starts, min(n_snippets, length(seizure_snippet_starts)), replace=false))
    control_snippet_starts = sort(sample(control_snippet_starts, n_snippets, replace=false))

    seizure_snippets_class_contribution = calc_class_contributions(eeg, boundary, contributions_desc; snippets_start_sec = seizure_snippet_starts, snippets_duration_s = snippets_duration_s, kwargs...)
    control_snippets_class_contribution = calc_class_contributions(eeg, boundary, contributions_desc; snippets_start_sec = control_snippet_starts, snippets_duration_s = snippets_duration_s, kwargs...)

    return (seizure=seizure_snippets_class_contribution, control=control_snippets_class_contribution)
end

function putative_signal_rms(arr; putative_signal_classes=2:14)
    arr = NamedDimsArray{(:motif_class, :time)}(arr)
    signal_arr = arr[putative_signal_classes,:]
    dropdims(sqrt.(mean(signal_arr .^ 2, dims=1)), dims=:motif_class)
end

function all_class_boxplots_control_vs_seizure(dct; show_outliers=true)
    control = dct[:control]; seizure = dct[:seizure]
    header = [offset_motif_numeral.(1:14)..., "case"]
    values = [control'; seizure']
    labels = [[1 for _ ∈ 1:size(control,:time)]...; [2 for _ ∈ 1:size(seizure,:time)]]
    labeled_values = hcat(values, labels)
    df = DataFrame(labeled_values, header)
    stacked_df = stack(df, 1:14)
    plt = data(stacked_df) * mapping(:case, :value => "A/N - 1", layout=:variable) * visual(BoxPlot; show_outliers=show_outliers)
    axis = (xticks = (1:2, ["control", "seizure"]), height=120, width=120)
    facet = (; linkyaxes = :none)
    draw(plt; axis, facet)
end

function boxplot_control_vs_seizure(dct)
    control = NamedDimsArray{(:time,)}(dct[:control]); 
    seizure = NamedDimsArray{(:time,)}(dct[:seizure])
    values = [control...; seizure...]
    labels = [[1 for _ ∈ control]...; [2 for _ ∈ seizure]]
    fig = Figure(); ax = Axis(fig[1,1])
    boxplot!(ax, labels, values; show_outliers=false)
    ax.xticks = (1:2, ["control", "seizure"])
    fig
end

function control_vs_seizure_all_class_statistics(seizure::NamedDimsArray{(:motif_class,:time)}, control::NamedDimsArray{(:motif_class,:time)})
    DataFrame(map(axes(seizure, :motif_class)) do class_i
        seizure_i = seizure[class_i,:]
        control_i = control[class_i,:]
        return (
            class = offset_motif_numeral(class_i),
            pvalue = pvalue(MannWhitneyUTest(seizure_i, control_i)),
            mean_effect = mean(seizure_i) - mean(control_i)
        )
    end)
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