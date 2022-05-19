using Makie

# ANY SCRIPT USING THIS MUST IMPORT A MAKIE
using DataFrames, AlgebraOfGraphics

function plot_eeg_traces(eeg::AbstractProcessedEEG; labels=get_channel_names(eeg), std_max=nothing, sample_rate=eeg.sample_rate, downsample_factor=1, layout=:row)
    arr = NamedDimsArray{(:channel, :time)}(get_signal(eeg))
    arr = if std_max !== nothing
        arr = copy(arr)
        for i_channel ∈ 1:size(arr,1)
            timeseries = arr[i_channel, :]
            std_i = std(dropmissing(timeseries))
            arr[i_channel, abs.(timeseries) .> (std_max * std_i)] .= NaN
        end
        arr
    else
        arr
    end
    time = get_times(eeg)[begin:downsample_factor:end]
    arr = arr[:,begin:downsample_factor:end]
    tbl = if labels === nothing
        Tables.table(arr', header=["$i" for i ∈ axes(arr, :channel)])
    else
        Tables.table(arr', header=labels)
    end
    df = DataFrame(tbl)
    df.time = time
    stacked_df = stack(df, axes(arr, :channel))
    rename!(stacked_df, :value => :signal, :variable => :channel)
    aog_data = data(stacked_df)
    plt = if layout == :row
        aog_data * mapping(:time, :signal, row=:channel) * visual(Lines)
    elseif layout == :layout
        aog_data * mapping(:time, :signal, layout=:channel) * visual(Lines)
    else
        error("bad layout kwarg")
    end
    return plt
end

function draw_eeg_traces(eeg::AbstractProcessedEEG; title=nothing, resolution, kwargs...)
    fig = Figure(resolution=resolution)
    times = get_times(eeg)
    for i_sig ∈ 1:size(eeg.signals, 1)
        ax = Axis(fig[i_sig,1])
        plot_contribution!(ax, times, eeg.signals[i_sig,:]; eeg=eeg)
        Label(fig[i_sig,2], eeg.labels[i_sig], tellwidth=true, tellheight=false, rotation=-pi/2)
    end
    if title !== nothing
        Label(fig[0,:], title, tellwidth=false)
    end
    fig
end

function plot_eeg_hists(maybe_arr::AbstractArray; labels=nothing)
    arr = NamedDimsArray{(:channel, :time)}(maybe_arr)
    tbl = if labels === nothing
        Tables.table(arr', header=["$i" for i ∈ axes(arr, :channel)])
    else
        Tables.table(arr', header=labels)
    end
    df = DataFrame(tbl)
    df.time = 1:nrow(df)
    stacked_df = stack(df, axes(arr, :channel))
    rename!(stacked_df, :value => :signal, :variable => :channel)
    aog_data = data(stacked_df)
    plt = aog_data * mapping(:signal, row=:channel) * AlgebraOfGraphics.histogram(bins=200)
    return plt
end

function draw_eeg_hists(arr::AbstractProcessedEEG; title=nothing, kwargs...)
    plt = plot_eeg_hists(arr; kwargs...)
    axis = (height=40, width=400)
    facet = (; linkyaxes = :none, linkxaxes=:none, grid=false, spines=false)
    fg = draw(plt; axis, facet)
    for axis_entry in fg.grid
        ax = axis_entry.axis
        tightlimits!(ax)
        #hidedecorations!(ax)
    end
    if title !== nothing
        Label(fg.figure[0,:], title, tellwidth=false)
    end
    fg
end

function plot_eeg_hists(eeg::AbstractProcessedEEG; kwargs...)
    plot_eeg_hists(eeg.signals; labels=eeg.labels, kwargs...)
end

plot_contributions(vec::AbstractVector; kwargs...) = plot_contributions(DataFrame([vec], :auto); kwargs...)
plot_contributions(arr::AbstractMatrix; kwargs...) = plot_contributions(DataFrame(arr', :auto); kwargs...)

noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")

function plot_contributions(times::AbstractVector, data::AbstractArray, eeg::AbstractProcessedEEG; title=nothing, resolution)
    n_rows = size(data, 1)
    fig = Figure(resolution = resolution, font = noto_sans)
    ax_plt_pairs = map(1:n_rows) do i_row
        ax = Axis(fig[i_row, 1], xlabel="", xgridvisible=false, ygridvisible=false)
        plt = plot_contribution!(ax, times, data[i_row,:]; 
            eeg=eeg
        )
        motif=offset_motif_numeral(i_row)
        Label(fig[i_row, 2], motif, tellheight=false, tellwidth=true, rotation=-pi/2)
        (ax, plt)
    end
    if title !== nothing
        Label(fig[0,:], title, tellwidth=false)
    end
    return fig
end

function plot_contribution(args...; title=nothing, resolution=(800,600),eeg=nothing, kwargs...)
    fig = Figure(resolution=resolution)
    ax = Axis(fig[1,1])
    cons_ax = Axis(fig[2,1])
    l = plot_contribution!(ax, args...; eeg=eeg, kwargs...)
    if title !== nothing
        Label(fig[0,:], title, tellwidth=false)
    end
    plot_reviewer_consensus!(cons_ax, eeg)
    linkxaxes!(ax, cons_ax)
    rowsize!(fig.layout, 3, Auto(0.3))
    return (fig, ax, l)
end

function plot_contribution!(ax, times::AbstractVector, data::AbstractVector; eeg=nothing)
    l = lines!(ax, times, data)
    tightlimits!(ax); hidespines!(ax)
    hidedecorations!(ax, ticklabels=false)
    if eeg !== nothing
        seizure_onsets = [on for (on,off) in eeg.seizure_annotations if on < off]
        seizure_offsets = [off for (on,off) in eeg.seizure_annotations if on < off]
        artifact_onsets = [on for (on,off) in eeg.artifact_annotations if on < off]
        artifact_offsets = [off for (on,off) in eeg.artifact_annotations if on < off]
        if !isempty(artifact_onsets)
            vspan!(ax, artifact_onsets, artifact_offsets, color=(:red, 0.2))
        end
        if !isempty(seizure_onsets)
            @warn "No seizures."
            vspan!(ax, seizure_onsets, seizure_offsets, color=(:blue, 0.2))
        end
    end
    l
end

function plot_reviewer_consensus!(ax, eeg::AbstractProcessedEEG)
    lines!(ax, get_times(eeg, sample_rate=1), eeg.seizure_reviewers_count)
    tightlimits!(ax); hidespines!(ax)
    hidedecorations!(ax, ticklabels=false)
end