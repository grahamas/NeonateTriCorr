
using DataFrames, CairoMakie, AlgebraOfGraphics

function plot_eeg_traces(maybe_arr::AbstractArray; labels=nothing)
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
    plt = aog_data * mapping(:time, :signal, row=:channel) * visual(Lines)
    return plt
end

function draw_eeg_traces(arr::Union{AbstractArray,ProcessedEEGv1}; kwargs...)
    plt = plot_eeg_traces(arr)
    axis = (height=60, width=800)
    facet = (; linkyaxes = :none, grid=false, spines=false)
    fg = draw(plt; axis, facet, kwargs...)
    for axis_entry in fg.grid
        ax = axis_entry.axis
        tightlimits!(ax)
        hidedecorations!(ax)
    end
    fg
end

function plot_eeg_traces(eeg::ProcessedEEGv1)
    plot_eeg_traces(eeg.signals; labels=eeg.labels)
end
