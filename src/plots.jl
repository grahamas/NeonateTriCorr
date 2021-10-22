
using DataFrames, CairoMakie, AlgebraOfGraphics

function plot_eeg_traces(maybe_arr::AbstractArray)
    arr = NamedDimsArray{(:channel, :time)}(maybe_arr)
    tbl = Tables.table(arr', header=["$i" for i ∈ axes(arr, :channel)])
    df = DataFrame(tbl)
    df.time_bin = 1:nrow(df)
    stacked_df = stack(df, axes(arr, :channel))
    rename!(stacked_df, :value => :signal, :variable => :channel)
    axis = (height=60, width=800)
    aog_data = data(stacked_df)
    plt = aog_data * mapping(:time_bin, :signal, row=:channel) * visual(Lines)
    return (plt, draw(plt; axis))
end
