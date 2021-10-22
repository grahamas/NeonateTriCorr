
using DataFrames, CairoMakie, AlgebraOfGraphics

function plot_eeg_traces(arr::AbstractArray)
    # TODO call NamedDimsArray constructor
    @assert size(arr,2) > size(arr,1) # assume first dim is channel
    tbl = Tables.table(arr', header=["$i" for i ∈ 1:19])
    df = DataFrame(tbl)
    df.time_bin = 1:nrow(df)
    stacked_df = stack(df, 1:19)
    rename!(stacked_df, :value => :signal, :variable => :channel)
    axis = (height=60, width=800)
    aog_data = data(stacked_df)
    plt = aog_data * mapping(:time_bin, :signal, row=:channel) * visual(Lines)
    return (plt, draw(plt; axis))
end
