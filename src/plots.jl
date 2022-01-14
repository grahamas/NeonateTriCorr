
using DataFrames, CairoMakie, AlgebraOfGraphics

function plot_eeg_traces(maybe_arr::AbstractArray; labels=nothing, std_max=nothing, sample_rate, downsample_factor=1, layout=:row)
    arr = NamedDimsArray{(:channel, :time)}(maybe_arr)
    arr = if std_max !== nothing
        arr = copy(arr)
        for i_channel ∈ 1:size(arr,1)
            timeseries = arr[i_channel, :]
            std_i = std(timeseries)
            arr[i_channel, abs.(timeseries) .> (std_max * std_i)] .= NaN
        end
        arr
    else
        arr
    end
    time = (0:downsample_factor:(size(arr,2)-1)) ./ sample_rate
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

function draw_eeg_traces(eeg::AbstractProcessedEEG; title=nothing, sample_rate=eeg.sample_rate, kwargs...)
    plt = plot_eeg_traces(eeg; sample_rate=sample_rate, kwargs...)
    axis = (height=50, width=500)
    facet = (; linkyaxes = :none, grid=false, spines=false)
    fg = draw(plt; axis, facet)
    for axis_entry in fg.grid
        ax = axis_entry.axis
        tightlimits!(ax)
        #hidedecorations!(ax)
    end
    if !isempty(eeg.seizure_annotations)
        onsets, offsets = [ann[1] for ann in eeg.seizure_annotations], [ann[2] for ann in eeg.seizure_annotations]
        vspan!.([ae.axis for ae in fg.grid], Ref(onsets), Ref(offsets), color=(:red, 0.2))
    end
    if !isempty(eeg.artifact_annotations)
        onsets, offsets = [ann[1] for ann in eeg.artifact_annotations], [ann[2] for ann in eeg.artifact_annotations]
        vspan!.([ae.axis for ae in fg.grid], Ref(onsets), Ref(offsets), color=(:blue, 0.2))
    end
    if title !== nothing
        Label(fg.figure[0,:], title, tellwidth=false)
    end
    fg
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

function plot_eeg_traces(eeg::AbstractProcessedEEG; kwargs...)
    plot_eeg_traces(eeg.signals; labels=eeg.labels, sample_rate=eeg.sample_rate, kwargs...)
end

function plot_eeg_hists(eeg::AbstractProcessedEEG; kwargs...)
    plot_eeg_hists(eeg.signals; labels=eeg.labels, kwargs...)
end

function plot_contributions(arr::AbstractArray; annotations=nothing, title=nothing, n_motif_classes=14)
    df = DataFrame(arr', :auto)
    @assert ncol(df) == n_motif_classes
    plt = data(df) * mapping(1:n_motif_classes, row=dims(1) => x -> roman_encode(x[1])) * visual(Lines)

    facet = (; linkyaxes = :none, grid=false)
    axis = (height=40, width=500, xlabel="")

    fg = draw(plt; axis, facet)
    axes = [ae.axis for ae in fg.grid]
   # xlabel!.(axes, "")
    tightlimits!.(axes)
    hidespines!.(axes)
    hidedecorations!.(axes, ticklabels=false)
    if annotations !== nothing
        onsets, offsets = calc_seizure_bounds(annotations)
        vspan!.(axes, Ref(onsets), Ref(offsets), color=(:red, 0.2))
    end
    if title !== nothing
        Label(fg.figure[0,:], title, tellwidth=false)
    end
    return fg
end