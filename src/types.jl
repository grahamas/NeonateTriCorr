using NamedDims, Dates

abstract type AbstractProcessedEEG end

struct ProcessedEEGv5{T,SIG<:NamedDimsArray{(:channel,:time),T},ANN_T,ANN<:Vector{NTuple{2,ANN_T}}} <: AbstractProcessedEEG
    signals::SIG
    labels::Vector{String}
    sample_rate::Int
    start::DateTime
    duration::Float64
    seizure_annotations::ANN
    artifact_annotations::ANN
end

function get_times(eeg::ProcessedEEGv5{T}) where T
    0.:1/eeg.sample_rate:(eeg.duration-1/eeg.sample_rate)
end
function get_signal(eeg::ProcessedEEGv5{T}) where T
    times = get_times(eeg)
    @show length(times)
    artifact_times = in_artifact.(times, Ref(eeg))
    arr = Array{Union{T,Missing}}(eeg.signals)
    arr[:, artifact_times] .= missing
    return arr
end

function in_tuple(sec, (start,stop)::NTuple{2})
    start <= sec < stop
end

function in_artifact(sec, eeg::ProcessedEEGv5)
    artifacts = eeg.artifact_annotations
    return any(in_tuple(sec, artifact) for artifact in artifacts)
end

function restrict_ranges(tups, start, stop)
    [(max(on,start),min(off,stop)) for (on,off) ∈ tups if on ∈ start..stop || off ∈ start..stop]
end

function snip(eeg::ProcessedEEGv5, duration)
    snip(eeg, 0, duration)
end

function snip(eeg::ProcessedEEGv5, snip_start_time::DateTime, snip_stop_time::DateTime)
    snip(eeg, snip_start_time - eeg.start, snip_stop_time - eeg.start)
end

function snip(eeg::ProcessedEEGv5, snip_start_sec::Int, snip_stop_sec::Int)
    snip_start_idx = snip_start_sec * eeg.sample_rate + 1
    snip_duration = snip_stop_sec - snip_start_sec
    if snip_duration > eeg.duration
        @warn "Requested $snip_start_sec thru $snip_stop_sec, greater than duration $(eeg.duration)"
    end
    snip_stop_idx = snip_duration * eeg.sample_rate
    ProcessedEEGv5(
        eeg.signals[:,snip_start_idx:snip_stop_idx],
        eeg.labels,
        eeg.sample_rate,
        eeg.start + Second(snip_start_sec),
        snip_duration |> float,
        restrict_ranges(eeg.seizure_annotations, snip_start_sec, snip_stop_sec),
        restrict_ranges(eeg.artifact_annotations, snip_start_sec, snip_stop_sec)
    )
end