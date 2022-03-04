using NamedDims, Dates

abstract type AbstractProcessedEEG end

struct ProcessedEEGv6{T,SIG<:NamedDimsArray{(:channel,:time),T},ANN_T,ANN<:Vector{NTuple{2,ANN_T}}} <: AbstractProcessedEEG
    signals::SIG
    labels::Vector{String}
    sample_rate::Int
    start_time::DateTime
    start_sec::Int
    duration::Float64
    seizure_annotations::ANN
    artifact_annotations::ANN
end

function get_times(eeg::ProcessedEEGv6; sample_rate=eeg.sample_rate)
    eeg.start_sec:1/sample_rate:(eeg.start_sec+eeg.duration-(1/sample_rate))
end

function set_artifacts_missing(signal::AbstractMatrix, eeg::ProcessedEEGv6{T}; sample_rate=eeg.sample_rate) where T
    times = get_times(eeg, sample_rate=sample_rate)
    artifact_times = in_artifact.(times, Ref(eeg))
    arr = Array{Union{T,Missing}}(signal)
    arr[:, artifact_times] .= missing
    return arr
end
function get_signal(eeg)
    set_artifacts_missing(eeg.signals, eeg)
end



function in_tuple(sec, (start,stop)::NTuple{2})
    start <= sec < stop
end

function in_artifact(sec, eeg::ProcessedEEGv6)
    artifacts = eeg.artifact_annotations
    return any(in_tuple(sec, artifact) for artifact in artifacts)
end

function restrict_ranges(tups, start, stop)
    [(max(on,start),min(off,stop)) for (on,off) ∈ tups if on ∈ start..stop || off ∈ start..stop]
end

function snip_start(eeg::ProcessedEEGv6, snip_start_sec::Int)
    snip(eeg, snip_start_sec, ceil(Int, eeg.duration))
end

function snip(eeg::ProcessedEEGv6, duration)
    snip(eeg, eeg.start_sec, duration)
end

function snip(eeg::ProcessedEEGv6, snip_start_sec::Int, snip_stop_sec::Int)
    @show snip_start_sec snip_stop_sec eeg.duration
    @show count(ismissing.(eeg.signals))
    snip_start_idx = snip_start_sec * eeg.sample_rate + 1
    snip_duration = snip_stop_sec - snip_start_sec
    if snip_duration > eeg.duration
        @warn "Requested $snip_start_sec thru $snip_stop_sec, greater than duration $(eeg.duration)"
    end
    snip_stop_idx = snip_duration * eeg.sample_rate
    ProcessedEEGv6(
        eeg.signals[:,snip_start_idx:snip_stop_idx],
        eeg.labels,
        eeg.sample_rate,
        eeg.start_time + Second(snip_start_sec),
        snip_start_sec,
        snip_duration |> float,
        restrict_ranges(eeg.seizure_annotations, snip_start_sec, snip_stop_sec),
        restrict_ranges(eeg.artifact_annotations, snip_start_sec, snip_stop_sec)
    )
end