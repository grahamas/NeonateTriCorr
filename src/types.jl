using NamedDims, Dates

abstract type AbstractProcessedEEG end

struct ProcessedEEGv7{T,SIG<:NamedDimsArray{(:channel,:time),T},ANN_T,ANN<:Vector{NTuple{2,ANN_T}}} <: AbstractProcessedEEG
    signals::SIG
    labels::Vector{String}
    sample_rate::Int
    start::Int
    duration::Float64
    seizure_annotations::ANN
    artifact_annotations::ANN
    seizure_reviewers_count::Vector{Int}
end

function get_times(eeg::ProcessedEEGv7; sample_rate=eeg.sample_rate)
    eeg.start:1/sample_rate:(eeg.start+eeg.duration-(1/sample_rate))
end

function set_artifacts_missing(signal::AbstractMatrix, eeg::ProcessedEEGv7{T}; sample_rate=eeg.sample_rate) where T
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

function in_artifact(sec, eeg::ProcessedEEGv7)
    artifacts = eeg.artifact_annotations
    return any(in_tuple(sec, artifact) for artifact in artifacts)
end

function restrict_ranges(tups, start, stop)
    [(max(on,start),min(off,stop)) for (on,off) ∈ tups if on ∈ start..stop || off ∈ start..stop]
end

function snip_start(eeg::ProcessedEEGv7, snip_start_sec::Int)
    snip(eeg, snip_start_sec, eeg.duration)
end

function snip(eeg::ProcessedEEGv7, duration)
    snip(eeg, eeg.start, duration)
end

function snip(eeg::ProcessedEEGv7, snip_start_sec::Int, snip_stop_sec::Int)
    snip_start_idx = snip_start_sec * eeg.sample_rate + 1
    snip_duration = snip_stop_sec - snip_start_sec
    if snip_stop_sec > (eeg.start + eeg.duration)
        @warn "Requested $snip_start_sec thru $snip_stop_sec, greater than duration $(eeg.duration)"
        snip_duration = (eeg.start + eeg.duration) - snip_start_sec
        snip_stop_sec = snip_start_sec + snip_duration
    end
    snip_stop_idx = floor(Int, snip_stop_sec * eeg.sample_rate)
    ProcessedEEGv7(
        eeg.signals[:,snip_start_idx:snip_stop_idx],
        eeg.labels,
        eeg.sample_rate,
        snip_start_sec,
        snip_duration |> float,
        restrict_ranges(eeg.seizure_annotations, snip_start_sec, snip_stop_sec),
        restrict_ranges(eeg.artifact_annotations, snip_start_sec, snip_stop_sec),
        eeg.seizure_reviewers_count[(snip_start_sec:snip_stop_sec) .+ 1]
    )
end