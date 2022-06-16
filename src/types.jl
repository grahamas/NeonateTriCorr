using TriCorrApplications
using NamedDims

abstract type AbstractProcessedEEG <: AbstractEEG end


struct ProcessedEEGv7{T,SIG<:NamedDimsArray{(:channel,:time),T},ANN_T,ANN_T2,ANN<:Vector{Tuple{ANN_T,ANN_T2}}} <: AbstractProcessedEEG
    signals::SIG
    labels::Vector{String}
    sample_rate::Int
    start::ANN_T #because annotations are whole secs
    duration::Float64
    seizure_annotations::ANN
    artifact_annotations::ANN
    seizure_reviewers_count::Vector{Int}
end

function get_signal_snippet(eeg::AbstractProcessedEEG, snippet_start_sec, snippet_stop_sec)
    i_start = floor(Int, (snippet_start_sec*eeg.sample_rate)+1)
    i_end = floor(Int, (snippet_stop_sec)*eeg.sample_rate)
    get_signal(eeg)[:, i_start:i_end]
end

function TriCorrApplications.get_channel_names(eeg::ProcessedEEGv7)
    eeg.labels
end

function TriCorrApplications.get_times(eeg::ProcessedEEGv7; sample_rate=eeg.sample_rate)
    times = eeg.start:1/sample_rate:(eeg.start+eeg.duration-(1/sample_rate))
    times
end

function set_artifacts_missing(signal::AbstractMatrix, eeg::ProcessedEEGv7{T}; sample_rate=eeg.sample_rate) where T
    times = get_times(eeg, sample_rate=sample_rate)
    artifact_times = in_artifact.(times, Ref(eeg))
    arr = Array{Union{T,Missing}}(signal)
    arr[:, artifact_times] .= missing
    return arr
end
function TriCorrApplications.get_signal(eeg::ProcessedEEGv7)
    sig = set_artifacts_missing(eeg.signals, eeg)
    sig
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

function snip(eeg::ProcessedEEGv7, snip_start_sec::Int, snip_stop_sec)
    @show snip_start_sec snip_stop_sec
    snip_start_idx = round(Int, snip_start_sec * eeg.sample_rate + 1)
    snip_duration = snip_stop_sec - snip_start_sec
    if snip_stop_sec > (eeg.start + eeg.duration)
        @warn "Requested $snip_start_sec thru $snip_stop_sec, greater than duration $(eeg.duration)"
        snip_duration = (eeg.start + eeg.duration) - snip_start_sec
        snip_stop_sec = snip_start_sec + snip_duration
    end
    snip_stop_idx = floor(Int, snip_stop_sec * eeg.sample_rate)
    @show snip_start_idx snip_stop_idx
    ProcessedEEGv7(
        eeg.signals[:,snip_start_idx:snip_stop_idx],
        eeg.labels,
        eeg.sample_rate,
        snip_start_sec,
        snip_duration |> float,
        restrict_ranges(eeg.seizure_annotations, snip_start_sec, snip_stop_sec),
        restrict_ranges(eeg.artifact_annotations, snip_start_sec, snip_stop_sec),
        eeg.seizure_reviewers_count[((snip_start_sec+1):floor(Int, snip_stop_sec))]
    )
end
