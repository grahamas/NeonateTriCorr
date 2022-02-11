using NamedDims, Dates

abstract type AbstractProcessedEEG end

struct ProcessedEEGv3{T,SIG<:NamedDimsArray{(:channel,:time),T}} <: AbstractProcessedEEG
    signals::SIG
    labels::Vector{String}
    sample_rate::Float64
    start::DateTime
    duration::Float64
    seizure_annotations::Vector{Tuple{Float64,Float64}}
    artifact_annotations::Vector{Tuple{Float64,Float64}}
end