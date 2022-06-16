
abstract type AbstractRollingCalculation end

function roll_through(t::Type{<:AbstractRollingCalculation}, vec::AbstractVector{T}, N::Int) where T
    rolling = initialize(t, vec, N)
    values = zeros(T, length(vec) - N + 1)
    values[1] = rolling.value
    for i ∈ (N+1):length(vec)
        values[i] = update!(rolling, vec[i])
    end
    return values
end


mutable struct RollingMean{T, V} <: AbstractRollingCalculation
    value::T
    N::Int
    xs::V
    i::Int
end
function update!(rmean::RollingMean{T}, x::T) where T
    rmean.value += (x - rmean.xs[rmean.i]) / rmean.N
    rmean.xs[rmean.i] = x
    rmean.i = (rmean.i % rmean.N) + 1
    return rmean.value
end
function initialize(::Type{<:RollingMean}, vec::V, N::Int) where {T, V<:AbstractVector{T}}
    RollingMean{T,V}(
        mean(vec[1:N]),
        N,
        copy(vec[1:N]),
        1
    )
end

mutable struct RollingStandardDeviation{T,V} <: AbstractRollingCalculation
    mean::RollingMean{T}
    value::T
    N::Int
    xs::V
    i::Int
end
function update!(rstd::RollingStandardDeviation{T}, x::T) where T
    old_mean = rstd.mean.value
    update!(rstd.mean, x)
    rstd.value += (x - rstd.xs[rstd.i]) * (x - rstd.mean.value + rstd.xs[rstd.i] - old_mean)
    rstd.xs[rstd.i] = x
    rstd.i = (rstd.i % rstd.N) + 1
    return sqrt(rstd.value / (rstd.N-1))
end
function initialize(::Type{<:RollingStandardDeviation}, vec::V, N::Int) where {T, V<:AbstractVector{T}}
    RollingStandardDeviation{T,V}(
        initialize(RollingMean, vec, N),
        std(vec[1:N]),
        N,
        copy(vec[1:N]),
        1
    )
end
