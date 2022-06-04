
function prep(func::Function, args...)
    func
end

function thrEpredStdNorm!(output, input, σ)
    # theoretically derived E predicated on standard normal inputs
    output .= input ./ σ
end
function prep(::typeof(thrEpredStdNorm!), raster_size, boundary, lag_extents)
    σ = estimate_std_of_standard_normals(raster_size, boundary, lag_extents)
    (o, i) -> thrEpredStdNorm!(o, i, σ)
end
