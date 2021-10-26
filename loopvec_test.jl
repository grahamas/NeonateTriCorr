using LoopVectorization

src = ones(19, 101)
axis_a, axis_b = axes(src)
max_a = 7; max_b = 9
range_a = -max_a:max_a
range_b = -max_b:max_b

target_arr = zeros(Float64, length(range_b), length(range_b))

padded_axis_a = (first(axis_a) .+ max_a):(last(axis_a) - max_a)
padded_axis_b = (first(axis_b) .+ max_b):(last(axis_b) - max_b)

@turbo for a1i ∈ eachindex(range_a), 
            a2i ∈ eachindex(range_a), 
            b1i ∈ eachindex(range_b), 
            b2i ∈ eachindex(range_b)
    a1 = range_a[a1i]; 
    a2 = range_a[a2i] 
    b1 = range_b[b1i]; b2 = range_b[b2i]
    contribution = target_arr[b1i, b2i]
    for i_a ∈ padded_axis_a, i_b ∈ padded_axis_b
        contribution += src[i_a, i_b] * src[i_a,i_b+b1] * src[i_a,i_b+b2]
    end
    target_arr[b1i, b2i] = contribution
end