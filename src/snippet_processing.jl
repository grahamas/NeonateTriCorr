
function extrema_01!(output, input)
    mn, mx = extrema(input)
    @. output = (input - mx) / mn
end

function zscore!(output, input, μ=mean(input), σ=std(input))
    @. output = (input - μ) / σ
end