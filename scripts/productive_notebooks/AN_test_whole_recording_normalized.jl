quickactivate("NeonateTriCorr")

include(scriptsdir("include_src.jl"))

eeg_snippet = load_EEG_snippet(datadir("exp_raw", "patient1_sz1_10secs.mat"))

bad_channels = [18]
good_channels = filter(x -> x ∉ bad_channels, axes(eeg_snippet, :channel))
n_time_points = size(eeg_snippet, :time)
n_snippets = 4
λ_max = (7,7)

function normalize_01!(arr)
    arr .-= minimum(arr)
    arr ./= maximum(arr)
    return arr
end

using TripleCorrelations
normalize_01!(eeg_snippet)
snippet_tri_corrs_7 = map(1:4) do i_snippet
    i_start = ((i_snippet-1)*n_time_points÷n_snippets)+1
    i_end = i_snippet*n_time_points÷n_snippets
    snippet = eeg_snippet[channel=good_channels, time=i_start:i_end]
    sequence_class_tricorr(snippet, λ_max...)
end

# https://www.rosettacode.org/wiki/Roman_numerals/Encode#Julia
function roman_encode(n::Integer)
    if n < 1 || n > 4999 throw(DomainError(n)) end
 
    DR = [["I", "X", "C", "M"] ["V", "L", "D", "MMM"]]
    rnum = ""
    for (omag, d) in enumerate(digits(n))
        if d == 0
            omr = ""
        elseif d <  4
            omr = DR[omag, 1] ^ d
        elseif d == 4
            omr = DR[omag, 1] * DR[omag, 2]
        elseif d == 5
            omr = DR[omag, 2]
        elseif d <  9
            omr = DR[omag, 2] * DR[omag, 1] ^ (d - 5)
        else
            omr = DR[omag, 1] * DR[omag + 1, 1]
        end
        rnum = omr * rnum
    end
    return rnum
end

log_tricorrs = snippet_tri_corrs_7 .|> x -> sign.(x) .* log10.(abs.(x))

fig = Figure(); fig[1,1]= ax = Axis(fig);
scas = scatter!.(ax, log_tricorrs)
ax.xlabel[] = "motif-class"
ax.ylabel[] = "sign(contribution) * log10(abs(contribution))"
ax.xticks = [1:5:14...]
ax.xtickformat[] = xs -> (offset_motif_numeral ∘ Int).(xs)
Legend(fig[1,2], scas, ["Segment $i" for i in 1:4])
display(fig)
save(plotsdir("log_class_contributions_raw_tricorr_$(round(Int,time())).png"), fig)

noise_snippet = randn(length(good_channels), size(eeg_snippet, :time)÷4)
normalize_01!(noise_snippet)
noise_class_cont = sequence_class_tricorr(noise_snippet, λ_max...)
noise_ratio = snippet_tri_corrs_7 .|> x -> x ./ noise_class_cont

fig = Figure(); fig[1,1]= ax = Axis(fig);
scas = scatter!.(ax, noise_ratio)
ax.xlabel[] = "motif-class"
ax.ylabel[] = "A/N"
ax.xticks = [1:5:14...]
ax.xtickformat[] = xs -> (offset_motif_numeral ∘ Int).(xs)
Legend(fig[1,2], scas, ["Segment $i" for i in 1:4])
display(fig)
save(plotsdir("noise_ratio_raw_tricorr_$(round(Int,time())).png"), fig)