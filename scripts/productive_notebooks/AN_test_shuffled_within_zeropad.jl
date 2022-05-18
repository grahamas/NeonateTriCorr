quickactivate("NeonateTriCorr")

include(scriptsdir("include_src.jl"))

eeg_snippet = load_EEG_snippet(datadir("exp_raw", "patient1_sz1_10secs.mat"))

bad_channels = [18]
good_channels = filter(x -> x ∉ bad_channels, axes(eeg_snippet, :channel))
n_time_points = size(eeg_snippet, :time)
n_snippets = 4
λ_max = (9,25)

function normalize_01!(arr)
    arr .-= minimum(arr)
    arr ./= maximum(arr)
    return arr
end

using TripleCorrelations, Random, Statistics
normalize_01!(eeg_snippet)
snippet_ans = map(1:4) do i_snippet
    i_start = ((i_snippet-1)*n_time_points÷n_snippets)+1
    i_end = i_snippet*n_time_points÷n_snippets
    snippet = eeg_snippet[channel=good_channels, time=i_start:i_end]
    actual = sequence_class_tricorr_zeropad(snippet, λ_max...)
    noise_control = mean(sequence_class_tricorr_zeropad(shuffle(snippet), λ_max...) for _ ∈ 1:2)
    actual ./ noise_control
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



fig = Figure(); fig[1,1]= ax = Axis(fig);
scas = scatter!.(ax, snippet_ans)
ax.xlabel[] = "motif-class"
ax.ylabel[] = "A/N"
ax.xticks = [1:5:14...]
ax.xtickformat[] = xs -> (offset_motif_numeral ∘ Int).(xs)
Legend(fig[1,2], scas, ["Segment $i" for i in 1:4])
display(fig)
save(plotsdir("an_shuffled_raw_tricorr_zeropad_$(λ_max)_$(round(Int,time())).png"), fig)