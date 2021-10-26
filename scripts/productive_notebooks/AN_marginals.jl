quickactivate("NeonateTriCorr")

include(scriptsdir("include_src.jl"))

file_keys = Dict(
    "pat50_ictal" => ("pat50_1570ICTAL.mat", "Yff", [2, 13, 14]),
    "pat50_interictal" => ("pat50_3100INTERICTAL.mat", "Yff", [2, 13, 14]),
    "pat1_onset" => ("patient1_sz1_10secs.mat", "snippetSave", [18])
)

PAT = "pat50_interictal"

eeg_snippet = load_EEG_snippet(datadir("exp_raw", file_keys[PAT][1]), file_keys[PAT][2])

bad_channels = file_keys[PAT][3]
good_channels = filter(x -> x ∉ bad_channels, axes(eeg_snippet, :channel))
n_time_points = size(eeg_snippet, :time)
n_snippets = 4
λ_max = (11,50)

function normalize_01!(arr)
    arr .-= minimum(arr)
    arr ./= maximum(arr)
    return arr
end
normalize_01!(eeg_snippet)

using TripleCorrelations, Random, Statistics
using ProgressMeter

function mean_nt(nts::AbstractArray{NT}) where {NAMES,NT<:NamedTuple{NAMES}}
    NamedTuple{NAMES}(
        reduce((x,y) -> x .+ y, values.(nts)) ./ length(nts)
    )
end
function divide_nt(nt1::NT, nt2::NT) where {NAMES, NT<:NamedTuple{NAMES}}
    NamedTuple{NAMES}(
        map((x,y) -> x ./ y, values(nt1), values(nt2))
    )
end


marginal_tricorr_an = @showprogress map(1:4) do i_snippet
    i_start = ((i_snippet-1)*n_time_points÷n_snippets)+1
    i_end = i_snippet*n_time_points÷n_snippets
    snippet = eeg_snippet[channel=good_channels, time=i_start:i_end]
    actual = marginal_tricorr_zeropad(snippet, λ_max...)
    noise_controls = [marginal_tricorr_zeropad(shuffle(snippet), λ_max...) for _ ∈ 1:2]
    noise_control = mean_nt(noise_controls)
    divide_nt(actual, noise_control)
end

let time_tricorr_an = [an.time for an in marginal_tricorr_an]
fig = Figure();
layout_coords = (Iterators.product(1:2, 1:2) |> collect)
for i ∈ 1:length(time_tricorr_an)
    ax = Axis(fig[layout_coords[i]...])
    surface!(ax, collect.(axes(time_tricorr_an[i]))..., parent(time_tricorr_an[i]))
    ax.aspect[] = DataAspect()
    ax.xlabel[] = "time lag"
    #ax.zlabel[] = "A/N"
    i_start = ((i-1)*n_time_points÷n_snippets)+1
    i_end = i*n_time_points÷n_snippets
    ax.title[] = "time = $i_start:$i_end"
end
supertitle = Label(fig[0,:], PAT)
display(fig)
save(plotsdir("$(PAT)_time_tricorr_ANshfl_zeropad_$(λ_max)_$(round(Int,time())).png"), fig)
end


# Try plotting all the motif-classes that can be plotted in 2D:
# II, III, IV, V, IV, VII, VIII
# Something could be done with others?

let space_time_tricorr_an = [an.space_time for an in marginal_tricorr_an]
fig = Figure();
layout_coords = (Iterators.product(1:2, 1:2) |> collect)
for i ∈ 1:length(space_time_tricorr_an)
    ax = Axis(fig[layout_coords[i]...])
    surface!(ax, collect.(axes(space_time_tricorr_an[i]))..., parent(space_time_tricorr_an[i]))
    ax.aspect[] = AxisAspect(1)
    ax.xlabel[] = "space lag"
    ax.ylabel[] = "time lag"
    #ax.zlabel[] = "A/N"
    i_start = ((i-1)*n_time_points÷n_snippets)+1
    i_end = i*n_time_points÷n_snippets
    ax.title[] = "time = $i_start:$i_end"
end
supertitle = Label(fig[0,:], PAT)
display(fig)
save(plotsdir("$(PAT)_space_time_tricorr_ANshfl_zeropad_$(λ_max)_$(round(Int,time())).png"), fig)
end


let space_tricorr_an = [an.space for an in marginal_tricorr_an]
fig = Figure();
layout_coords = (Iterators.product(1:2, 1:2) |> collect)
for i ∈ 1:length(space_tricorr_an)
    ax = Axis(fig[layout_coords[i]...])
    surface!(ax, collect.(axes(space_tricorr_an[i]))..., parent(space_tricorr_an[i]))
    ax.aspect[] = DataAspect()
    ax.xlabel[] = "space lag"
    #ax.zlabel[] = "A/N"
    i_start = ((i-1)*n_time_points÷n_snippets)+1
    i_end = i*n_time_points÷n_snippets
    ax.title[] = "time = $i_start:$i_end"
end
supertitle = Label(fig[0,:], PAT)
display(fig)
save(plotsdir("$(PAT)_space_tricorr_ANshfl_zeropad_$(λ_max)_$(round(Int,time())).png"), fig)
end


