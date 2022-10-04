using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, xticklabelsize=20, yticklabelsize=20, linecolor=:black, linewidth=4, labelsize=36)
set_theme!(font_theme)


using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity
using ROCAnalysis

include(scriptsdir("include_src.jl"))

let signals_reduction_name = "maxany",
    patients_considered = [9,12,47,62],
    resolution=(2000,1600);

signal_types = ["tricorr", "aEEG"]

sr_dict = Dict(
    "aEEG" => signals_reduction_name,
    "tricorr" => "$(signals_reduction_name)abs"
)
colors = Dict(
    "tricorr" => :blue,
    "aEEG" => :red
)

save_dir = plotsdir("fig4_signals_$(Dates.now())")
mkpath(save_dir)

epoch_s = common_params[:epoch_s]

signals = Dict(
    pat_num => Dict(
        sig => load_signal(; common_params..., analysis_particular_params[sig]..., patient_num=pat_num) for sig in signal_types
    ) for pat_num in patients_considered
)

eegs = Dict(
    pat_num => load_helsinki_eeg(pat_num; common_params..., excluded_artifact_grades=Int[1]) for pat_num in patients_considered
    ) 

fig = Figure(resolution=resolution, figure_padding=(10,30,10,30))

axes = Array{Any}(undef, 2, 3, 2)

map(eachindex(axes)) do idx
    axes[idx] = Axis(fig)
end

eeg_nums = [
    9 62 12;
    9 47 12
]

function snip_interval((on, off), eeg, time, signal, min_len=1000)
    dur = off - on
    wide_on = max(0, on - (min_len - dur)/2)
    wide_off = min(off + (min_len - (off - wide_on)), time[end])
    time_idxs = wide_on .<= time .<= wide_off
    snipped_time = time[time_idxs]
    return (snip(eeg, snipped_time[begin], snipped_time[end]), time[time_idxs], signal[time_idxs])
end
function snip_seizure(eeg, time, signal, sz_num, min_len=1000)
    interval = eeg.seizure_annotations[sz_num]
    snip_interval(interval, eeg, time, signal, min_len)
end
function snip_artifact(eeg, time, signal, art_num, min_len=1000)
    interval = eeg.artifact_annotations[art_num]
    snip_interval(interval, eeg, time, signal, min_len)
end

patients_snipping = [
    ((a...) -> snip_seizure(a..., 1)) ((a...) -> snip_artifact(a..., 1)) ((a...) -> snip_interval((3885, 3938), a...));
    ((a...) -> snip_seizure(a..., 2)) ((a...) -> snip_artifact(a..., 2)) ((a...) -> snip_interval((2000, 2053), a...))
]



for idx ∈ CartesianIndices(eeg_nums)
    eeg_num = eeg_nums[idx]
    confusion_entry = GridLayout()
    snip_fn = patients_snipping[idx]
    eeg = eegs[eeg_num]
for (i_st, st) ∈ enumerate(signal_types)
    raw_signal = signals[eeg_num][st]
    rolled_signal = roll_signals(raw_signal; common_params..., analysis_particular_params[st]...)
    not_missings = .!ismissing.(rolled_signal[1,:])
    rolled_signal .-= mean(rolled_signal[:, not_missings], dims=:time)
    rolled_signal ./= std(rolled_signal[:, not_missings], dims=:time)
    reduce_signals_fn = get_reduce_signals_fn(sr_dict[st])
    reduced_signal = reduce_signals_fn(rolled_signal)

    dt = analysis_particular_params[st][:snippets_duration_s]
    time = 0.:dt:((length(reduced_signal)-1)*dt)

    snipped_eeg, snipped_time, snipped_signal = snip_fn(eeg, time, reduced_signal)

    TriCorrApplications.plot_contribution!(axes[idx, i_st], snipped_eeg, snipped_time, snipped_signal; color=colors[st])
    confusion_entry[i_st, 1] = axes[idx, i_st]
end
hidexdecorations!(axes[idx, 1])

fig[(Tuple(idx) .+ (2,2))...] = confusion_entry
end

tightlimits!.(axes, Ref(Left()), Ref(Right()))
ylims!.(axes[:,:,1], Ref((0, 7)))
ylims!.(axes[:,:,2], Ref((-1.,5.5)))

hidedecorations!.(axes; ticklabels=false, label=false)
hideydecorations!.(axes[:, 2:3, :])

axes[1,1,1].ylabel[] = "tricorr"
axes[2,1,1].ylabel[] = "tricorr"

axes[1,1,2].ylabel[] = "aEEG"
axes[2,1,2].ylabel[] = "aEEG"

axes[2,1,2].xlabel[] = "time"

fig[1,3:5] = Label(fig, "clinician determination", tellwidth=false)
fig[2,3:5] = [Label(fig, "seizure", tellwidth=false)
    Label(fig, "artifact", tellwidth=false)
    Label(fig, "nothing", tellwidth=false)]

fig[3:4,1] = Label(fig, "tricorr detector determination", tellheight=false, rotation=pi/2)
fig[3,2] = Label(fig, "detection", tellheight=false, rotation=pi/2)
fig[4,2] = Label(fig, "no detection", tellheight=false, rotation=pi/2)

label_A = fig[3,3,TopLeft()] = Label(fig, "A", font=noto_sans_bold, textsize=56, halign=:left, padding=(0,0,0,0))
label_B = fig[3,4,TopLeft()] = Label(fig, "B", font=noto_sans_bold, textsize=56, halign=:left, padding=(0,0,0,0))
label_C = fig[3,5,TopLeft()] = Label(fig, "C", font=noto_sans_bold, textsize=56, halign=:left, padding=(0,0,0,0))
label_D = fig[4,3,TopLeft()] = Label(fig, "D", font=noto_sans_bold, textsize=56, halign=:left)
label_E = fig[4,4,TopLeft()] = Label(fig, "E", font=noto_sans_bold, textsize=56, halign=:left)
label_F = fig[4,5,TopLeft()] = Label(fig, "F", font=noto_sans_bold, textsize=56, halign=:left)

save(joinpath(save_dir, "fig4_detection_traces.$(ext)"), fig)

fig

end