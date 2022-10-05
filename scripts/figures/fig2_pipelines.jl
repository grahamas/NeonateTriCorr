using DrWatson
@quickactivate "NeonateTriCorr"

## Saves contributions timeseries to datadir()/exp_pro/
## Excludes seconds marked as artifactual by

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity
using FFTW

include(scriptsdir("include_src.jl"))

noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")
font_theme = Theme(fontsize=36, Lines=Theme(color=:black), linewidth=4)
set_theme!(font_theme)

let signal_types = ["aEEG", "tricorr"],
    aEEG_signals_reduction_name = "meanall",
    tricorr_signals_reduction_name = "$(aEEG_signals_reduction_name)abs",
    patients_considered = [13],
    standardization="within",
    resolution=(4000,2000);

save_dir = plotsdir("contributions_plots_std$(standardization)_$(Dates.now())")
mkpath(save_dir)

signals_reduction_names = Dict(
    "aEEG" => aEEG_signals_reduction_name,
    "tricorr" => tricorr_signals_reduction_name
)

params = Dict(st => merge(common_params, analysis_particular_params[st], Dict(:signals_reduction_name => signals_reduction_names[st], :standardization => standardization)) for st ∈ signal_types)

epoch_s = first(params)[2][:epoch_s]
@assert all(params[st][:epoch_s] == epoch_s for st ∈ signal_types)

for patient_num ∈ patients_considered
    fig = Figure(resolution=resolution, axis=(labelsize=36,))
    n_st = (length(signal_types))
    cons_ax = Axis(fig)
    for (i_st, st) ∈ enumerate(signal_types)
        p = params[st]
        eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure=p[:min_reviewers_per_seizure], excluded_artifact_grades=Int[1], discretization_s=p[:discretization_s])
        signal = load_signal(; p..., 
            signal_type=st, patient_num=patient_num
        )
        times = get_times(eeg, sample_rate=1/p[:snippets_duration_s])
        rolled_signal = roll_signals(signal; p...)
        not_missings = .!ismissing.(rolled_signal[1,:])
        rolled_signal .-= mean(rolled_signal[:, not_missings], dims=:time)
        rolled_signal ./= std(rolled_signal[:, not_missings], dims=:time)
        fig[1:10,i_st] = signal_plts = plot_contributions!(fig, eeg, times, rolled_signal; get_label=st, epoch_s=p[:epoch_s], consensus_plot=true)


    end

    label_A = fig[1,1,TopLeft()] = Label(fig, "A", font=noto_sans_bold, textsize=56, halign=:left)

    
    save(joinpath(save_dir, "both_contributions_patient$(patient_num)_$(Dates.now()).png"), fig)

    fig
end

end