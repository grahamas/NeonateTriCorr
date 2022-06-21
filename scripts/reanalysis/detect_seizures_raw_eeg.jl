using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

min_reviewers_per_seizure = 3
all_channel_results_df = load(datadir("channel_results_df_reviewers$(min_reviewers_per_seizure).jld2"))["channel_results_df"]
snippets_duration_s = 1
alert_grace_s=60

save_dir = plotsdir("raw_eeg_1Hz_reviewers$(min_reviewers_per_seizure)_$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, [47]) do PAT#[1:15..., 19,31,44,47,50,62]) do PAT
    rolling_window_len = 60
    min_dist_to_seizure = 30

    @info "Loading EEG $(PAT)..."
    eeg = load_helsinki_eeg(PAT; min_reviewers_per_seizure = min_reviewers_per_seizure)
    raw_signal = get_signal(eeg)
    channels = 1:size(raw_signal,1)

    # Get signal in 1 Hz
    signal = mapreduce(hcat, eeg.sample_rate:eeg.sample_rate:size(raw_signal,2)) do i
        mean(raw_signal[:,(i-eeg.sample_rate+1):i]; dims=2)
    end

    @info "done."

    @info "Calculating signal traces..."
    θ = 3
    channel_results_df = filter(:patient => p -> p == PAT, all_channel_results_df)
    max_significance_μ_weights = motif_weights_based_on_pvalues(channel_results_df, :Δμ, 5)
    detections_mean = apply_rolling_deviation_window(signal, mean, rolling_window_len) * max_significance_μ_weights

    max_significance_σ_weights = motif_weights_based_on_pvalues(channel_results_df, :Δσ, 5)
    detections_std = apply_rolling_deviation_window(signal, std, rolling_window_len) * max_significance_μ_weights

    @info "done. Plotting..."



    fig = Figure()
    ax_mean = Axis(fig[1,1:2]; ylabel = "mean μ (five most significant channels)")
    ax_std = Axis(fig[2,1:2]; ylabel = "mean σ (five most significant channels)")
    ax_rev = Axis(fig[3,1:2]; ylabel = "# reviewers", xlabel = "time")
    TriCorrApplications.plot_reviewer_consensus!(ax_rev, eeg)
    TriCorrApplications.plot_contribution!(ax_mean, eeg, get_times(eeg, sample_rate=1), detections_mean)
    TriCorrApplications.plot_contribution!(ax_std, eeg, get_times(eeg, sample_rate=1), detections_std)
    hlines!(ax_mean, [θ], color=:red, linestyle=:dash)

    fig[:,end+1] = roc_column = GridLayout()
    # save(joinpath(save_dir, "signals_patient$(PAT)_reviewers$(min_reviewers_per_seizure).png"), fig)


    @info "Calculating ROC curve..."
    signals = Dict(
        "σ" => (:Δσ, std),
        "μ" => (:Δμ, mean)
    )
    roc_drws = map(enumerate(["μ", "σ"])) do (i, signal_type)
        roc_data = calculate_patient_ROC(PAT, all_channel_results_df, signals[signal_type]..., rolling_window_len; alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, n_most_significant_motifs=5, n_θs=100, signal=signal, signal_times=signal_times, min_reviewers_per_seizure=min_reviewers_per_seizure)
        @info "done. Now plotting."
        
        if !isnothing(roc_data)
            roc_plt = data(roc_data) * mapping((:false_positives,:gt_negative) => ((f, gt) -> f / gt) => "FPR", (:true_positives,:gt_positive)=> ((t, gt) -> t / gt) => "TPR") * visual(Lines, color=:blue, linewidth=5)
            roc_drw = draw!(roc_column[i,1], roc_plt, axis=(title="Patient $PAT ($(signal_type) signal)", limits=((0.,1.),(0.,1.))))

            roc_drw
        end
    end
    save(joinpath(save_dir, "roc_patient$(PAT)_signals_reviewers$(min_reviewers_per_seizure).png"), fig)
end