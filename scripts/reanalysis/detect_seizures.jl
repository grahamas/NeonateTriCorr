using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

## Saves contributions timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

all_motif_results_df = load(datadir("motif_results_df.jld2"))["motif_results_df"]
snippets_duration = 1
save_dir = plotsdir("tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(snippets_duration)_lagextents8x25_helsinkiEEG_$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, [1:15..., 19,31,44,47,50,62]) do PAT
    rolling_window_len = 60
    min_dist_to_seizure = 30
    min_reviewers_per_seizure = 2

    target_match_str = "tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(snippets_duration)_lagextents8x25_helsinkiEEG$(PAT)_"

    jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    contributions = jld_dict["contributions"]

    # seizure_bounds, consensus = load_helsinki_seizure_annotations(PAT; min_reviewers_per_seizure=min_reviewers_per_seizure)

    @info "Loading EEG $(PAT)..."
    eeg = load_helsinki_eeg(PAT; min_reviewers_per_seizure = min_reviewers_per_seizure)
    @info "done."

    @info "Calculating signal traces..."
    θ = 3
    motif_results_df = filter(:patient => p -> p == PAT, all_motif_results_df)
    max_significance_μ_weights = motif_weights_based_on_pvalues(motif_results_df, :Δμ, 5)
    detections_mean = apply_rolling_deviation_window(contributions, mean, rolling_window_len) * max_significance_μ_weights
    @info "done. Plotting..."

    fig_mean = Figure()
    ax_mean = Axis(fig_mean[1,1])
    ax_mean_rev = Axis(fig_mean[2,1])
    TriCorrApplications.plot_reviewer_consensus!(ax_mean_rev, eeg)
    TriCorrApplications.plot_contribution!(ax_mean, eeg, get_times(eeg, sample_rate=1), detections_mean)
    hlines!(ax_mean, [θ], color=:red, linestyle=:dash)
    save(joinpath(save_dir, "signal_mean_patient$(PAT)_nrev$(min_reviewers_per_seizure).png"), fig_mean)
    @info "done. Now σ..."

    θ = 3
    max_significance_σ_weights = motif_weights_based_on_pvalues(motif_results_df, :Δσ, 5)
    detections_std = apply_rolling_deviation_window(contributions, std, rolling_window_len) * max_significance_μ_weights

    fig_std = Figure()
    ax_std = Axis(fig_std[1,1])
    ax_std_rev = Axis(fig_std[2,1])
    TriCorrApplications.plot_reviewer_consensus!(ax_mean_rev, eeg)
    TriCorrApplications.plot_contribution!(ax_std, eeg, get_times(eeg, sample_rate=1), detections_std)
    hlines!(ax_std, [θ], color=:red, linestyle=:dash)
    fig_std
    save(joinpath(save_dir, "signal_std_patient$(PAT)_nrev$(min_reviewers_per_seizure).png"), fig_std)
    @info "done."


    @info "Calculating ROC curve..."
    signals = Dict(
        "σ" => (:Δσ, std),
        "μ" => (:Δμ, mean)
    )
    roc_drws = map(["μ", "σ"]) do signal_type
        roc_data = calculate_patient_ROC(PAT, all_motif_results_df, signals[signal_type]..., rolling_window_len; alert_grace=60, n_most_significant_motifs=5, n_θs=100, target_match_str=target_match_str, min_reviewers_per_seizure=min_reviewers_per_seizure)
        @info "done. Now plotting."
        
        if !isnothing(roc_data)
            roc_plt = data(roc_data) * mapping((:false_positives,:gt_negative) => ((f, gt) -> f / gt) => "FPR", (:true_positives,:gt_positive)=> ((t, gt) -> t / gt) => "TPR") * visual(Lines, color=:blue, linewidth=5)
            roc_drw = draw(roc_plt, axis=(title="Patient $PAT ($(signal_type) signal)", limits=((0.,1.),(0.,1.))))

            save(joinpath(save_dir, "roc_patient$(PAT)_signal$(signal_type)_reviewers$(min_reviewers_per_seizure).png"), roc_drw)

            roc_drw
        end
    end
end