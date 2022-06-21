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

min_reviewers_per_seizure = 3
all_motif_results_df = load(datadir("motif_results_df_reviewers$(min_reviewers_per_seizure).jld2"))["motif_results_df"]
snippets_duration_s = 1
alerts_grace_s = 60
contributions_spec = "tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(snippets_duration_s)_lagextents8x25_helsinkiEEG"
save_dir = plotsdir("motif0_$(contributions_spec)_reviewers$(min_reviewers_per_seizure)_$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, [1:15..., 19,31,44,47,50,62]) do PAT
    rolling_window_len = 60
    min_dist_to_seizure = 30

    target_match_str = "$(contributions_spec)$(PAT)_"

    jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    contributions = jld_dict["contributions"]

    # seizure_bounds, consensus = load_helsinki_seizure_annotations(PAT; min_reviewers_per_seizure=min_reviewers_per_seizure)

    @info "Loading EEG $(PAT)..."
    eeg = load_helsinki_eeg(PAT; min_reviewers_per_seizure = min_reviewers_per_seizure)
    signal_times = get_times(eeg, sample_rate=snippets_duration_s)
    @info "done."

    @info "Calculating signal traces..."
    θ = 3
    motif_results_df = filter(:patient => p -> p == PAT, all_motif_results_df)
    max_significance_μ_weights = motif0_weight_based_effect(motif_results_df, :Δμ, 5)
    detections_mean = apply_rolling_deviation_window(contributions, mean, rolling_window_len) * max_significance_μ_weights

    max_significance_σ_weights = motif0_weight_based_effect(motif_results_df, :Δσ, 5)
    detections_std = apply_rolling_deviation_window(contributions, std, rolling_window_len) * max_significance_μ_weights

    @info "done. Plotting..."



    fig = Figure()
    ax_mean = Axis(fig[1,1:2]; ylabel = "mean μ (five most significant motifs)")
    ax_std = Axis(fig[2,1:2]; ylabel = "mean σ (five most significant motifs)")
    ax_rev = Axis(fig[3,1:2]; ylabel = "# reviewers", xlabel = "time")
    TriCorrApplications.plot_reviewer_consensus!(ax_rev, eeg)
    TriCorrApplications.plot_contribution!(ax_mean, eeg, signal_times, detections_mean)
    TriCorrApplications.plot_contribution!(ax_std, eeg, signal_times, detections_std)
    hlines!(ax_mean, [θ], color=:red, linestyle=:dash)

    fig[:,end+1] = roc_column = GridLayout()
    # save(joinpath(save_dir, "signals_patient$(PAT)_reviewers$(min_reviewers_per_seizure).png"), fig)


    @info "Calculating ROC curve..."
    signals = Dict(
        "σ" => (:Δσ, std),
        "μ" => (:Δμ, mean)
    )
    roc_drws = map(enumerate(["μ", "σ"])) do (i, signal_type)
        roc_data = calculate_patient_ROC(PAT, all_motif_results_df, signals[signal_type]..., rolling_window_len; alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, n_most_significant_motifs=5, n_θs=100, signal=contributions, signal_times=signal_times, min_reviewers_per_seizure=min_reviewers_per_seizure)
        @info "done. Now plotting."

        @show roc_data
        if isnothing(roc_data)
            @show "No ROC for Patient $PAT"
        else
            roc_plt = data(roc_data) * mapping((:false_positives,:gt_negative) => ((f, gt) -> f / gt) => "FPR", (:true_positives,:gt_positive)=> ((t, gt) -> t / gt) => "TPR") * visual(Lines, color=:blue, linewidth=5)
            roc_drw = draw!(roc_column[i,1], roc_plt, axis=(title="Patient $PAT ($(signal_type) signal)", limits=((0.,1.),(0.,1.))))

            roc_drw
        end
    end
    save(joinpath(save_dir, "roc_patient$(PAT)_signals_reviewers$(min_reviewers_per_seizure).png"), fig)
end