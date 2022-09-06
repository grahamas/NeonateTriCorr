function calculate_epoch_differences(target_match_str::AbstractString; parent_session_id=nothing, min_reviewers_per_seizure, kwargs...)
    jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    contributions = jld_dict["contributions"]
    save_dir = if isnothing(parent_session_id)
        plotsdir("$(target_match_str)nrev$(min_reviewers_per_seizure)_$(Dates.now())")
    else
        plotsdir("$(target_match_str)nrev$(min_reviewers_per_seizure)_$(parent_session_id)")

    calculate_epoch_differences(contributions; 
        min_reviewers_per_seizure=min_reviewers_per_seizure, 
        save_dir=save_dir, kwargs...
    )
end

function calculate_epoch_differences(contributions;
        save_dir,
        excluded_artifact_grades = [1],
        patient_num,
        eeg = load_helsinki_eeg(patient_num; 
            excluded_artifact_grades=excluded_artifact_grades
        ),
        min_snippets_for_comparison = 80,
        min_dist_to_seizure = 30,
        snippets_duration_s = 1
    )
    mkpath(save_dir)

    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(PAT; min_reviewers_per_seizure=min_reviewers_per_seizure)

    rolling_window = 60
    rolling_std_z = mapreduce(hcat, 1:14) do motif_num
        stds = [std(contributions[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(contributions,:time))]
        (stds .- mean(skipmissing(stds))) ./ std(skipmissing(stds))
    end
    rolling_mean_z = mapreduce(hcat, 1:14) do motif_num
        means = [mean(contributions[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(contributions,:time))]
        (means .- mean(skipmissing(means))) ./ std(skipmissing(means))
    end

    fig_std = plot_contributions(eeg, get_times(eeg, sample_rate=1)[rolling_window:end], rolling_std_z'; title="Rolling σ (zscored per motif; window = $(rolling_window))", resolution=(800,1400))
    save(joinpath(save_dir, "contributions_standard_deviations_zscored_window$(rolling_window)_pat$(PAT).png"), fig_std)

    fig_mean = plot_contributions(eeg, get_times(eeg, sample_rate=1)[rolling_window:end], rolling_mean_z'; title="Rolling μ (zscored per motif; window = $(rolling_window))", resolution=(800,1400))
    save(joinpath(save_dir, "contributions_means_zscored_window$(rolling_window)_pat$(PAT).png"), fig_mean)

    rolling_window = 60
    rolling_std = mapreduce(hcat, 1:14) do motif_num
        [std(contributions[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(contributions,:time))]
    end
    rolling_mean = mapreduce(hcat, 1:14) do motif_num
        [mean(contributions[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(contributions,:time))]
    end

    fig_std = plot_contributions(eeg, get_times(eeg, sample_rate=1)[rolling_window:end], rolling_std'; title="Rolling σ (window = $(rolling_window))", resolution=(800,1400))
    save(joinpath(save_dir, "contributions_standard_deviations_window$(rolling_window)_pat$(PAT).png"), fig_std)

    fig_mean = plot_contributions(eeg, get_times(eeg, sample_rate=1)[rolling_window:end], rolling_mean'; title="Rolling μ (window = $(rolling_window))", resolution=(800,1400))
    save(joinpath(save_dir, "contributions_means_window$(rolling_window)_pat$(PAT).png"), fig_mean)

    if isempty(seizure_bounds)
        @warn "Patient $(PAT): no seizures."
        return []
    end

    control_bounds = calc_control_bounds(eeg, snippets_duration_s, min_dist_to_seizure)
    control_snippets = mapreduce((x,y) -> x .|| y, control_bounds) do (on, off)
        on .<= (signal_times .+ 1) .< off
    end
    seizure_snippets = mapreduce((x,y) -> x .|| y, seizure_bounds) do (on, off)
        on .<= (signal_times .+ 1) .< off
    end


    if (length(control_snippets) < min_snippets_for_comparison)
        @warn "Not enough control snippets ($(length(control_snippets))) in Patient $PAT"
        return []
    end
    if (length(seizure_snippets) < min_snippets_for_comparison) 
        @warn "Not enough seizure snippets ($(length(seizure_snippets))) in Patient $PAT"
        return []
    end

    control_contributions = contributions[:, control_snippets]
    seizure_contributions = contributions[:, seizure_snippets]

    if sum(.!ismissing.(seizure_contributions)) < (14 * min_snippets_for_comparison)
        @warn "Not enough non-missing seizure snippets ($(length(seizure_snippets))) in Patient $PAT"
        return []
    end

    motif_JSDistances_KSpvalues_effect = map(1:14) do motif_num
        motif_control_contributions = control_contributions[motif_num, :] |> skipmissing |> collect
        motif_seizure_contributions = seizure_contributions[motif_num, :] |> skipmissing |> collect

        (JS=estimate_JSDistance(motif_control_contributions, motif_seizure_contributions; show_plots=false),
        p=pvalue(ApproximateTwoSampleKSTest(motif_control_contributions, motif_seizure_contributions)),
        Δμ=mean(motif_seizure_contributions) - mean(motif_control_contributions),
        Δσ=std(motif_seizure_contributions) - std(motif_control_contributions),
        patient=PAT,
        motif=offset_motif_numeral(motif_num))
    end # map(1:14)

    fig_kdes = plot_estimated_distributions(control_contributions, seizure_contributions; title="Patient = $(PAT)")
    save(joinpath(save_dir, "estimated_distributions_pat$(PAT).png"), fig_kdes)


    motif_JSDistances_KSpvalues_effect

    end # mapreduce(PATs)

    motif_results_df = DataFrame(motif_results)
    fig_significance = draw_significances_plot!(motif_results_df)
    fig_Δμ = draw_Δμ_plot(motif_results_df)
    fig_Δσ = draw_Δσ_plot(motif_results_df)

    @save datadir("motif_results_df_reviewers$(min_reviewers_per_seizure).jld2") motif_results_df
    
    return motif_results_df
end
