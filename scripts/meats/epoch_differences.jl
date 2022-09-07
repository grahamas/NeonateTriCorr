
function calculate_epoch_differences_across_patients(signal_type, patients; session_id=Dates.now(), params...)

    save_root = make_epochdiff_stem(signal_type; params...)
    save_dir = plotsdir("$(save_root)$(session_id)")

    results_df = mapreduce(vcat, patients) do patient_num
        target_match_str = make_signal_stem(signal_type; params..., patient_num=patient_num)
        maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
        contributions = if isnothing(maybe_dict)
            @error "No signal named $(target_match_str)"
        else
            maybe_dict["contributions"]
        end
        calculate_epoch_differences_single_patient(contributions;
            save_dir = save_dir,
            params..., patient_num=patient_num
        )
    end

    @assert !isempty(results_df)
    results_df = DataFrame(results_df)
    fig_significance = draw_significances_plot!(results_df)
    fig_Δμ = draw_Δμ_plot(results_df)
    fig_Δσ = draw_Δσ_plot(results_df)

    @save datadir("$(save_root)$(session_id).jld2") results_df

    return (results_df, fig_significance, fig_Δμ, fig_Δσ)
end


function calculate_epoch_differences_single_patient(contributions;
        save_dir,
        excluded_artifact_grades = [1],
        patient_num,
        eeg = load_helsinki_eeg(patient_num; 
            excluded_artifact_grades=excluded_artifact_grades
        ),
        min_reviewers_per_seizure,
        min_snippets_for_comparison,
        min_dist_to_seizure,
        snippets_duration_s,
        unused_params...
    )
    mkpath(save_dir)

    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

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
    save(joinpath(save_dir, "contributions_standard_deviations_zscored_window$(rolling_window)_pat$(patient_num).png"), fig_std)

    fig_mean = plot_contributions(eeg, get_times(eeg, sample_rate=1)[rolling_window:end], rolling_mean_z'; title="Rolling μ (zscored per motif; window = $(rolling_window))", resolution=(800,1400))
    save(joinpath(save_dir, "contributions_means_zscored_window$(rolling_window)_pat$(patient_num).png"), fig_mean)

    rolling_window = 60
    rolling_std = mapreduce(hcat, 1:14) do motif_num
        [std(contributions[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(contributions,:time))]
    end
    rolling_mean = mapreduce(hcat, 1:14) do motif_num
        [mean(contributions[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(contributions,:time))]
    end

    fig_std = plot_contributions(eeg, get_times(eeg, sample_rate=1)[rolling_window:end], rolling_std'; title="Rolling σ (window = $(rolling_window))", resolution=(800,1400))
    save(joinpath(save_dir, "contributions_standard_deviations_window$(rolling_window)_pat$(patient_num).png"), fig_std)

    fig_mean = plot_contributions(eeg, get_times(eeg, sample_rate=1)[rolling_window:end], rolling_mean'; title="Rolling μ (window = $(rolling_window))", resolution=(800,1400))
    save(joinpath(save_dir, "contributions_means_window$(rolling_window)_pat$(patient_num).png"), fig_mean)

    if isempty(seizure_bounds)
        @warn "Patient $(patient_num): no seizures."
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
        @warn "Not enough control snippets ($(length(control_snippets))) in Patient $patient_num"
        return []
    end
    if (length(seizure_snippets) < min_snippets_for_comparison) 
        @warn "Not enough seizure snippets ($(length(seizure_snippets))) in Patient $patient_num"
        return []
    end

    control_contributions = contributions[:, control_snippets]
    seizure_contributions = contributions[:, seizure_snippets]

    if sum(.!ismissing.(seizure_contributions)) < (14 * min_snippets_for_comparison)
        @warn "Not enough non-missing seizure snippets ($(length(seizure_snippets))) in Patient $patient_num"
        return []
    end

    motif_JSDistances_KSpvalues_effect = map(1:14) do motif_num
        motif_control_contributions = control_contributions[motif_num, :] |> skipmissing |> collect
        motif_seizure_contributions = seizure_contributions[motif_num, :] |> skipmissing |> collect

        (JS=estimate_JSDistance(motif_control_contributions, motif_seizure_contributions; show_plots=false),
        p=pvalue(ApproximateTwoSampleKSTest(motif_control_contributions, motif_seizure_contributions)),
        Δμ=mean(motif_seizure_contributions) - mean(motif_control_contributions),
        Δσ=std(motif_seizure_contributions) - std(motif_control_contributions),
        patient=patient_num,
        motif=offset_motif_numeral(motif_num))
    end # map(1:14)

    fig_kdes = plot_estimated_distributions(control_contributions, seizure_contributions; title="Patient = $(patient_num)")
    save(joinpath(save_dir, "estimated_distributions_pat$(patient_num).png"), fig_kdes)


    return motif_JSDistances_KSpvalues_effect

end
