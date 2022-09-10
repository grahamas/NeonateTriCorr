
function calculate_epoch_differences_across_patients(signal_type, patients; get_channel_label=nothing, session_id=Dates.now(), params...)

    save_root = make_epochdiff_stem(signal_type; params...)
    save_dir = plotsdir("$(save_root)$(session_id)")

    results_df = mapreduce(vcat, patients) do patient_num
        target_match_str = make_signal_stem(signal_type; params..., patient_num=patient_num)
        maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
        signal = if isnothing(maybe_dict)
            @error "No signal named $(target_match_str)"
        else
            if signal_type ∈ keys(maybe_dict)
                if signal_type == "aEEG"
                    aEEG_lower_margin(maybe_dict[signal_type])
                else
                    maybe_dict[signal_type]
                end
            elseif signal_type == "tricorr"
                maybe_dict["contributions"]
            else
                @error "Unsupported signal type: $signal_type"
            end
        end
        if get_channel_label |> isnothing
            calculate_epoch_differences_single_patient(signal;
                save_dir = save_dir,
                params..., patient_num=patient_num
            )
        else
            calculate_epoch_differences_single_patient(signal;
                save_dir = save_dir, get_channel_label = get_channel_label,
                params..., patient_num=patient_num
            )
        end
    end



    @assert !isempty(results_df)
    results_df = DataFrame(results_df)
    all_channel_labels = if isnothing(get_channel_label)
        unique(results_df.channel)
    else
        1:length(unique(results_df.channel)) .|> get_channel_label
    end
    fig_significance = draw_significances_plot!(results_df; all_channel_labels=all_channel_labels)
    fig_Δμ = draw_Δμ_plot(results_df; all_channel_labels=all_channel_labels)
    fig_Δσ = draw_Δσ_plot(results_df; all_channel_labels=all_channel_labels)

    @save datadir("$(save_root)$(session_id).jld2") results_df

    return (results_df, fig_significance, fig_Δμ, fig_Δσ)
end


function calculate_epoch_differences_single_patient(signal::NamedDimsArray;
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
        rolling_window_s,
        get_channel_label = i -> eeg.labels[i],
        unused_params...
    )
    signal = NamedDimsArray{(:_, :time)}(signal)
    n_channels = size(signal, 1)
    mkpath(save_dir)

    @show "snippets duration: $snippets_duration_s"
    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    rolling_window = rolling_window_s ÷ snippets_duration_s
    rolling_std_z = mapreduce(hcat, 1:n_channels) do i_signal
        stds = [std(signal[i_signal,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,:time))]
        (stds .- mean(skipmissing(stds))) ./ std(skipmissing(stds))
    end
    rolling_mean_z = mapreduce(hcat, 1:n_channels) do i_signal
        means = [mean(signal[i_signal,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,:time))]
        (means .- mean(skipmissing(means))) ./ std(skipmissing(means))
    end

    fig_std = plot_contributions(eeg, signal_times[rolling_window:end], rolling_std_z'; title="Rolling σ (zscored per channel; window = $(rolling_window))", resolution=(800,n_channels*100), get_label=get_channel_label)
    save(joinpath(save_dir, "standard_deviations_zscored_window$(rolling_window)_pat$(patient_num).png"), fig_std)

    fig_mean = plot_contributions(eeg, signal_times[rolling_window:end], rolling_mean_z'; title="Rolling μ (zscored per channel; window = $(rolling_window))", resolution=(800,n_channels*100), get_label=get_channel_label)
    save(joinpath(save_dir, "means_zscored_window$(rolling_window)_pat$(patient_num).png"), fig_mean)

    # SAME AS ABOVE, but not zscored
    rolling_std = mapreduce(hcat, 1:n_channels) do i_signal
        [std(signal[i_signal,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,:time))]
    end
    rolling_mean = mapreduce(hcat, 1:n_channels) do i_signal
        [mean(signal[i_signal,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,:time))]
    end

    fig_std = plot_contributions(eeg, signal_times[rolling_window:end], rolling_std'; title="Rolling σ (window = $(rolling_window))", resolution=(800,n_channels*100), get_label=get_channel_label)
    save(joinpath(save_dir, "standard_deviations_window$(rolling_window)_pat$(patient_num).png"), fig_std)

    fig_mean = plot_contributions(eeg, signal_times[rolling_window:end], rolling_mean'; title="Rolling μ (window = $(rolling_window))", resolution=(800,n_channels*100), get_label=get_channel_label)
    save(joinpath(save_dir, "means_window$(rolling_window)_pat$(patient_num).png"), fig_mean)

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

    control_signal = signal[:, control_snippets]
    seizure_signal = signal[:, seizure_snippets]

    if sum(.!ismissing.(seizure_signal)) < (n_channels * min_snippets_for_comparison)
        @warn "Not enough non-missing seizure snippets ($(length(seizure_snippets))) in Patient $patient_num"
        return []
    end

    channel_JSDistances_KSpvalues_effect = map(1:n_channels) do i_signal
        channel_control_signal = control_signal[i_signal, :] |> skipmissing |> collect
        channel_seizure_signal = seizure_signal[i_signal, :] |> skipmissing |> collect

        (JS=estimate_JSDistance(channel_control_signal, channel_seizure_signal; show_plots=false),
        p=pvalue(ApproximateTwoSampleKSTest(channel_control_signal, channel_seizure_signal)),
        Δμ=mean(channel_seizure_signal) - mean(channel_control_signal),
        Δσ=std(channel_seizure_signal) - std(channel_control_signal),
        patient=patient_num,
        channel=get_channel_label(i_signal))
    end # map(1:n_channels)

    all_channel_labels = 1:n_channels .|> get_channel_label

    fig_kdes = plot_estimated_distributions(control_signal, seizure_signal; title="Patient = $(patient_num)", all_channel_labels=all_channel_labels)
    save(joinpath(save_dir, "estimated_distributions_pat$(patient_num).png"), fig_kdes)


    return (channel_JSDistances_KSpvalues_effect)

end
