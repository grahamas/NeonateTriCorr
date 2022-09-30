
function load_or_calculate_epoch_differences_across_patients(signal_type, patients; session_id=Dates.now(), params...)
    save_root = make_epochdiff_stem(signal_type; params...)
    maybe_jld = load_most_recent_jld2(save_root, datadir())
    results_df = if isnothing(maybe_jld)
        results_df = calculate_epoch_differences_across_patients(signal_type, patients; session_id=session_id, params...)
        @save datadir("$(save_root)$(session_id).jld2") results_df
        results_df
    else
        maybe_jld["results_df"]
    end
end


function calculate_epoch_differences_across_patients(signal_type, patients; get_channel_label=nothing, session_id=Dates.now(), params...)

    save_root = make_epochdiff_stem(signal_type; params...)
    save_dir = plotsdir("$(save_root)$(session_id)")

    results_df = mapreduce(vcat, patients) do patient_num
        target_match_str = make_signal_stem(signal_type; params..., patient_num=patient_num)
        maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
        signal = if isnothing(maybe_dict)
            @error "No signal named $(target_match_str)"
        else
            get_signal_from_dct_fn(signal_type)(maybe_dict)
        end
        if get_channel_label |> isnothing
            calculate_epoch_differences_single_patient(signal;
                save_dir = save_dir,
                params..., patient_num=patient_num, signal_type=signal_type
            )
        else
            calculate_epoch_differences_single_patient(signal;
                save_dir = save_dir, get_channel_label = get_channel_label,
                params..., patient_num=patient_num, signal_type=signal_type
            )
        end
    end
    return results_df
end

function plot_epoch_differences_across_patients(signal_type, patients; session_id=Dates.now(), get_channel_label=nothing, patient_res_width=20, params...)
    results_df = load_or_calculate_epoch_differences_across_patients(signal_type, patients; session_id=session_id, get_channel_label=get_channel_label, params...)

    @assert !isempty(results_df)
    results_df = DataFrame(results_df)
    all_channel_labels = if isnothing(get_channel_label)
        unique(results_df.channel)
    else
        1:length(unique(results_df.channel)) .|> get_channel_label
    end
    fig_p = (resolution=(2600,2100),)
    axis_p = (xticklabelrotation=pi/2,)
    fig_significance = draw_significances_plot!(results_df; all_channel_labels=all_channel_labels, figure=fig_p, axis=axis_p)
    fig_Δμ = draw_Δμ_plot(results_df; all_channel_labels=all_channel_labels, figure=fig_p, axis=axis_p)
    fig_Δσ = draw_Δσ_plot(results_df; all_channel_labels=all_channel_labels, figure=fig_p, axis=axis_p)

    return (results_df, fig_significance, fig_Δμ, fig_Δσ)
end


function calculate_epoch_differences_single_patient(signal::NamedDimsArray;
        save_dir,
        excluded_artifact_grades,
        patient_num,
        discretization_s,
        eeg = load_helsinki_eeg(patient_num; 
            excluded_artifact_grades=excluded_artifact_grades,
            discretization_s=discretization_s
        ),
        min_reviewers_per_seizure,
        min_snippets_for_comparison,
        snippets_duration_s,
        get_channel_label = i -> eeg.labels[i],
        params...
    )
    signal = NamedDimsArray{(:_, :time)}(signal)
    n_channels = size(signal, 1)
    mkpath(save_dir)

    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure, discretization_s=discretization_s)

    rolled_signal = roll_signals(signal; snippets_duration_s=snippets_duration_s, params...)

    standardized_signal = standardize_signals!(rolled_signal; params...)

    fig = plot_contributions(eeg, signal_times, standardized_signal; title="Rolling $(string(params[:window_fn])) (zscored $(params[:standardization]); window = $(params[:rolling_window_s])s)", resolution=(800,n_channels*100), get_label=params[:signal_type])
    save(joinpath(save_dir, "rolling_$(string(params[:window_fn]))_window$(params[:rolling_window_s])_std$(params[:standardization])_pat$(patient_num).png"), fig)

    if isempty(seizure_bounds)
        @warn "Patient $(patient_num): no seizures."
        return []
    end

    epochized_seizure_bounds = epochize_bounds(seizure_bounds, first(signal_times), last(signal_times); params...)
    epochized_artifact_bounds = epochize_bounds(eeg.artifact_annotations, first(signal_times), last(signal_times); params...)
    noncontrol_bounds = merge_bounds(epochized_seizure_bounds, epochized_artifact_bounds)
    control_bounds = invert_bounds(noncontrol_bounds, 0, eeg.duration)

    if isempty(control_bounds)
        @warn "Patient $(patient_num): no control snippets"
        return []
    end

    control_snippets = mapreduce((x,y) -> x .|| y, control_bounds) do (on, off)
        on .<= (signal_times .+ 1) .< off
    end
    seizure_snippets = mapreduce((x,y) -> x .|| y, seizure_bounds) do (on, off)
        on .<= (signal_times .+ 1) .< off
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
