add_nts(nt1::Nothing, nt2::NamedTuple) = nt2
add_nts(nt1::NamedTuple, nt2::Nothing) = nt1
function add_nts(nt1::NT, nt2::NT) where {L, NT <: NamedTuple{L}}
    NamedTuple(
        key => nt1[key] + nt2[key] for key in L
    )
end

function detect_all_patients_seizures(patients_considered; signal_type,
        save_dir,
        excluded_artifact_grades, min_reviewers_per_seizure, 
        snippets_duration_s, task_name, signal_from_dct_fn, 
        signals_reduction_name, signals_reduction_params,
        alert_grace_s, n_θs,
        roc_resolution=(800,600), remaining_params...
    )
    reduce_signals_fn = get_reduce_signals_fn(signals_reduction_name)

    # Step 1: Standardize signals across all patients

    signals = map(patients_considered) do patient_num
        target_match_str = make_signal_stem(signal_type; 
            excluded_artifact_grades=excluded_artifact_grades,
            min_reviewers_per_seizure=min_reviewers_per_seizure,
            snippets_duration_s=snippets_duration_s,
            patient_num=patient_num,
            remaining_params...
        )
        maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
        @assert !isnothing(maybe_dict) "No data like $target_match_str"
        signal_from_dct_fn(maybe_dict)
    end

    # need to validate that all signals are in same order

    eegs = map(patients_considered) do patient_num
        load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades)
    end
    if signal_type == "aEEG"
        channel_labels = eegs[1].labels
        @assert length(channel_labels) == length(unique(channel_labels))
        if !all(Ref(channel_labels) .== (eegs .|> (eeg -> eeg.labels)))
            channel_labels = eegs[1].labels |> sort
            @assert all(Ref(channel_labels |> sort) .== (eegs .|> (eeg -> sort(eeg.labels))))
            for i_signal ∈ eachindex(signals)
                eeg = eegs[i_signal]
                resort_idxs = [findfirst(eeg.labels .== labels) for labels ∈ channel_labels]
                signals[i_signal] = signals[i_signal][resort_idxs, :]
            end
        end
    end

    cat_signals = cat(signals..., dims=:time)
    signal_means = mean(cat_signals, dims=:time) # FIXME doesn't handle missings
    signal_stds = std(cat_signals, dims=:time)

    for idx ∈ eachindex(signals)
        signals[idx] .-= signal_means
        signals[idx] ./= signal_stds
    end

    # Step 2: Calculate ROC across all patients

    signal_times = map(eegs) do eeg
        get_times(eeg, sample_rate=1/snippets_duration_s)
    end
    seizure_bounds = map(patients_considered) do patient_num
        bounds, _ = load_helsinki_seizure_annotations(patient_num; 
            min_reviewers_per_seizure=min_reviewers_per_seizure
        )
        bounds
    end

    rolled_signals = map(signals) do signal
        reduce_signals_fn(signal; snippets_duration_s, signals_reduction_params...)
    end

    min_θ, max_θ = calculate_threshold_range(rolled_signals)
    @show min_θ max_θ

    function detect_all_patients(θ)
        roc_vals = nothing
        clean_negatives = (
            gt_negative_clean=0,
            false_positives_clean=0
        )
        for (signal, bounds, times) in zip(rolled_signals, seizure_bounds, signal_times)
            single_patient_roc_vals = evaluate_detection_posseizure_negalerts(bounds, signal .>= θ, times; snippets_duration_s=snippets_duration_s, alert_grace_s=alert_grace_s)
            if isempty(bounds)
                clean_negatives = add_nts(clean_negatives, (
                    gt_negative_clean=single_patient_roc_vals.gt_negative,
                    false_positives_clean=single_patient_roc_vals.false_positives
                ))
            end
            roc_vals = add_nts(roc_vals, single_patient_roc_vals)
        end
        return merge(roc_vals, clean_negatives)
    end

    non_seizure_hours = mapreduce(+, eegs) do eeg
        seizures_with_grace_period = add_grace_to_truth_bounds(
            eeg.seizure_annotations, 
            get_times(eeg, sample_rate=1/snippets_duration_s), 
            alert_grace_s
        )
        seizure_and_artifact_bounds = merge_bounds(seizures_with_grace_period, eeg.artifact_annotations)
        calc_non_seizure_hours(eeg, seizure_and_artifact_bounds)
    end

    non_seizure_hours_clean = mapreduce(+, eegs) do eeg
        if isempty(eeg.seizure_annotations)
            eeg.duration / (60 * 60)
        else
            0
        end
    end

    all_patient_ROC_df = calculate_ROC(detect_all_patients, range(min_θ, max_θ, length=n_θs))

    auc = calculate_AUC(all_patient_ROC_df)
    all_patient_fig = Figure(resolution=roc_resolution)
    plot_seizure_detection_ROC!(all_patient_fig[1,1], all_patient_ROC_df; non_seizure_hours=non_seizure_hours, title="Patients $(patients_considered) (AUC = $auc)")
    save(joinpath(save_dir, "$(task_name)_$(length(patients_considered))patients_roc_nrev$(min_reviewers_per_seizure).png"), all_patient_fig)

    # Where the FPR is only from non-seizing patients
    n_clean_patients = sum(isempty.(seizure_bounds))
    clean_patient_ROC_df = copy(all_patient_ROC_df)
    clean_patient_ROC_df.gt_negative = clean_patient_ROC_df.gt_negative_clean
    clean_patient_ROC_df.false_positives = clean_patient_ROC_df.false_positives_clean
    clean_auc = calculate_AUC(clean_patient_ROC_df)
    clean_patient_fig = Figure(resolution=roc_resolution)
    plot_seizure_detection_ROC!(clean_patient_fig[1,1], clean_patient_ROC_df; non_seizure_hours=non_seizure_hours_clean, title="Patients $(patients_considered); FPR from $(n_clean_patients) clean patients (AUC = $clean_auc)")
    save(joinpath(save_dir, "$(task_name)_$(length(patients_considered))_FPRclean_patients_roc_nrev$(min_reviewers_per_seizure).png"), clean_patient_fig)

    # # Step 3: Calculate ROC for every standardized patient

    for (patient_num, eeg, signal, times, bounds) ∈ zip(patients_considered, eegs, signals, signal_times, seizure_bounds)
        fig = plot_μ_and_σ_signals_and_roc(signal, times, bounds; analysis_eeg=eeg, alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, 
        signals_reduction_name=signals_reduction_name, signals_reduction_params=signals_reduction_params, title="Patient $(patient_num)", remaining_params...)

        save(joinpath(save_dir, "$(task_name)_roc_patient$(patient_num)_nrev$(min_reviewers_per_seizure).png"), fig)
    end

    return (all_patient_ROC_df, rolled_signals, signals)
end
