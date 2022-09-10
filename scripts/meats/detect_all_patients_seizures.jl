function detect_all_patients_seizures(patients_considered; save_dir,
    excluded_artifact_grades, min_reviewers_per_seizure, snippets_duration_s,
    task_name, signal_from_dct_fn, reduce_signals_fn, signals_reduction_params, remaining_params...
)
    
    # Step 1: Standardize signals across all patients

    signals = map(patients_considered) do patient_num
        target_match_str = make_signal_stem("tricorr"; 
            excluded_artifact_grades=excluded_artifact_grades,
            min_reviewers_per_seizure=min_reviewers_per_seizure,
            snippets_duration_s=snippets_duration_s,
            patient_num=patient_num,
            remaining_params...
        )
        maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
        if isnothing(maybe_dict)
            @error "No data like $target_match_str"
        else
            signal_from_dct_fn(maybe_dict)
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

    eegs = map(patients_considered) do patient_num
        load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades)
    end
    signal_times = map(eegs) do eeg
        get_times(eeg, sample_rate=1/snippets_duration_s)
    end
    seizure_bounds = map(patients_considered) do patient_num
        load_helsinki_seizure_annotations(patient_num; 
            min_reviewers_per_seizure=min_reviewers_per_seizure
        )
    end

    rolled_signals = map(signals) do signal
        reduce_signals_fn(signal; snippets_duration_s, signals_reduction_params...)
    end

    min_θ, max_θ = calculate_threshold_range(rolled_signal, signal_times,truth_bounds; alert_grace_s=alert_grace_s)

    


    # Step 3: Calculate ROC for every standardized patient

    for patient_num ∈ patients_considered
        eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades)
        signals = standardized_signals[:,:,patient_num] FIXME
        signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

        seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

        fig = plot_μ_and_σ_signals_and_roc(signals, signal_times, seizure_bounds; analysis_eeg=eeg, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)", remaining_params...)

        save(joinpath(save_dir, "$(task_name)_roc_patient$(patient_num)_reviewers$(min_reviewers_per_seizure).png"), fig)
    end
end