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
        snippets_duration_s, signal_from_dct_fn = get_signal_from_dct_fn(signal_type), 
        signals_reduction_name,
        epoch_s, n_θs, 
        discretization_s,
        task_name = "$(signal_type)$(signals_reduction_name)",
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
        sig = signal_from_dct_fn(maybe_dict)
        discretize_missings!(sig, discretization_s ÷ snippets_duration_s)
        sig
    end

    # need to validate that all signals are in same order

    eegs = map(patients_considered) do patient_num
        load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades, discretization_s=discretization_s)
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

    # Step 2: Calculate ROC across all patients

    signal_times = map(eegs) do eeg
        get_times(eeg, sample_rate=1/snippets_duration_s)
    end
    seizure_bounds = map(patients_considered) do patient_num
        bounds, _ = load_helsinki_seizure_annotations(patient_num; 
            min_reviewers_per_seizure=min_reviewers_per_seizure, discretization_s=discretization_s
        )
        bounds
    end

    rolled_signals = roll_signals.(signals; snippets_duration_s=snippets_duration_s, remaining_params...)

    cat_signals = cat(rolled_signals..., dims=:time)
    cat_not_missings = .!ismissing.(cat_signals[1,:])
    signal_means = mean(cat_signals[:,cat_not_missings], dims=:time)
    signal_stds = std(cat_signals[:,cat_not_missings], dims=:time)

    for idx ∈ eachindex(rolled_signals)
        rolled_signals[idx] .-= signal_means
        rolled_signals[idx] ./= signal_stds
    end

    reduced_signals = reduce_signals_fn.(rolled_signals)

    non_seizure_hours = mapreduce(+, seizure_bounds, signal_times) do bounds, times
        calculate_negepoch_non_seizure_hours(bounds, times, epoch_s, snippets_duration_s)
    end

    min_θ, max_θ = calculate_threshold_range(reduced_signals)

    function detect_all_patients(θ)
        roc_vals = nothing
        clean_negatives = (
            gt_negative_clean=0,
            false_positives_clean=0
        )
        for (signal, bounds, times) in zip(reduced_signals, seizure_bounds, signal_times)
            single_patient_roc_vals = evaluate_detection_posseizure_negepoch(bounds, signal .>= θ, times; snippets_duration_s=snippets_duration_s, epoch_s=epoch_s)
            if isempty(bounds)
                clean_negatives = add_nts(clean_negatives, (
                    gt_negative_clean=single_patient_roc_vals.gt_negative,
                    false_positives_clean=single_patient_roc_vals.false_positives
                ))
            end
            roc_vals = add_nts(roc_vals, single_patient_roc_vals)
        end
        non_seizure_hours = mapreduce(+, seizure_bounds, signal_times) do bounds, times
            calculate_negepoch_non_seizure_hours(bounds, times, epoch_s, snippets_duration_s)
        end
        return merge(roc_vals, clean_negatives)
    end

    all_patient_ROC_df = calculate_ROC(detect_all_patients, range(min_θ, max_θ, length=n_θs))
    save_multipatient_ROC_df(patients_considered, all_patient_ROC_df; 
        signal_type=signal_type, 
        excluded_artifact_grades=excluded_artifact_grades, 
        min_reviewers_per_seizure=min_reviewers_per_seizure, 
        snippets_duration_s=snippets_duration_s, task_name=task_name, 
        signal_from_dct_fn=signal_from_dct_fn, 
        signals_reduction_name=signals_reduction_name, 
        non_seizure_hours=non_seizure_hours,
        epoch_s=epoch_s, n_θs=n_θs, remaining_params...
    )

    auc = calculate_AUC(all_patient_ROC_df)
    all_patient_fig = Figure(resolution=roc_resolution)
    plot_seizure_detection_ROC!(all_patient_fig[1,1], all_patient_ROC_df; epoch_s=epoch_s, non_seizure_hours=non_seizure_hours, title="Patients $(patients_considered) (AUC = $auc)")
    save(joinpath(save_dir, "$(task_name)_$(length(patients_considered))patients_roc_nrev$(min_reviewers_per_seizure).png"), all_patient_fig)

    # # Where the FPR is only from non-seizing patients
    # n_clean_patients = sum(isempty.(seizure_bounds))
    # clean_patient_ROC_df = copy(all_patient_ROC_df)
    # clean_patient_ROC_df.gt_negative = clean_patient_ROC_df.gt_negative_clean
    # clean_patient_ROC_df.false_positives = clean_patient_ROC_df.false_positives_clean
    # clean_auc = calculate_AUC(clean_patient_ROC_df)
    # clean_patient_fig = Figure(resolution=roc_resolution)
    # plot_seizure_detection_ROC!(clean_patient_fig[1,1], clean_patient_ROC_df; non_seizure_hours=non_seizure_hours_clean, title="Patients $(patients_considered); FPR from $(n_clean_patients) clean patients (AUC = $clean_auc)")
    # save(joinpath(save_dir, "$(task_name)_$(length(patients_considered))_FPRclean_patients_roc_nrev$(min_reviewers_per_seizure).png"), clean_patient_fig)

    # # Step 3: Calculate and plot ROC for every standardized patient

    for (patient_num, eeg, signal, times, bounds) ∈ zip(patients_considered, eegs, signals, signal_times, seizure_bounds)
        fig = plot_μ_and_σ_signals_and_roc(signal, times, bounds; analysis_eeg=eeg, epoch_s=epoch_s, snippets_duration_s=snippets_duration_s, 
        signals_reduction_name=signals_reduction_name, title="Patient $(patient_num) (within standardized)", remaining_params...)

        save(joinpath(save_dir, "$(task_name)_roc_patient$(patient_num)_nrev$(min_reviewers_per_seizure).png"), fig)
    end

    return all_patient_ROC_df
end

function evaluate_clinician_FPR(patients_considered; 
    excluded_artifact_grades,
    epoch_s, remaining_params...
)

    non_seizure_hours = 0
    eegs_and_results = map(patients_considered) do patient_num
        eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = 3, excluded_artifact_grades=excluded_artifact_grades)
        truth_bounds, _ = load_helsinki_seizure_annotations(patient_num; 
            min_reviewers_per_seizure=3, discretization_s=discretization_s
        )
        any_clinician_bounds, _ = load_helsinki_seizure_annotations(patient_num; 
            min_reviewers_per_seizure=1, discretization_s=discretization_s
        )
        times = get_times(eeg, sample_rate=1.)
        detections = zeros(Bool, size(times))
        for (sr,sp) ∈ any_clinician_bounds
            detections[sr:sp] .= true
        end
        results = evaluate_detection_posseizure_negepoch(truth_bounds, detections, times; snippets_duration_s=1, epoch_s=epoch_s)
        return (eeg, results)
    end

    eegs = first.(eegs_and_results)
    all_patient_results = reduce(add_nts, last.(eegs_and_results))

    return (eegs, all_patient_results)

end
