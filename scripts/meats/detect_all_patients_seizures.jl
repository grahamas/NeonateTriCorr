add_nts(nt1::Nothing, nt2::NamedTuple) = nt2
add_nts(nt1::NamedTuple, nt2::Nothing) = nt1
function add_nts(nt1::NT, nt2::NT) where {L, NT <: NamedTuple{L}}
    NamedTuple(
        key => nt1[key] + nt2[key] for key in L
    )
end

function load_signal(; discretization_s, snippets_duration_s, signal_type, signal_from_dct_fn = get_signal_from_dct_fn(signal_type), params...)
    target_match_str = make_signal_stem(signal_type; snippets_duration_s=snippets_duration_s, params...)
    maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    @assert !isnothing(maybe_dict) "No data like $target_match_str"
    sig = signal_from_dct_fn(maybe_dict)
    discretize_missings!(sig, discretization_s ÷ snippets_duration_s)
    sig
end

function detect_all_patients_seizures(patients_considered; signal_type, 
        save_dir,
        excluded_artifact_grades, min_reviewers_per_seizure, 
        snippets_duration_s, signal_from_dct_fn = get_signal_from_dct_fn(signal_type), 
        signals_reduction_name,
        epoch_s, 
        discretization_s,
        task_name = "$(signal_type)$(signals_reduction_name)",
        roc_resolution=(800,600),
        calculate_targets_fn,
        remaining_params...
    )
    reduce_signals_fn = get_reduce_signals_fn(signals_reduction_name)

    # Step 1: Standardize signals across all patients

    signals = map(patients_considered) do patient_num
        load_signal(;
            signal_type = signal_type,
            excluded_artifact_grades=excluded_artifact_grades,
            min_reviewers_per_seizure=min_reviewers_per_seizure,
            snippets_duration_s=snippets_duration_s,
            patient_num=patient_num,
            discretization_s=discretization_s,
            epoch_s=epoch_s,
            remaining_params...
        )
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

    standardize_signals!(rolled_signals; remaining_params...)

    reduced_signals = reduce_signals_fn.(rolled_signals)

    non_seizure_hours = mapreduce(+, seizure_bounds, signal_times) do bounds, times
        calculate_negepoch_non_seizure_hours(bounds, times, epoch_s, snippets_duration_s)
    end

    tars_and_nons = map(reduced_signals, signal_times, seizure_bounds) do signal, times, bounds
         targets, non_targets = calculate_targets_fn(signal, times, bounds; epoch_s=epoch_s, snippets_duration_s=snippets_duration_s)
         (targets, non_targets)
    end
    targets, non_targets = reduce(tars_and_nons) do (tar1, non1), (tar2, non2)
        (vcat(tar1, tar2), vcat(non1, non2))
    end

    r = roc(targets, non_targets)
    save_multipatient_ROC(patients_considered, all_patient_ROC_df; 
        signal_type=signal_type, 
        excluded_artifact_grades=excluded_artifact_grades, 
        min_reviewers_per_seizure=min_reviewers_per_seizure, 
        snippets_duration_s=snippets_duration_s, task_name=task_name, 
        signal_from_dct_fn=signal_from_dct_fn, 
        signals_reduction_name=signals_reduction_name, 
        non_seizure_hours=non_seizure_hours,
        calculate_targets_fn=calculate_targets_fn,
        epoch_s=epoch_s, remaining_params...
    )


    fig, plt, ax = plot_roc_fphr(r; resolution=roc_resolution)
    save(joinpath(save_dir, "$(task_name)_$(length(patients_considered))patients_roc_nrev$(min_reviewers_per_seizure).png"), fig)

    return r
end

function evaluate_clinician_FPR(patients_considered; 
    excluded_artifact_grades,
    epoch_s, discretization_s, remaining_params...
)
    eegs_and_results = map(patients_considered) do patient_num
        eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = 3, excluded_artifact_grades=excluded_artifact_grades, discretization_s=15)
        truth_bounds, _ = load_helsinki_seizure_annotations(patient_num; 
            min_reviewers_per_seizure=3, discretization_s=discretization_s
        )
        any_clinician_bounds, _ = load_helsinki_seizure_annotations(patient_num; 
            min_reviewers_per_seizure=1, discretization_s=discretization_s
        )
        times = get_times(eeg, sample_rate=1.)
        detections = zeros(Bool, size(times))
        for (sr,sp) ∈ any_clinician_bounds
            detections[sr .<= times .< sp] .= true
        end
        results = evaluate_detection_posseizure_negepoch(truth_bounds, detections, times; snippets_duration_s=1, epoch_s=epoch_s)
        return (eeg, results)
    end

    eegs = first.(eegs_and_results)
    all_patient_results = reduce(add_nts, last.(eegs_and_results))

    return (eegs, all_patient_results)

end
