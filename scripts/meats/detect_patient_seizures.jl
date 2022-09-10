function detect_patient_seizures(patient_num; save_dir,
    excluded_artifact_grades, min_reviewers_per_seizure, snippets_duration_s,
    task_name, signal_from_dct_fn, remaining_params...
)
    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=excluded_artifact_grades)

    target_match_str = make_signal_stem("tricorr"; 
        excluded_artifact_grades=excluded_artifact_grades,
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        snippets_duration_s=snippets_duration_s,
        remaining_params...
    )
    maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    signals = if isnothing(maybe_dict)
        @error "No data like $target_match_str"
    else
        signal_from_dct_fn(maybe_dict)
    end
    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    fig = plot_μ_and_σ_signals_and_roc(signals, signal_times, seizure_bounds; analysis_eeg=eeg, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)", remaining_params...)

    save(joinpath(save_dir, "$(task_name)_roc_patient$(patient_num)_reviewers$(min_reviewers_per_seizure).png"), fig)
end