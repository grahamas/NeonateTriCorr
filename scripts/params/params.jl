
patients_all = 1:79
patients_artifact_annotated = [1:15..., 19,31,44,47,50,62,75]
patients_unannotated = setdiff(patients_all, patients_artifact_annotated) 

common_params = Dict(
    :lag_extents => (8,25),
    :min_reviewers_per_seizure => 3,
    :excluded_artifact_grades => Int[],
    :epoch_s => 60,
    :rolling_window_s => 60,
    :window_fn => mean,
    :n_θs => 100,
    :discretization_s => 15,
    :standardization => "within",
    :min_snippets_for_comparison => 20,
    :calculate_targets_fn => pospatient_negepoch
)

analysis_particular_params = Dict(
    "tricorr" => Dict(
        :snippets_duration_s => 1,
        :preproc! => TripleCorrelations.zscore!, 
        :postproc! => TripleCorrelations.identity!,
        :assumption => IndStdNormal(), :conditioned_on => None(),
        :signal_type => "tricorr"
    ),
    "aEEG" => Dict(
        :lowpass_freq => 0.31,
        :snippets_duration_s => 15,
        :lower_margin_perc => 0.09,
        :upper_margin_perc => 0.93,
        :signal_type => "aEEG"
    )
)