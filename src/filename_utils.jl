function fn2str(func::Function)
    strip(string(func), ['!'])
end
function obj2str(obj)
    strip(string(obj), ['(', ')'])
end

function artifacts_str(excluded_artifact_grades)
    artifacts = if isempty(excluded_artifact_grades) 
        "_artifacts"
    else
        ""
    end
end

function make_signal_stem(signal_type; excluded_artifact_grades, preproc!, postproc!, assumption, conditioned_on, snippets_duration_s, lag_extents, patient_num, unused_params...)
    if !isempty(unused_params)
        @warn "Signal filename provided unused parameters: $unused_params"
    end
    artifacts = artifacts_str(excluded_artifact_grades)
    "$(signal_type)$(artifacts)_$(fn2str(preproc!))_$(fn2str(postproc!))_$(obj2str(assumption))_$(obj2str(conditioned_on))_snippets$(snippets_duration_s)_lagextents$(lag_extents[1])x$(lag_extents[2])_helsinkiEEG$(patient_num)_"
end

function make_tricorr_stem(; params...)
    make_signal_stem("tricorr"; params...)
end

function make_epochdiff_stem(signal_type; min_reviewers_per_seizure, min_dist_to_seizure, params...)
    signal_stem = make_signal_stem("$(signal_type)epochdiff"; params...)
    "$(signal_stem)nrev$(min_reviewers_per_seizure)_szdist$(min_dist_to_seizure)_"
end

function make_detection_stem(signal_type, reduction_type; alert_grace_s, rolling_window_s, params...)
    epochdiff_stem = make_epochdiff_stem("$(signal_type)$(reduction_type)detect"; params...)
    "$(epochdiff_stem)grace$(alert_grace_s)_window$(rolling_window_s)_"
end