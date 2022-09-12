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

function make_aEEG_stem(; excluded_artifact_grades, patient_num="", lowpass_freq, snippets_duration_s, lower_margin_perc, upper_margin_perc, unused_params...)
    if !isempty(unused_params)
        @warn "Signal filename provided unused parameters: $unused_params"
    end
    artifacts = artifacts_str(excluded_artifact_grades)
    "aEEG$(artifacts)_snippets$(snippets_duration_s)_lowpass$(lowpass_freq)_lmargin$(lower_margin_perc)_umargin$(upper_margin_perc)_helsinkiEEG$(patient_num)_"
end

function make_tricorr_stem(; excluded_artifact_grades, preproc!, postproc!, assumption, conditioned_on, snippets_duration_s, lag_extents, patient_num="", unused_params...)
    if !isempty(unused_params)
        @warn "Signal filename provided unused parameters: $unused_params"
    end
    artifacts = artifacts_str(excluded_artifact_grades)
    "tricorr$(artifacts)_$(fn2str(preproc!))_$(fn2str(postproc!))_$(obj2str(assumption))_$(obj2str(conditioned_on))_snippets$(snippets_duration_s)_lagextents$(lag_extents[1])x$(lag_extents[2])_helsinkiEEG$(patient_num)_"
end

function make_signal_stem(signal_type; params...)
    if signal_type == "tricorr"
        make_tricorr_stem(; params...)
    elseif signal_type == "aEEG"
        make_aEEG_stem(; params...)
    else
        @error "Unsupported signal type: $signal_type"
    end
end

function make_epochdiff_stem(signal_type; min_reviewers_per_seizure, min_dist_to_seizure, params...)
    signal_stem = make_signal_stem(signal_type; params...)
    "epochdiff$(signal_stem)nrev$(min_reviewers_per_seizure)_szdist$(min_dist_to_seizure)_"
end

function make_detection_stem(signal_type, reduction_type; alert_grace_s, signals_reduction_params, params...)
    epochdiff_stem = make_epochdiff_stem(signal_type; params...)
    reductionstr = make_reduction_str(reduction_type; signals_reduction_params...)
    "detect$(reduction_type)$(epochdiff_stem)grace$(alert_grace_s)_$(reductionstr)_"
end

function make_reduction_str(reduction_type; params...)
    if reduction_type == "maxany"
        make_reduction_maxany_str(; params...)
    elseif reduction_type == "meansignificant"
        make_reduction_meansignificant_str(; params...)
    elseif reduction_type == "meanall"
        make_reduction_meanall_str(; params...)
    else
        @error "Unsupported reduction type: $(reduction_type)"
    end
end

function make_reduction_maxany_str(; rolling_window_s, unused...)
    "window$(rolling_window_s)"
end
make_reduction_meanall_str = make_reduction_maxany_str
function make_reduction_meansignificant_str(; rolling_window_s, n_signals_used, unused...)
    "window$(rolling_window_s)_nsigs$(n_signals_used)"
end

function convert_keys_to_strings(dct)
    #shallow
    Dict(
        string(key) => val for (key, val) ∈ pairs(dct)
    )
end