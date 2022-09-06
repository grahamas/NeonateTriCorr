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

function make_filename_stem(basename; excluded_artifact_grades, preproc!, postproc!, assumption, conditioned_on, snippets_duration_s, lag_extents, patient_num)
    artifacts = artifacts_str(excluded_artifact_grades)
    "$(basename)$(artifacts)_$(fn2str(preproc!))_$(fn2str(postproc!))_$(obj2str(assumption))_$(obj2str(conditioned_on))_snippets$(snippets_duration_s)_lagextents$(lag_extents[1])x$(lag_extents[2])_helsinkiEEG$(patient_num)_"
end