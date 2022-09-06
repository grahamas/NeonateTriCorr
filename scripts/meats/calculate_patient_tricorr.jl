function calculate_patient_tricorr(patient_num;
        excluded_artifact_grades=[1],
        eeg = load_helsinki_eeg(patient_num; excluded_artifact_grades=excluded_artifact_grades),
        snippets_duration_s, 
        preproc!, postproc!,
        assumption, conditioned_on,
        lag_extents, 
        plot_traces=true,
        force_recalculate_contributions=false,
        parent_session_id=nothing
    )

    unique_id = if !isnothing(parent_session_id)
        parent_session_id
    else
        Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")
    end
    basename = if isempty(excluded_artifact_grades) 
        "tricorr_artifacts"
    else
        "tricorr"
    end
    target_match_str = make_filename_stem("tricorr"; 
        excluded_artifact_grades=excluded_artifact_grades,
        preproc! = preproc!, postproc! = postproc!,
        assumption=assumption, conditioned_on=conditioned_on,
        lag_extents=lag_extents, patient_num=patient_num
    )
    session_name = "$(target_match_str)$(unique_id)"
    @info "Running: $session_name"

    if preproc! == TripleCorrelations.zscore!
        preproc! = (o,i) -> TripleCorrelations.zscore!(o,i,mean(i),std(i))
    end
    maybe_jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    contributions = if isnothing(maybe_jld_dict) || force_recalculate_contributions
        @info "Calculating triple correlation..."
        calc_class_contributions(eeg, Periodic(), 
                preproc!, postproc!,
                assumption, conditioned_on
                ;
                lag_extents = lag_extents,
                n_motif_classes = 14,
                snippets_duration_s=snippets_duration_s
        )
        save(datadir("exp_pro", "$(session_name).jld2"), 
            Dict(
                "contributions" => contributions,
                "excluded_artifact_grades" => excluded_artifact_grades,
                "snippets_duration_s"=> snippets_duration_s,
                "lag_extents" => lag_extents,
                "assumption" => assumption,
                "conditioned_on" => conditioned_on,
                "preproc_str" => fn2str(preproc!),
                "postproc_str" => fn2str(postproc!)
            )
        )
    else
        @info "Using cached triple correlation!"
        contributions = pop!(maybe_jld_dict, "contributions")
        parameters = Dict(
            "excluded_artifact_grades" => excluded_artifact_grades,
            "snippets_duration_s"=> snippets_duration_s,
            "lag_extents" => lag_extents,
            "assumption" => assumption,
            "conditioned_on" => conditioned_on,
            "preproc_str" => fn2str(preproc!),
            "postproc_str" => fn2str(postproc!)
        )
        @assert maybe_jld_dict == parameters
        return contributions
    end

    if plot_traces
        @info "Plotting EEG and TriCorr motif-class contributions..."
        plots_subdir = plotsdir(session_name)
        mkpath(plots_subdir)

        eeg_fig = draw_eeg_traces(eeg; title = "EEG (Patient $patient_num)", resolution=(1000,1600))
        save(joinpath(plots_subdir, "pat$(patient_num)_eeg_traces.png"), eeg_fig)

        times = get_times(eeg, sample_rate=snippets_duration_s)
        conts_fig = plot_contributions(eeg, times, contributions; title="Motif Contributions (Patient $patient_num)", resolution=(1000,1600))
        save(joinpath(plots_subdir, "pat$(patient_num)_contributions.png"), conts_fig)
    end

    return contributions
end