function calculate_patient_aEEG(patient_num;
        excluded_artifact_grades,
        eeg = load_helsinki_eeg(patient_num; excluded_artifact_grades=excluded_artifact_grades),
        plot_traces=true,
        force_recalculate_aEEG=false,
        parent_session_id=nothing,
        snippets_duration_s,
        params...
    )
    params = Dict(params..., :excluded_artifact_grades => excluded_artifact_grades, :patient_num => patient_num, :snippets_duration_s => snippets_duration_s)
    str_params = convert_keys_to_strings(params)

    unique_id = if !isnothing(parent_session_id)
        parent_session_id
    else
        Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")
    end
    target_match_str = make_aEEG_stem(; params...)
    session_name = "$(target_match_str)$(unique_id)"
    @info "Running: $session_name"

    maybe_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    aEEG = if isnothing(maybe_dict) || force_recalculate_aEEG
        @info "Calculating aEEG..."
        aEEG = calculate_aEEG(eeg; params...)
        save(datadir("exp_pro", "$(session_name).jld2"), Dict(
            str_params...,
            "aEEG" => aEEG,
        )
        )
        aEEG
    else
        @info "Using cached aEEG!"
        aEEG = pop!(maybe_dict, "aEEG")
        @warn "popping unused variable"
        pop!(maybe_dict, "min_dist_to_seizure")
        @assert maybe_dict == str_params "Expected $(maybe_dict); got $(str_params)"
        aEEG
    end

    @warn "Only taking lower margin."
    @show size(aEEG)
    @show typeof(aEEG)
    aEEG = aEEG_lower_margin(aEEG)

    if plot_traces
        @info "Plotting EEG and aEEG..."
        plots_subdir = plotsdir(session_name)
        mkpath(plots_subdir)

        # eeg_fig = draw_eeg_traces(eeg; title = "EEG (Patient $patient_num)", resolution=(1000,1600))
        # save(joinpath(plots_subdir, "pat$(patient_num)_eeg.png"), eeg_fig)

        times = get_times(eeg, sample_rate=1/snippets_duration_s)
        @show length(times)
        @show size(aEEG) typeof(aEEG)
        conts_fig = plot_contributions(eeg, times, aEEG; title="aEEG (Patient $patient_num)", resolution=(1000,1600))
        save(joinpath(plots_subdir, "pat$(patient_num)_aEEG.png"), conts_fig)
    end

    return aEEG
end