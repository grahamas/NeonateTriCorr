using CSV
using DataFrames
using Dates
using DSP
using Statistics

struct ArtifactDetectionMethod
    name::String
    description::String
    score_fn::Function
end

function artifact_labeled_patients(; csv_path=scriptsdir("helsinki_artifacts.csv"))
    df = CSV.read(csv_path, DataFrame)
    sort(unique(Int.(df[:, "Patient #"])))
end

function has_present_throughout_artifact(text)
    text = lowercase(text)
    occursin("present throughout", text) ||
        occursin("throughout record", text) ||
        occursin("throughout recording", text)
end

function bounded_artifact_tuple(start, duration, recording_duration; artifact_buffer=1)
    on, off = artifact_tuple(start, duration, artifact_buffer)
    (max(Float64(on), 0.0), min(Float64(off), Float64(recording_duration)))
end

function merge_artifact_bounds!(bounds)
    sort!(bounds)
    if isempty(bounds)
        return bounds
    end

    merged = Tuple{Float64,Float64}[]
    push!(merged, first(bounds))
    for (on, off) in bounds[begin+1:end]
        prev_on, prev_off = merged[end]
        if on <= prev_off
            merged[end] = (prev_on, max(prev_off, off))
        else
            push!(merged, (on, off))
        end
    end
    empty!(bounds)
    append!(bounds, merged)
    bounds
end

function load_helsinki_artifact_ground_truth(
        eeg_num,
        artifact_grades=[1, 2];
        csv_path=scriptsdir("helsinki_artifacts.csv"),
        start_time::Time=Time(EDF.read(datadir("exp_raw", "helsinki", "eeg$(eeg_num).edf")).header.start),
        recording_duration,
        artifact_buffer=1
    )
    df = CSV.read(csv_path, DataFrame)
    subset!(df, "Patient #" => ByRow(==(eeg_num)))

    artifact_bounds = Tuple{Float64,Float64}[]
    for row in eachrow(df)
        grade = parse_grade(row.Text)
        if grade ∉ artifact_grades
            continue
        end

        if has_present_throughout_artifact(row.Text)
            push!(artifact_bounds, (0.0, Float64(recording_duration)))
        else
            start = parse_start(row.Time, start_time)
            duration = parse_artifact_duration(row.Duration)
            on, off = bounded_artifact_tuple(start, duration, recording_duration; artifact_buffer=artifact_buffer)
            if off > on
                push!(artifact_bounds, (on, off))
            end
        end
    end

    merge_artifact_bounds!(artifact_bounds)
    artifact_bounds
end

function artifact_bin_bounds(recording_duration, bin_s)
    bin_starts = collect(0.0:bin_s:Float64(recording_duration))
    if isempty(bin_starts) || bin_starts[end] < recording_duration
        push!(bin_starts, Float64(recording_duration))
    end
    if length(bin_starts) == 1
        return Float64[], Float64[]
    end
    (bin_starts[begin:end-1], bin_starts[begin+1:end])
end

function artifact_sample_ranges(eeg::AbstractProcessedEEG, bin_starts, bin_stops)
    map(zip(bin_starts, bin_stops)) do (bin_start, bin_stop)
        start_idx = floor(Int, bin_start * eeg.sample_rate) + 1
        stop_idx = min(floor(Int, bin_stop * eeg.sample_rate), size(eeg.signals, 2))
        start_idx:stop_idx
    end
end

function bounds_to_bin_labels(bounds, bin_starts, bin_stops)
    map(zip(bin_starts, bin_stops)) do (bin_start, bin_stop)
        any(bounds) do (on, off)
            on < bin_stop && off > bin_start
        end
    end
end

function seizure_event_ids_for_bins(seizure_bounds, bin_starts, bin_stops)
    event_ids = zeros(Int, length(bin_starts))
    for (event_id, (on, off)) in enumerate(seizure_bounds)
        for idx in eachindex(bin_starts)
            if on < bin_stops[idx] && off > bin_starts[idx]
                event_ids[idx] = event_id
            end
        end
    end
    event_ids
end

function finite_or_nan(x)
    isfinite(x) ? x : NaN
end

function robust_center_and_scale(xs)
    center = median(xs)
    scale = 1.4826 * median(abs.(xs .- center))
    if !isfinite(scale) || scale <= eps(Float64)
        scale = std(xs)
    end
    if !isfinite(scale) || scale <= eps(Float64)
        scale = 1.0
    end
    (center, scale)
end

function robust_channel_stats(signals::AbstractMatrix)
    centers = zeros(Float64, size(signals, 1))
    scales = zeros(Float64, size(signals, 1))
    for channel in axes(signals, 1)
        centers[channel], scales[channel] = robust_center_and_scale(vec(signals[channel, :]))
    end
    centers, scales
end

function score_amplitude_robust_z(eeg::AbstractProcessedEEG, sample_ranges)
    signals = Matrix(eeg.signals)
    centers, scales = robust_channel_stats(signals)
    scores = zeros(Float64, length(sample_ranges))
    for (bin_idx, sample_range) in enumerate(sample_ranges)
        score = 0.0
        for channel in axes(signals, 1)
            center = centers[channel]
            scale = scales[channel]
            for sample_idx in sample_range
                score = max(score, abs((signals[channel, sample_idx] - center) / scale))
            end
        end
        scores[bin_idx] = score
    end
    scores
end

function score_jump_robust_z(eeg::AbstractProcessedEEG, sample_ranges)
    signals = Matrix(eeg.signals)
    differences = diff(signals, dims=2)
    centers, scales = robust_channel_stats(differences)
    scores = zeros(Float64, length(sample_ranges))
    for (bin_idx, sample_range) in enumerate(sample_ranges)
        if length(sample_range) < 2
            scores[bin_idx] = 0.0
            continue
        end

        diff_range = first(sample_range):min(last(sample_range) - 1, size(differences, 2))
        score = 0.0
        for channel in axes(differences, 1)
            center = centers[channel]
            scale = scales[channel]
            for sample_idx in diff_range
                score = max(score, abs((differences[channel, sample_idx] - center) / scale))
            end
        end
        scores[bin_idx] = score
    end
    scores
end

function periodogram_bandpower(segment, fs, bands)
    p = periodogram(segment; fs=fs)
    frequencies = freq(p)
    powers = power(p)
    total = 0.0
    for (low, high) in bands
        for idx in eachindex(frequencies)
            if low <= frequencies[idx] < high
                total += powers[idx]
            end
        end
    end
    total
end

function score_frequency_ratio(eeg::AbstractProcessedEEG, sample_ranges; numerator_bands, denominator_bands)
    signals = Matrix(eeg.signals)
    scores = zeros(Float64, length(sample_ranges))
    for (bin_idx, sample_range) in enumerate(sample_ranges)
        score = -Inf
        for channel in axes(signals, 1)
            segment = vec(signals[channel, sample_range])
            numerator = periodogram_bandpower(segment, eeg.sample_rate, numerator_bands)
            denominator = periodogram_bandpower(segment, eeg.sample_rate, denominator_bands)
            score = max(score, log((numerator + eps(Float64)) / (denominator + eps(Float64))))
        end
        scores[bin_idx] = finite_or_nan(score)
    end
    scores
end

function score_high_frequency_ratio(eeg::AbstractProcessedEEG, sample_ranges)
    score_frequency_ratio(eeg, sample_ranges;
        numerator_bands=[(30.0, 45.0), (55.0, 70.0)],
        denominator_bands=[(1.0, 30.0)]
    )
end

function score_low_frequency_ratio(eeg::AbstractProcessedEEG, sample_ranges)
    score_frequency_ratio(eeg, sample_ranges;
        numerator_bands=[(0.5, 4.0)],
        denominator_bands=[(4.0, 30.0)]
    )
end

function score_flatline_dropout(eeg::AbstractProcessedEEG, sample_ranges)
    signals = Matrix(eeg.signals)
    bin_stds = zeros(Float64, size(signals, 1), length(sample_ranges))
    for (bin_idx, sample_range) in enumerate(sample_ranges)
        for channel in axes(signals, 1)
            bin_stds[channel, bin_idx] = length(sample_range) > 1 ? std(signals[channel, sample_range]) : 0.0
        end
    end

    centers = zeros(Float64, size(signals, 1))
    scales = zeros(Float64, size(signals, 1))
    for channel in axes(signals, 1)
        centers[channel], scales[channel] = robust_center_and_scale(vec(bin_stds[channel, :]))
    end

    scores = zeros(Float64, length(sample_ranges))
    for bin_idx in eachindex(scores)
        score = 0.0
        for channel in axes(signals, 1)
            score = max(score, max(0.0, (centers[channel] - bin_stds[channel, bin_idx]) / scales[channel]))
        end
        scores[bin_idx] = score
    end
    scores
end

function artifact_detection_methods()
    [
        ArtifactDetectionMethod(
            "amplitude_robust_z",
            "Maximum absolute robust z-score across channels and samples.",
            score_amplitude_robust_z
        ),
        ArtifactDetectionMethod(
            "jump_robust_z",
            "Maximum absolute robust z-score of first differences.",
            score_jump_robust_z
        ),
        ArtifactDetectionMethod(
            "high_frequency_ratio",
            "Log high-frequency power ratio: (30-45 Hz plus 55-70 Hz) / 1-30 Hz.",
            score_high_frequency_ratio
        ),
        ArtifactDetectionMethod(
            "low_frequency_ratio",
            "Log low-frequency power ratio: 0.5-4 Hz / 4-30 Hz.",
            score_low_frequency_ratio
        ),
        ArtifactDetectionMethod(
            "flatline_dropout",
            "Robust low-variance outlier score from per-bin channel standard deviation.",
            score_flatline_dropout
        )
    ]
end

function artifact_score_rows_for_patient(
        patient_num;
        methods=artifact_detection_methods(),
        bin_s=15,
        artifact_grades=[1, 2],
        min_reviewers_per_seizure=3,
        artifact_csv_path=scriptsdir("helsinki_artifacts.csv")
    )
    eeg = load_helsinki_eeg(patient_num;
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        excluded_artifact_grades=Int[]
    )
    artifact_bounds = load_helsinki_artifact_ground_truth(patient_num, artifact_grades;
        csv_path=artifact_csv_path,
        recording_duration=eeg.duration
    )
    seizure_bounds = eeg.seizure_annotations
    bin_starts, bin_stops = artifact_bin_bounds(eeg.duration, bin_s)
    sample_ranges = artifact_sample_ranges(eeg, bin_starts, bin_stops)
    artifact_truth = bounds_to_bin_labels(artifact_bounds, bin_starts, bin_stops)
    seizure_truth = bounds_to_bin_labels(seizure_bounds, bin_starts, bin_stops)
    seizure_event_ids = seizure_event_ids_for_bins(seizure_bounds, bin_starts, bin_stops)

    rows = NamedTuple[]
    for method in methods
        scores = method.score_fn(eeg, sample_ranges)
        append!(rows, map(eachindex(scores)) do bin_idx
            (
                patient=patient_num,
                method=method.name,
                method_description=method.description,
                bin_index=bin_idx,
                bin_start=bin_starts[bin_idx],
                bin_stop=bin_stops[bin_idx],
                artifact_truth=artifact_truth[bin_idx],
                seizure_truth=seizure_truth[bin_idx],
                clean_seizure_truth=seizure_truth[bin_idx] && !artifact_truth[bin_idx],
                seizure_event_id=seizure_event_ids[bin_idx],
                score=scores[bin_idx]
            )
        end)
    end

    DataFrame(rows)
end

ratio_or_nan(num, denom) = denom == 0 ? NaN : num / denom

function artifact_confusion_metrics(predicted_artifact, artifact_truth)
    true_positives = count(predicted_artifact .& artifact_truth)
    false_positives = count(predicted_artifact .& .!artifact_truth)
    true_negatives = count(.!predicted_artifact .& .!artifact_truth)
    false_negatives = count(.!predicted_artifact .& artifact_truth)
    precision = ratio_or_nan(true_positives, true_positives + false_positives)
    recall = ratio_or_nan(true_positives, true_positives + false_negatives)
    specificity = ratio_or_nan(true_negatives, true_negatives + false_positives)
    f1 = (isnan(precision) || isnan(recall) || precision + recall == 0) ? NaN : 2 * precision * recall / (precision + recall)
    (
        true_positives=true_positives,
        false_positives=false_positives,
        true_negatives=true_negatives,
        false_negatives=false_negatives,
        precision=precision,
        recall=recall,
        specificity=specificity,
        f1=f1
    )
end

function seizure_as_artifact_metrics(predicted_artifact, artifact_truth, seizure_truth, seizure_event_ids)
    seizure_bins = count(seizure_truth)
    clean_seizure_truth = seizure_truth .& .!artifact_truth
    clean_seizure_bins = count(clean_seizure_truth)
    predicted_artifact_bins = count(predicted_artifact)
    seizure_bins_flagged_as_artifact = count(predicted_artifact .& seizure_truth)
    clean_seizure_bins_flagged_as_artifact = count(predicted_artifact .& clean_seizure_truth)
    seizure_events = Set(filter(!=(0), seizure_event_ids))
    clean_seizure_events = Set(seizure_event_ids[clean_seizure_truth])
    delete!(clean_seizure_events, 0)
    clean_flagged_event_ids = Set(seizure_event_ids[predicted_artifact .& clean_seizure_truth])
    delete!(clean_flagged_event_ids, 0)

    (
        seizure_bins=seizure_bins,
        seizure_bins_flagged_as_artifact=seizure_bins_flagged_as_artifact,
        seizure_bin_flag_fraction=ratio_or_nan(seizure_bins_flagged_as_artifact, seizure_bins),
        clean_seizure_bins=clean_seizure_bins,
        clean_seizure_bins_flagged_as_artifact=clean_seizure_bins_flagged_as_artifact,
        clean_seizure_bin_flag_fraction=ratio_or_nan(clean_seizure_bins_flagged_as_artifact, clean_seizure_bins),
        seizure_events=length(seizure_events),
        clean_seizure_events=length(clean_seizure_events),
        clean_seizure_events_flagged_as_artifact=length(clean_flagged_event_ids),
        clean_seizure_event_flag_fraction=ratio_or_nan(length(clean_flagged_event_ids), length(clean_seizure_events)),
        predicted_artifact_bins=predicted_artifact_bins,
        predicted_artifact_bins_overlapping_seizure=seizure_bins_flagged_as_artifact,
        predicted_artifact_seizure_overlap_fraction=ratio_or_nan(seizure_bins_flagged_as_artifact, predicted_artifact_bins)
    )
end

function auc_from_scores(scores, truth)
    finite_idxs = findall(isfinite.(scores))
    scores = scores[finite_idxs]
    truth = truth[finite_idxs]
    positives = count(truth)
    negatives = length(truth) - positives
    if positives == 0 || negatives == 0
        return NaN
    end

    order = sortperm(scores, rev=true)
    tp = 0
    fp = 0
    prev_tpr = 0.0
    prev_fpr = 0.0
    auc = 0.0
    idx = firstindex(order)
    while idx <= lastindex(order)
        score = scores[order[idx]]
        group_tp = 0
        group_fp = 0
        while idx <= lastindex(order) && scores[order[idx]] == score
            if truth[order[idx]]
                group_tp += 1
            else
                group_fp += 1
            end
            idx += 1
        end
        tp += group_tp
        fp += group_fp
        tpr = tp / positives
        fpr = fp / negatives
        auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2
        prev_tpr = tpr
        prev_fpr = fpr
    end
    auc
end

function auprc_from_scores(scores, truth)
    finite_idxs = findall(isfinite.(scores))
    scores = scores[finite_idxs]
    truth = truth[finite_idxs]
    positives = count(truth)
    if positives == 0
        return NaN
    end

    order = sortperm(scores, rev=true)
    tp = 0
    fp = 0
    prev_recall = 0.0
    auprc = 0.0
    idx = firstindex(order)
    while idx <= lastindex(order)
        score = scores[order[idx]]
        group_tp = 0
        group_fp = 0
        while idx <= lastindex(order) && scores[order[idx]] == score
            if truth[order[idx]]
                group_tp += 1
            else
                group_fp += 1
            end
            idx += 1
        end
        tp += group_tp
        fp += group_fp
        recall = tp / positives
        precision = tp / (tp + fp)
        auprc += precision * (recall - prev_recall)
        prev_recall = recall
    end
    auprc
end

function threshold_values(scores; n_thresholds=100)
    finite_scores = collect(skipmissing(scores[isfinite.(scores)]))
    if isempty(finite_scores) || n_thresholds <= 0
        return Float64[]
    end

    min_score, max_score = extrema(finite_scores)
    pad = max(eps(Float64), eps(Float64) * max(abs(min_score), abs(max_score)))
    if n_thresholds == 1
        return [median(finite_scores)]
    elseif n_thresholds == 2
        return [max_score + pad, min_score - pad]
    end

    interior_threshold_count = n_thresholds - 2
    quantile_points = interior_threshold_count == 1 ? [0.5] : range(0.0, 1.0, length=interior_threshold_count)
    score_values = if length(unique(finite_scores)) <= interior_threshold_count
        unique(finite_scores)
    else
        quantile(finite_scores, quantile_points)
    end
    sort(unique(vcat(max_score + pad, score_values, min_score - pad)), rev=true)
end

function top_fraction_predictions(scores; fraction=0.05)
    finite_idxs = findall(isfinite.(scores))
    predictions = falses(length(scores))
    if isempty(finite_idxs) || fraction <= 0
        return predictions
    end

    n_predicted = max(1, ceil(Int, fraction * length(finite_idxs)))
    ordered_local = sortperm(scores[finite_idxs], rev=true)
    predicted_idxs = finite_idxs[ordered_local[1:min(n_predicted, length(ordered_local))]]
    predictions[predicted_idxs] .= true
    predictions
end

function offset_event_ids(patient_df::DataFrame, event_offset)
    event_ids = copy(patient_df.seizure_event_id)
    positive_ids = event_ids .> 0
    event_ids[positive_ids] .+= event_offset
    event_ids
end

function method_metric_row(method_df::AbstractDataFrame, predicted_artifact; method, operating_point, threshold=NaN)
    artifact_truth = Vector{Bool}(method_df.artifact_truth)
    seizure_truth = Vector{Bool}(method_df.seizure_truth)
    event_ids = Vector{Int}(method_df.seizure_event_id)
    merge(
        (
            method=method,
            operating_point=operating_point,
            threshold=threshold,
            n_bins=nrow(method_df),
            artifact_bins=count(artifact_truth),
            nonartifact_bins=count(.!artifact_truth)
        ),
        artifact_confusion_metrics(predicted_artifact, artifact_truth),
        seizure_as_artifact_metrics(predicted_artifact, artifact_truth, seizure_truth, event_ids)
    )
end

function threshold_sweep_for_method(method_df::AbstractDataFrame; n_thresholds=100)
    rows = NamedTuple[]
    scores = Vector{Float64}(method_df.score)
    method = first(method_df.method)
    for threshold in threshold_values(scores; n_thresholds=n_thresholds)
        predicted_artifact = isfinite.(scores) .& (scores .>= threshold)
        push!(rows, method_metric_row(method_df, predicted_artifact;
            method=method,
            operating_point="threshold_sweep",
            threshold=threshold
        ))
    end
    DataFrame(rows)
end

function fixed_top_fraction_patient_rows(score_df::DataFrame; fraction=0.05)
    rows = NamedTuple[]
    for patient_df in groupby(score_df, [:patient, :method])
        scores = Vector{Float64}(patient_df.score)
        predicted_artifact = top_fraction_predictions(scores; fraction=fraction)
        threshold = any(predicted_artifact) ? minimum(scores[predicted_artifact]) : NaN
        push!(rows, merge(
            (patient=first(patient_df.patient),),
            method_metric_row(patient_df, predicted_artifact;
                method=first(patient_df.method),
                operating_point="top_$(round(Int, 100 * fraction))pct_per_patient",
                threshold=threshold
            )
        ))
    end
    DataFrame(rows)
end

function fixed_top_fraction_summary_rows(score_df::DataFrame; fraction=0.05)
    rows = NamedTuple[]
    for method_df in groupby(score_df, :method)
        predicted = falses(nrow(method_df))
        for patient in unique(method_df.patient)
            patient_mask = method_df.patient .== patient
            local_predictions = top_fraction_predictions(Vector{Float64}(method_df.score[patient_mask]); fraction=fraction)
            predicted[patient_mask] .= local_predictions
        end
        scores = Vector{Float64}(method_df.score)
        threshold = any(predicted) ? minimum(scores[predicted]) : NaN
        push!(rows, method_metric_row(method_df, predicted;
            method=first(method_df.method),
            operating_point="top_$(round(Int, 100 * fraction))pct_per_patient",
            threshold=threshold
        ))
    end
    DataFrame(rows)
end

function add_score_auc_columns!(summary_df::DataFrame, score_df::DataFrame)
    aurocs = Dict{String,Float64}()
    auprcs = Dict{String,Float64}()
    for method_df in groupby(score_df, :method)
        method = first(method_df.method)
        scores = Vector{Float64}(method_df.score)
        truth = Vector{Bool}(method_df.artifact_truth)
        aurocs[method] = auc_from_scores(scores, truth)
        auprcs[method] = auprc_from_scores(scores, truth)
    end
    summary_df.auroc = [aurocs[row.method] for row in eachrow(summary_df)]
    summary_df.auprc = [auprcs[row.method] for row in eachrow(summary_df)]
    summary_df
end

function plot_artifact_roc_curves(threshold_df::DataFrame; resolution=(1200, 800))
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1]; xlabel="False positive rate", ylabel="Artifact recall", title="Artifact ROC sweep")
    for method_df in groupby(threshold_df, :method)
        sorted_df = DataFrame(method_df)
        sort!(sorted_df, [:false_positives, :true_positives])
        fpr = sorted_df.false_positives ./ (sorted_df.false_positives .+ sorted_df.true_negatives)
        tpr = sorted_df.recall
        valid = .!isnan.(fpr) .& .!isnan.(tpr)
        lines!(ax, fpr[valid], tpr[valid], label=first(sorted_df.method), linewidth=3)
    end
    axislegend(ax; position=:rb)
    fig
end

function plot_artifact_summary_bars(summary_df::DataFrame; resolution=(1400, 800))
    methods = summary_df.method
    x = collect(1:length(methods))
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1];
        ylabel="Fraction",
        title="Fixed top-5% artifact predictions",
        xticks=(x, methods),
        xticklabelrotation=pi / 5
    )
    barplot!(ax, x .- 0.25, summary_df.precision; width=0.22, label="precision")
    barplot!(ax, x, summary_df.recall; width=0.22, label="recall")
    barplot!(ax, x .+ 0.25, summary_df.clean_seizure_bin_flag_fraction; width=0.22, label="clean seizure flagged")
    ylims!(ax, 0, 1)
    axislegend(ax; position=:rt)
    fig
end

function plot_artifact_seizure_tradeoff(threshold_df::DataFrame; resolution=(1200, 800))
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1];
        xlabel="Artifact recall",
        ylabel="Clean seizure bins flagged as artifact",
        title="Artifact detection vs seizure contamination"
    )
    for method_df in groupby(threshold_df, :method)
        valid = .!isnan.(method_df.recall) .& .!isnan.(method_df.clean_seizure_bin_flag_fraction)
        lines!(ax, method_df.recall[valid], method_df.clean_seizure_bin_flag_fraction[valid],
            label=first(method_df.method), linewidth=3)
    end
    axislegend(ax; position=:rb)
    fig
end

function save_artifact_detection_plots(summary_df, threshold_df, plot_dir)
    mkpath(plot_dir)
    save(joinpath(plot_dir, "artifact_roc_sweep.png"), plot_artifact_roc_curves(threshold_df))
    save(joinpath(plot_dir, "artifact_fixed_top5_summary.png"), plot_artifact_summary_bars(summary_df))
    save(joinpath(plot_dir, "artifact_seizure_tradeoff.png"), plot_artifact_seizure_tradeoff(threshold_df))
end

function run_artifact_detection_survey(;
        patients=artifact_labeled_patients(),
        methods=artifact_detection_methods(),
        bin_s=15,
        artifact_grades=[1, 2],
        min_reviewers_per_seizure=3,
        n_thresholds=100,
        fixed_fraction=0.05,
        output_root=datadir("exp_pro", "artifact_detection"),
        plot_root=plotsdir("artifact_detection_$(Dates.now())"),
        artifact_csv_path=scriptsdir("helsinki_artifacts.csv"),
        save_outputs=true
    )
    score_dfs = DataFrame[]
    event_offset = 0
    for patient in patients
        @info "Scoring artifact methods" patient
        patient_df = artifact_score_rows_for_patient(patient;
            methods=methods,
            bin_s=bin_s,
            artifact_grades=artifact_grades,
            min_reviewers_per_seizure=min_reviewers_per_seizure,
            artifact_csv_path=artifact_csv_path
        )
        local_event_count = maximum(patient_df.seizure_event_id)
        patient_df.seizure_event_id = offset_event_ids(patient_df, event_offset)
        event_offset += local_event_count
        push!(score_dfs, patient_df)
    end

    score_df = reduce(vcat, score_dfs)
    threshold_df = reduce(vcat, [
        threshold_sweep_for_method(DataFrame(method_df); n_thresholds=n_thresholds)
        for method_df in groupby(score_df, :method)
    ])
    patient_df = fixed_top_fraction_patient_rows(score_df; fraction=fixed_fraction)
    summary_df = fixed_top_fraction_summary_rows(score_df; fraction=fixed_fraction)
    add_score_auc_columns!(summary_df, score_df)

    if save_outputs
        session_id = Dates.now()
        output_dir = joinpath(output_root, string(session_id))
        mkpath(output_dir)
        CSV.write(joinpath(output_dir, "artifact_detection_summary.csv"), summary_df)
        CSV.write(joinpath(output_dir, "artifact_detection_per_patient.csv"), patient_df)
        CSV.write(joinpath(output_dir, "artifact_detection_threshold_sweep.csv"), threshold_df)
        CSV.write(joinpath(output_dir, "artifact_detection_bin_scores.csv"), score_df)
        save_artifact_detection_plots(summary_df, threshold_df, plot_root)
        @info "Saved artifact detection survey" output_dir plot_root
        return (
            summary_df=summary_df,
            patient_df=patient_df,
            threshold_df=threshold_df,
            score_df=score_df,
            output_dir=output_dir,
            plot_dir=plot_root
        )
    else
        return (
            summary_df=summary_df,
            patient_df=patient_df,
            threshold_df=threshold_df,
            score_df=score_df,
            output_dir=nothing,
            plot_dir=nothing
        )
    end
end
