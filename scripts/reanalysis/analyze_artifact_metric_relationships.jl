#!/usr/bin/env julia

# Analyze how already-measured metric columns relate to artifact_truth.
#
# This script deliberately uses only Julia stdlibs. The project CSV/DataFrames
# stack is the normal choice, but keeping this dependency-free makes it usable
# even when the local package precompile cache is unhappy.

using Dates
using LinearAlgebra
using Printf
using Random
using Statistics

const ID_COLUMNS = Set([
    "patient",
    "channel",
    "channel_label",
    "bin_index",
    "bin_start",
    "bin_stop",
    "sample_start_index",
    "sample_stop_index",
    "sample_count",
    "artifact_truth",
    "seizure_truth",
    "clean_seizure_truth",
    "seizure_event_id",
])

const SUMMARY_COLUMNS = [
    :feature,
    :family,
    :n,
    :n_finite,
    :artifact_finite,
    :other_finite,
    :artifact_prevalence,
    :auc_high_means_artifact,
    :directional_auc,
    :rank_biserial_signed,
    :direction,
    :artifact_median,
    :other_median,
    :median_delta_artifact_minus_other,
    :robust_shift_iqr,
    :artifact_q25,
    :artifact_q75,
    :other_q25,
    :other_q75,
    :patient_auc_count,
    :patient_directional_auc_median,
    :patient_direction_match_fraction,
]

const REGRESSION_COLUMNS = [
    :feature,
    :family,
    :coefficient_standardized,
    :abs_coefficient_standardized,
    :direction,
    :center,
    :scale,
    :univariate_directional_auc,
    :univariate_direction,
    :univariate_patient_direction_match_fraction,
]

function usage()
    println("""
    Usage:
      julia scripts/reanalysis/analyze_artifact_metric_relationships.jl \\
        --whole-hvg-csv PATH \\
        [--whole-hvg-large-csv PATH] \\
        [--output-dir PATH] \\
        [--large-sample-rows N] \\
        [--seed N]
    """)
end

function parse_args(args)
    parsed = Dict{String,String}()
    idx = 1
    while idx <= length(args)
        arg = args[idx]
        if !startswith(arg, "--")
            error("Unexpected positional argument: $(arg)")
        end
        key = arg[3:end]
        if key in ("help", "h")
            parsed[key] = "true"
            idx += 1
        else
            idx == length(args) && error("Missing value for $(arg)")
            parsed[key] = args[idx + 1]
            idx += 2
        end
    end
    parsed
end

split_csv_line(line::AbstractString) = split(chomp(line), ',')

function read_header(path::AbstractString)
    open(path, "r") do io
        return String.(split_csv_line(readline(io)))
    end
end

function metric_feature_indices(header; prefix="hvg_")
    indices = Int[]
    names = String[]
    for (idx, name) in enumerate(header)
        if startswith(name, prefix) && !(name in ID_COLUMNS)
            push!(indices, idx)
            push!(names, name)
        end
    end
    indices, names
end

function parse_bool(value)
    text = lowercase(String(value))
    if text == "true" || text == "1"
        return true
    elseif text == "false" || text == "0"
        return false
    end
    error("Cannot parse Bool from $(value)")
end

function parse_float_or_nan(value)
    isempty(value) && return NaN
    parsed = tryparse(Float64, value)
    isnothing(parsed) ? NaN : parsed
end

function metric_family(feature::AbstractString)
    for suffix in ("_mean", "_std", "_max", "_q95")
        if endswith(feature, suffix)
            return feature[1:end-length(suffix)]
        end
    end
    feature
end

function load_metric_columns(path::AbstractString, feature_indices, feature_names)
    header = read_header(path)
    patient_idx = findfirst(==("patient"), header)
    truth_idx = findfirst(==("artifact_truth"), header)
    isnothing(patient_idx) && error("No patient column in $(path)")
    isnothing(truth_idx) && error("No artifact_truth column in $(path)")

    patients = Int[]
    truth = Bool[]
    columns = [Float64[] for _ in feature_names]
    total_rows = 0

    open(path, "r") do io
        readline(io)
        for line in eachline(io)
            total_rows += 1
            fields = split_csv_line(line)
            push!(patients, parse(Int, fields[patient_idx]))
            push!(truth, parse_bool(fields[truth_idx]))
            for (out_idx, field_idx) in enumerate(feature_indices)
                push!(columns[out_idx], parse_float_or_nan(fields[field_idx]))
            end
        end
    end

    (; feature_names, patients, truth, columns, total_rows, sampled=false)
end

function sample_metric_columns(
        path::AbstractString,
        feature_indices,
        feature_names;
        sample_rows::Integer,
        seed::Integer
    )
    header = read_header(path)
    patient_idx = findfirst(==("patient"), header)
    truth_idx = findfirst(==("artifact_truth"), header)
    isnothing(patient_idx) && error("No patient column in $(path)")
    isnothing(truth_idx) && error("No artifact_truth column in $(path)")

    rng = MersenneTwister(seed)
    patients = Vector{Int}(undef, sample_rows)
    truth = Vector{Bool}(undef, sample_rows)
    columns = [Vector{Float64}(undef, sample_rows) for _ in feature_names]
    total_rows = 0
    kept_rows = 0

    open(path, "r") do io
        readline(io)
        for line in eachline(io)
            total_rows += 1
            slot = 0
            if total_rows <= sample_rows
                kept_rows = total_rows
                slot = total_rows
            else
                replacement = rand(rng, 1:total_rows)
                if replacement <= sample_rows
                    slot = replacement
                end
            end

            slot == 0 && continue

            fields = split_csv_line(line)
            patients[slot] = parse(Int, fields[patient_idx])
            truth[slot] = parse_bool(fields[truth_idx])
            for (out_idx, field_idx) in enumerate(feature_indices)
                columns[out_idx][slot] = parse_float_or_nan(fields[field_idx])
            end
        end
    end

    if kept_rows < sample_rows
        resize!(patients, kept_rows)
        resize!(truth, kept_rows)
        foreach(column -> resize!(column, kept_rows), columns)
    end

    (; feature_names, patients, truth, columns, total_rows, sampled=true)
end

function finite_subset(values, truth, target)
    out = Float64[]
    for idx in eachindex(values)
        value = values[idx]
        if truth[idx] == target && isfinite(value)
            push!(out, value)
        end
    end
    out
end

function quantile_from_sorted(sorted_values, q)
    n = length(sorted_values)
    n == 0 && return NaN
    n == 1 && return sorted_values[1]
    pos = 1 + (n - 1) * q
    lo = floor(Int, pos)
    hi = ceil(Int, pos)
    lo == hi && return sorted_values[lo]
    sorted_values[lo] + (pos - lo) * (sorted_values[hi] - sorted_values[lo])
end

function safe_quantile(values, q)
    isempty(values) && return NaN
    sorted_values = sort(values)
    quantile_from_sorted(sorted_values, q)
end

function median_or_nan(values)
    safe_quantile(values, 0.5)
end

function auc_from_scores(values, truth)
    idxs = Int[]
    for idx in eachindex(values)
        isfinite(values[idx]) && push!(idxs, idx)
    end
    positives = count(idx -> truth[idx], idxs)
    negatives = length(idxs) - positives
    if positives == 0 || negatives == 0
        return NaN
    end

    sort!(idxs; by=idx -> values[idx])
    rank_sum_positive = 0.0
    i = 1
    while i <= length(idxs)
        j = i
        value = values[idxs[i]]
        while j < length(idxs) && values[idxs[j + 1]] == value
            j += 1
        end
        average_rank = (i + j) / 2
        for k in i:j
            if truth[idxs[k]]
                rank_sum_positive += average_rank
            end
        end
        i = j + 1
    end

    u_positive = rank_sum_positive - positives * (positives + 1) / 2
    u_positive / (positives * negatives)
end

function patient_directional_stats(values, truth, patients, direction)
    patient_aucs = Float64[]
    for patient in sort(unique(patients))
        patient_values = Float64[]
        patient_truth = Bool[]
        for idx in eachindex(values)
            if patients[idx] == patient && isfinite(values[idx])
                push!(patient_values, values[idx])
                push!(patient_truth, truth[idx])
            end
        end
        any(patient_truth) || continue
        any(.!patient_truth) || continue
        auc = auc_from_scores(patient_values, patient_truth)
        isfinite(auc) && push!(patient_aucs, auc)
    end

    if isempty(patient_aucs) || direction == "unavailable"
        return 0, NaN, NaN
    end

    directional = direction == "higher_in_artifact" ? patient_aucs : 1 .- patient_aucs
    (
        length(patient_aucs),
        median_or_nan(directional),
        count(>(0.5), directional) / length(directional)
    )
end

function summarize_feature(feature, values, truth, patients)
    finite_count = count(isfinite, values)
    artifact_values = finite_subset(values, truth, true)
    other_values = finite_subset(values, truth, false)
    artifact_prevalence = count(truth) / length(truth)
    auc = auc_from_scores(values, truth)

    direction = if !isfinite(auc)
        "unavailable"
    elseif auc >= 0.5
        "higher_in_artifact"
    else
        "lower_in_artifact"
    end
    directional_auc = !isfinite(auc) ? NaN : max(auc, 1 - auc)

    artifact_q25 = safe_quantile(artifact_values, 0.25)
    artifact_median = safe_quantile(artifact_values, 0.5)
    artifact_q75 = safe_quantile(artifact_values, 0.75)
    other_q25 = safe_quantile(other_values, 0.25)
    other_median = safe_quantile(other_values, 0.5)
    other_q75 = safe_quantile(other_values, 0.75)
    pooled_iqr = median_or_nan([
        artifact_q75 - artifact_q25,
        other_q75 - other_q25,
    ])
    median_delta = artifact_median - other_median
    robust_shift = isfinite(pooled_iqr) && pooled_iqr > 0 ? median_delta / pooled_iqr : NaN
    patient_auc_count, patient_auc_median, patient_match_fraction =
        patient_directional_stats(values, truth, patients, direction)

    (
        feature=feature,
        family=metric_family(feature),
        n=length(values),
        n_finite=finite_count,
        artifact_finite=length(artifact_values),
        other_finite=length(other_values),
        artifact_prevalence=artifact_prevalence,
        auc_high_means_artifact=auc,
        directional_auc=directional_auc,
        rank_biserial_signed=isfinite(auc) ? 2 * auc - 1 : NaN,
        direction=direction,
        artifact_median=artifact_median,
        other_median=other_median,
        median_delta_artifact_minus_other=median_delta,
        robust_shift_iqr=robust_shift,
        artifact_q25=artifact_q25,
        artifact_q75=artifact_q75,
        other_q25=other_q25,
        other_q75=other_q75,
        patient_auc_count=patient_auc_count,
        patient_directional_auc_median=patient_auc_median,
        patient_direction_match_fraction=patient_match_fraction,
    )
end

function summarize_table(table)
    rows = [
        summarize_feature(feature, table.columns[idx], table.truth, table.patients)
        for (idx, feature) in enumerate(table.feature_names)
    ]
    sort!(
        rows;
        by = row -> (
            isfinite(row.directional_auc) ? -row.directional_auc : Inf,
            isfinite(row.patient_directional_auc_median) ? -row.patient_directional_auc_median : Inf,
            row.feature,
        )
    )
    rows
end

function best_by_family(summary_rows)
    best = Dict{String,NamedTuple}()
    for row in summary_rows
        if !haskey(best, row.family) ||
                row.directional_auc > best[row.family].directional_auc
            best[row.family] = row
        end
    end
    rows = collect(values(best))
    sort!(
        rows;
        by = row -> (
            isfinite(row.directional_auc) ? -row.directional_auc : Inf,
            row.family,
        )
    )
    rows
end

function robust_center_and_scale(values)
    finite_values = Float64[value for value in values if isfinite(value)]
    isempty(finite_values) && return 0.0, 1.0
    center = median_or_nan(finite_values)
    deviations = abs.(finite_values .- center)
    scale = 1.4826 * median_or_nan(deviations)
    if !isfinite(scale) || scale <= eps(Float64)
        scale = std(finite_values)
    end
    if !isfinite(scale) || scale <= eps(Float64)
        scale = 1.0
    end
    center, scale
end

function standardized_matrix(table; z_clip=20.0)
    n = length(table.truth)
    p = length(table.feature_names)
    X = Matrix{Float64}(undef, n, p)
    centers = Vector{Float64}(undef, p)
    scales = Vector{Float64}(undef, p)

    for feature_idx in 1:p
        values = table.columns[feature_idx]
        center, scale = robust_center_and_scale(values)
        centers[feature_idx] = center
        scales[feature_idx] = scale
        for row_idx in 1:n
            value = values[row_idx]
            standardized = isfinite(value) ? (value - center) / scale : 0.0
            if !isfinite(standardized)
                standardized = 0.0
            end
            X[row_idx, feature_idx] = clamp(standardized, -z_clip, z_clip)
        end
    end

    X, centers, scales
end

function weighted_linear_ridge(table; summary_rows, lambda=10.0)
    X, centers, scales = standardized_matrix(table)
    n, p = size(X)
    y = Float64.(table.truth)
    positives = count(table.truth)
    negatives = n - positives
    positive_weight = positives == 0 ? 1.0 : n / (2 * positives)
    negative_weight = negatives == 0 ? 1.0 : n / (2 * negatives)
    weights = [table.truth[idx] ? positive_weight : negative_weight for idx in 1:n]
    sqrt_weights = sqrt.(weights)

    design = Matrix{Float64}(undef, n, p + 1)
    design[:, 1] .= 1.0
    design[:, 2:end] .= X
    for col_idx in axes(design, 2)
        @inbounds for row_idx in axes(design, 1)
            design[row_idx, col_idx] *= sqrt_weights[row_idx]
        end
    end

    weighted_y = y .* sqrt_weights
    normal_matrix = transpose(design) * design
    for coef_idx in 2:(p + 1)
        normal_matrix[coef_idx, coef_idx] += lambda
    end
    rhs = transpose(design) * weighted_y
    coefficients = normal_matrix \ rhs

    scores = coefficients[1] .+ X * coefficients[2:end]
    score_auc = auc_from_scores(scores, table.truth)
    weighted_mse = sum(weights .* (scores .- y).^2) / sum(weights)
    summary_by_feature = Dict(row.feature => row for row in summary_rows)

    coefficient_rows = NamedTuple[]
    for feature_idx in 1:p
        feature = table.feature_names[feature_idx]
        coef = coefficients[feature_idx + 1]
        summary = summary_by_feature[feature]
        push!(coefficient_rows, (
            feature=feature,
            family=metric_family(feature),
            coefficient_standardized=coef,
            abs_coefficient_standardized=abs(coef),
            direction=coef >= 0 ? "higher_predicts_artifact" : "lower_predicts_artifact",
            center=centers[feature_idx],
            scale=scales[feature_idx],
            univariate_directional_auc=summary.directional_auc,
            univariate_direction=summary.direction,
            univariate_patient_direction_match_fraction=summary.patient_direction_match_fraction,
        ))
    end

    sort!(
        coefficient_rows;
        by = row -> (
            -row.abs_coefficient_standardized,
            row.feature,
        )
    )

    model_row = (
        n=n,
        feature_count=p,
        lambda=lambda,
        artifact_prevalence=positives / n,
        score_auc_high_means_artifact=score_auc,
        score_directional_auc=isfinite(score_auc) ? max(score_auc, 1 - score_auc) : NaN,
        weighted_mse=weighted_mse,
        intercept=coefficients[1],
    )

    coefficient_rows, model_row
end

function csv_escape(value)
    if value isa AbstractFloat
        if isnan(value)
            return "NaN"
        elseif isinf(value)
            return value > 0 ? "Inf" : "-Inf"
        else
            return @sprintf("%.12g", value)
        end
    end
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    text
end

function write_rows(path, rows; columns=SUMMARY_COLUMNS)
    open(path, "w") do io
        println(io, join(String.(columns), ","))
        for row in rows
            values = [csv_escape(getproperty(row, column)) for column in columns]
            println(io, join(values, ","))
        end
    end
end

function write_model_meta(path, row)
    columns = [
        :n,
        :feature_count,
        :lambda,
        :artifact_prevalence,
        :score_auc_high_means_artifact,
        :score_directional_auc,
        :weighted_mse,
        :intercept,
    ]
    write_rows(path, [row]; columns)
end

function markdown_table(rows; columns=[:feature, :directional_auc, :direction, :artifact_median, :other_median, :patient_direction_match_fraction], n=12)
    shown = rows[1:min(n, length(rows))]
    lines = String[]
    push!(lines, "| " * join(String.(columns), " | ") * " |")
    push!(lines, "| " * join(fill("---", length(columns)), " | ") * " |")
    for row in shown
        values = [csv_escape(getproperty(row, column)) for column in columns]
        push!(lines, "| " * join(values, " | ") * " |")
    end
    join(lines, "\n")
end

function write_meta(path; label, input_csv, table, sample_rows, seed)
    open(path, "w") do io
        println(io, "label,input_csv,loaded_rows,total_source_rows,sampled,sample_rows,seed,feature_count,artifact_rows,other_rows")
        println(io, join(csv_escape.([
            label,
            input_csv,
            length(table.truth),
            table.total_rows,
            table.sampled,
            isnothing(sample_rows) ? "" : sample_rows,
            seed,
            length(table.feature_names),
            count(table.truth),
            count(.!table.truth),
        ]), ","))
    end
end

function analyze_one(path, output_dir, label; sample_rows=nothing, seed=7)
    header = read_header(path)
    feature_indices, feature_names = metric_feature_indices(header)
    isempty(feature_names) && error("No hvg_ feature columns found in $(path)")

    @info "Loading metric table" label path features=length(feature_names) sample_rows
    table = isnothing(sample_rows) ?
        load_metric_columns(path, feature_indices, feature_names) :
        sample_metric_columns(path, feature_indices, feature_names; sample_rows, seed)

    @info "Summarizing artifact associations" label loaded_rows=length(table.truth)
    summary_rows = summarize_table(table)
    family_rows = best_by_family(summary_rows)
    @info "Fitting joint ridge linear probability model" label loaded_rows=length(table.truth)
    regression_rows, regression_meta = weighted_linear_ridge(table; summary_rows)
    suffix = isnothing(sample_rows) ? "exact" : "sampled"

    write_rows(joinpath(output_dir, "$(label)_artifact_metric_summary_$(suffix).csv"), summary_rows)
    write_rows(joinpath(output_dir, "$(label)_artifact_metric_family_best_$(suffix).csv"), family_rows)
    write_rows(
        joinpath(output_dir, "$(label)_artifact_metric_ridge_linear_coefficients_$(suffix).csv"),
        regression_rows;
        columns=REGRESSION_COLUMNS,
    )
    write_model_meta(
        joinpath(output_dir, "$(label)_artifact_metric_ridge_linear_model_$(suffix).csv"),
        regression_meta,
    )
    write_meta(
        joinpath(output_dir, "$(label)_artifact_metric_summary_$(suffix)_meta.csv");
        label,
        input_csv=path,
        table,
        sample_rows,
        seed,
    )
    summary_rows, family_rows, regression_rows, regression_meta, table
end

function main()
    args = parse_args(ARGS)
    if haskey(args, "help") || !haskey(args, "whole-hvg-csv")
        usage()
        return
    end

    output_dir = get(
        args,
        "output-dir",
        joinpath("data", "exp_pro", "artifact_metric_relationships", string(Dates.now())),
    )
    mkpath(output_dir)
    seed = parse(Int, get(args, "seed", "7"))
    large_sample_rows = parse(Int, get(args, "large-sample-rows", "150000"))

    exact_rows, exact_family_rows, exact_regression_rows, exact_regression_meta, _ = analyze_one(
        args["whole-hvg-csv"],
        output_dir,
        "whole_hvg_15s";
        seed,
    )

    large_rows = nothing
    large_family_rows = nothing
    large_regression_rows = nothing
    large_regression_meta = nothing
    if haskey(args, "whole-hvg-large-csv")
        large_rows, large_family_rows, large_regression_rows, large_regression_meta, _ = analyze_one(
            args["whole-hvg-large-csv"],
            output_dir,
            "whole_hvg_1s";
            sample_rows=large_sample_rows,
            seed,
        )
    end

    readme_path = joinpath(output_dir, "artifact_metric_relationships_readme.md")
    open(readme_path, "w") do io
        println(io, "# Artifact Metric Relationship Analysis")
        println(io)
        println(io, "- Generated: $(Dates.now())")
        println(io, "- Whole-HVG exact input: `$(args["whole-hvg-csv"])`")
        if haskey(args, "whole-hvg-large-csv")
            println(io, "- Whole-HVG sampled input: `$(args["whole-hvg-large-csv"])`")
            println(io, "- 1 s sample rows: $(large_sample_rows)")
        end
        println(io)
        println(io, "## Top Whole-HVG 15 s Metrics")
        println(io)
        println(io, markdown_table(exact_rows))
        println(io)
        println(io, "## Top Whole-HVG 15 s Metric Families")
        println(io)
        println(io, markdown_table(exact_family_rows))
        println(io)
        println(io, "## Joint Ridge Linear Model, 15 s")
        println(io)
        println(io, "- Directional score AUC: $(csv_escape(exact_regression_meta.score_directional_auc))")
        println(io, "- Weighted MSE: $(csv_escape(exact_regression_meta.weighted_mse))")
        println(io)
        println(io, markdown_table(
            exact_regression_rows;
            columns=[
                :feature,
                :coefficient_standardized,
                :direction,
                :univariate_directional_auc,
                :univariate_patient_direction_match_fraction,
            ],
        ))
        if !isnothing(large_rows)
            println(io)
            println(io, "## Top Whole-HVG 1 s Metrics")
            println(io)
            println(io, markdown_table(large_rows))
            println(io)
            println(io, "## Top Whole-HVG 1 s Metric Families")
            println(io)
            println(io, markdown_table(large_family_rows))
            println(io)
            println(io, "## Joint Ridge Linear Model, 1 s Sample")
            println(io)
            println(io, "- Directional score AUC: $(csv_escape(large_regression_meta.score_directional_auc))")
            println(io, "- Weighted MSE: $(csv_escape(large_regression_meta.weighted_mse))")
            println(io)
            println(io, markdown_table(
                large_regression_rows;
                columns=[
                    :feature,
                    :coefficient_standardized,
                    :direction,
                    :univariate_directional_auc,
                    :univariate_patient_direction_match_fraction,
                ],
            ))
        end
    end

    println("Wrote artifact metric relationship analysis to $(output_dir)")
    println()
    println("Top 15 s whole-HVG metrics:")
    println(markdown_table(exact_rows; n=8))
    if !isnothing(large_rows)
        println()
        println("Top sampled 1 s whole-HVG metrics:")
        println(markdown_table(large_rows; n=8))
    end
    println()
    println("Top 15 s ridge linear coefficients:")
    println(markdown_table(
        exact_regression_rows;
        columns=[
            :feature,
            :coefficient_standardized,
            :direction,
            :univariate_directional_auc,
            :univariate_patient_direction_match_fraction,
        ],
        n=8,
    ))
    if !isnothing(large_regression_rows)
        println()
        println("Top sampled 1 s ridge linear coefficients:")
        println(markdown_table(
            large_regression_rows;
            columns=[
                :feature,
                :coefficient_standardized,
                :direction,
                :univariate_directional_auc,
                :univariate_patient_direction_match_fraction,
            ],
            n=8,
        ))
    end
end

main()
