using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using Random
using Statistics

const STABLE_ARTIFACT_FEATURES = Symbol[
    :hvg_bin_edge_span_mean,
    :hvg_bin_edge_span_q95,
    :hvg_degree_imbalance_std,
    :hvg_degree_abs_robust_z_mean,
    :hvg_random_law_surprisal_abs_robust_z_mean,
    :hvg_bin_global_clustering_coefficient,
]

const ID_COLUMNS = Symbol[
    :patient,
    :channel,
    :channel_label,
    :bin_index,
    :bin_start,
    :bin_stop,
    :artifact_truth,
    :seizure_truth,
    :clean_seizure_truth,
    :seizure_event_id,
]

function parse_cli_args(args)
    parsed = Dict{String,String}()
    idx = 1
    while idx <= length(args)
        arg = args[idx]
        startswith(arg, "--") || throw(ArgumentError("Unexpected positional argument: $(arg)"))
        key = arg[3:end]
        idx == length(args) && throw(ArgumentError("Missing value for $(arg)"))
        parsed[key] = args[idx + 1]
        idx += 2
    end
    parsed
end

function parse_int_list(text)
    [parse(Int, strip(part)) for part in split(text, ',') if !isempty(strip(part))]
end

function robust_center_and_scale(values)
    finite_values = Float64[Float64(value) for value in values if !ismissing(value) && isfinite(Float64(value))]
    isempty(finite_values) && return 0.0, 1.0
    center = median(finite_values)
    scale = 1.4826 * median(abs.(finite_values .- center))
    if !isfinite(scale) || scale <= eps(Float64)
        scale = std(finite_values)
    end
    if !isfinite(scale) || scale <= eps(Float64)
        scale = 1.0
    end
    center, scale
end

function selected_input_columns(feature_names)
    selected = Set(String.(feature_names))
    union!(selected, String.(ID_COLUMNS))
    (_, name) -> String(name) in selected
end

function standardized_matrix(df::AbstractDataFrame, feature_names; z_clip=20.0)
    n = nrow(df)
    p = length(feature_names)
    Z = Matrix{Float64}(undef, n, p)
    centers = Float64[]
    scales = Float64[]

    for (feature_idx, feature_name) in enumerate(feature_names)
        values = df[!, feature_name]
        center, scale = robust_center_and_scale(values)
        push!(centers, center)
        push!(scales, scale)
        for row_idx in 1:n
            value = values[row_idx]
            value = ismissing(value) ? NaN : Float64(value)
            standardized = isfinite(value) ? (value - center) / scale : 0.0
            if !isfinite(standardized)
                standardized = 0.0
            end
            Z[row_idx, feature_idx] = clamp(standardized, -z_clip, z_clip)
        end
    end

    Z, centers, scales
end

sqdist_row_center(X, row_idx, centers, center_idx) =
    sum((X[row_idx, feature_idx] - centers[center_idx, feature_idx])^2 for feature_idx in axes(X, 2))

function weighted_index_sample(rng, weights)
    total = sum(weights)
    if !isfinite(total) || total <= 0
        return rand(rng, eachindex(weights))
    end
    threshold = rand(rng) * total
    running = 0.0
    for idx in eachindex(weights)
        running += weights[idx]
        if running >= threshold
            return idx
        end
    end
    lastindex(weights)
end

function kmeanspp_init(rng, X, k)
    n, p = size(X)
    centers = Matrix{Float64}(undef, k, p)
    first_idx = rand(rng, 1:n)
    centers[1, :] .= X[first_idx, :]
    min_d2 = [sqdist_row_center(X, row_idx, centers, 1) for row_idx in 1:n]

    for center_idx in 2:k
        row_idx = weighted_index_sample(rng, min_d2)
        centers[center_idx, :] .= X[row_idx, :]
        for idx in 1:n
            min_d2[idx] = min(min_d2[idx], sqdist_row_center(X, idx, centers, center_idx))
        end
    end

    centers
end

function assign_clusters!(assignments, distances, X, centers)
    changed = 0
    total_sse = 0.0
    for row_idx in axes(X, 1)
        best_cluster = 1
        best_dist = sqdist_row_center(X, row_idx, centers, 1)
        for cluster in 2:size(centers, 1)
            dist = sqdist_row_center(X, row_idx, centers, cluster)
            if dist < best_dist
                best_dist = dist
                best_cluster = cluster
            end
        end
        if assignments[row_idx] != best_cluster
            changed += 1
            assignments[row_idx] = best_cluster
        end
        distances[row_idx] = best_dist
        total_sse += best_dist
    end
    changed, total_sse
end

function update_centers!(rng, centers, X, assignments)
    fill!(centers, 0.0)
    counts = zeros(Int, size(centers, 1))
    for row_idx in axes(X, 1)
        cluster = assignments[row_idx]
        counts[cluster] += 1
        for feature_idx in axes(X, 2)
            centers[cluster, feature_idx] += X[row_idx, feature_idx]
        end
    end
    for cluster in axes(centers, 1)
        if counts[cluster] == 0
            centers[cluster, :] .= X[rand(rng, axes(X, 1)), :]
        else
            centers[cluster, :] ./= counts[cluster]
        end
    end
    counts
end

function fit_kmeans(rng, X, k; max_iter=60, restarts=3, tol_changed=0)
    best = nothing
    for restart in 1:restarts
        centers = kmeanspp_init(rng, X, k)
        assignments = zeros(Int, size(X, 1))
        distances = zeros(Float64, size(X, 1))
        sse = Inf
        counts = zeros(Int, k)
        iterations = 0
        for iter in 1:max_iter
            changed, sse = assign_clusters!(assignments, distances, X, centers)
            counts = update_centers!(rng, centers, X, assignments)
            iterations = iter
            changed <= tol_changed && break
        end
        _, sse = assign_clusters!(assignments, distances, X, centers)
        if isnothing(best) || sse < best.sse
            best = (
                centers=copy(centers),
                sample_assignments=copy(assignments),
                sample_distances=copy(distances),
                sample_counts=copy(counts),
                sse=sse,
                iterations=iterations,
                restart=restart,
            )
        end
    end
    best
end

function assign_all(X, centers)
    assignments = zeros(Int, size(X, 1))
    distances = zeros(Float64, size(X, 1))
    _, sse = assign_clusters!(assignments, distances, X, centers)
    assignments, distances, sse
end

choose_fit_indices(rng, n, fit_rows) =
    n <= fit_rows ? collect(1:n) : randperm(rng, n)[1:fit_rows]

function comb2(x)
    x < 2 && return 0.0
    x * (x - 1) / 2
end

function adjusted_rand_index(assignments, truth, k)
    n = length(assignments)
    cluster_counts = zeros(Int, k)
    label_counts = zeros(Int, 2)
    contingency = zeros(Int, k, 2)
    for idx in 1:n
        cluster = assignments[idx]
        label = truth[idx] ? 2 : 1
        cluster_counts[cluster] += 1
        label_counts[label] += 1
        contingency[cluster, label] += 1
    end

    index = sum(comb2, contingency)
    cluster_pairs = sum(comb2, cluster_counts)
    label_pairs = sum(comb2, label_counts)
    total_pairs = comb2(n)
    expected = cluster_pairs * label_pairs / total_pairs
    max_index = (cluster_pairs + label_pairs) / 2
    denom = max_index - expected
    abs(denom) <= eps(Float64) ? 0.0 : (index - expected) / denom
end

function normalized_mutual_information(assignments, truth, k)
    n = length(assignments)
    cluster_counts = zeros(Int, k)
    label_counts = zeros(Int, 2)
    contingency = zeros(Int, k, 2)
    for idx in 1:n
        cluster = assignments[idx]
        label = truth[idx] ? 2 : 1
        cluster_counts[cluster] += 1
        label_counts[label] += 1
        contingency[cluster, label] += 1
    end

    mutual_information = 0.0
    for cluster in 1:k, label in 1:2
        count = contingency[cluster, label]
        count == 0 && continue
        mutual_information += count / n * log((count * n) / (cluster_counts[cluster] * label_counts[label]))
    end
    h_cluster = -sum(count == 0 ? 0.0 : count / n * log(count / n) for count in cluster_counts)
    h_label = -sum(count == 0 ? 0.0 : count / n * log(count / n) for count in label_counts)
    denom = sqrt(h_cluster * h_label)
    denom <= eps(Float64) ? 0.0 : mutual_information / denom
end

function confusion_metrics(predicted, truth)
    tp = count(predicted .& truth)
    fp = count(predicted .& .!truth)
    tn = count(.!predicted .& .!truth)
    fn = count(.!predicted .& truth)
    precision = tp + fp == 0 ? NaN : tp / (tp + fp)
    recall = tp + fn == 0 ? NaN : tp / (tp + fn)
    specificity = tn + fp == 0 ? NaN : tn / (tn + fp)
    f1 = isfinite(precision) && isfinite(recall) && precision + recall > 0 ?
        2 * precision * recall / (precision + recall) :
        NaN
    balanced_accuracy = isfinite(recall) && isfinite(specificity) ?
        (recall + specificity) / 2 :
        NaN
    (
        true_positives=tp,
        false_positives=fp,
        true_negatives=tn,
        false_negatives=fn,
        precision=precision,
        recall=recall,
        specificity=specificity,
        f1=f1,
        balanced_accuracy=balanced_accuracy,
        predicted_artifact_fraction=(tp + fp) / length(truth),
    )
end

function cluster_detail_rows(assignments, truth, distances, k, artifact_prevalence)
    rows = NamedTuple[]
    for cluster in 1:k
        mask = assignments .== cluster
        n_cluster = count(mask)
        artifacts = count(truth[mask])
        rate = n_cluster == 0 ? NaN : artifacts / n_cluster
        push!(rows, (
            cluster=cluster,
            n=n_cluster,
            fraction=n_cluster / length(assignments),
            artifact_count=artifacts,
            other_count=n_cluster - artifacts,
            artifact_rate=rate,
            artifact_enrichment=rate / artifact_prevalence,
            mean_squared_distance=n_cluster == 0 ? NaN : mean(distances[mask]),
        ))
    end
    sort!(rows; by=row -> (-row.artifact_rate, -row.n, row.cluster))
    rows
end

function best_artifact_cluster_selection(assignments, truth, detail_rows)
    best = nothing
    selected = Set{Int}()
    for row in detail_rows
        push!(selected, row.cluster)
        predicted = [assignment in selected for assignment in assignments]
        metrics = confusion_metrics(predicted, truth)
        candidate = merge((
            selected_cluster_count=length(selected),
            selected_clusters=join(sort(collect(selected)), ";"),
            selected_artifact_rate_minimum=row.artifact_rate,
        ), metrics)
        if isnothing(best) || candidate.f1 > best.f1
            best = candidate
        end
    end
    best
end

function write_feature_scaling(path, feature_names, centers, scales)
    DataFrame(
        feature=String.(feature_names),
        robust_center=centers,
        robust_scale=scales,
    ) |> df -> CSV.write(path, df)
end

function write_readme(path; input_csv, output_dir, feature_names, summary_df, best_details_df, fit_rows, restarts, max_iter)
    open(path, "w") do io
        println(io, "# Stable Whole-HVG Artifact Clustering")
        println(io)
        println(io, "- Generated: $(Dates.now())")
        println(io, "- Input CSV: `$(input_csv)`")
        println(io, "- Output dir: `$(output_dir)`")
        println(io, "- Fit sample rows: $(fit_rows)")
        println(io, "- Restarts: $(restarts)")
        println(io, "- Max iterations: $(max_iter)")
        println(io)
        println(io, "## Features")
        println(io)
        for feature in feature_names
            println(io, "- `$(feature)`")
        end
        println(io)
        println(io, "## Summary")
        println(io)
        println(io, "| k | ARI | NMI | best F1 | precision | recall | selected clusters | max cluster artifact rate |")
        println(io, "| ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |")
        for row in eachrow(summary_df)
            println(io,
                "| $(row.k) | $(@sprintf("%.4f", row.adjusted_rand_index)) | " *
                "$(@sprintf("%.4f", row.normalized_mutual_information)) | " *
                "$(@sprintf("%.4f", row.best_f1)) | $(@sprintf("%.4f", row.best_precision)) | " *
                "$(@sprintf("%.4f", row.best_recall)) | $(row.best_selected_clusters) | " *
                "$(@sprintf("%.4f", row.max_cluster_artifact_rate)) |"
            )
        end
        println(io)
        println(io, "## Most Artifact-Enriched Clusters")
        println(io)
        println(io, "| k | cluster | n | artifact rate | enrichment |")
        println(io, "| ---: | ---: | ---: | ---: | ---: |")
        for row in eachrow(best_details_df)
            println(io,
                "| $(row.k) | $(row.cluster) | $(row.n) | " *
                "$(@sprintf("%.4f", row.artifact_rate)) | $(@sprintf("%.2f", row.artifact_enrichment)) |"
            )
        end
    end
end

function main()
    args = parse_cli_args(ARGS)
    input_csv = get(args, "input-csv", "data/exp_pro/whole_hvg_anomaly/2026-06-08T19:57:17.670/whole_hvg_bin_scores.csv")
    output_dir = get(args, "output-dir", joinpath("data", "exp_pro", "artifact_clustering", string(Dates.now())))
    fit_rows = parse(Int, get(args, "fit-rows", "100000"))
    seed = parse(Int, get(args, "seed", "11"))
    restarts = parse(Int, get(args, "restarts", "3"))
    max_iter = parse(Int, get(args, "max-iter", "60"))
    ks = parse_int_list(get(args, "ks", "2,3,4,6,8,12,16"))
    feature_names = isdefined(Main, :ARTIFACT_CLUSTER_FEATURES) ?
        Symbol.(collect(ARTIFACT_CLUSTER_FEATURES)) :
        STABLE_ARTIFACT_FEATURES

    mkpath(output_dir)
    @info "Loading clustering input" input_csv features=feature_names
    df = CSV.read(input_csv, DataFrame; select=selected_input_columns(feature_names))
    truth = Bool.(df.artifact_truth)
    artifact_prevalence = count(truth) / length(truth)
    @info "Standardizing features" rows=nrow(df) artifact_prevalence
    Z, centers, scales = standardized_matrix(df, feature_names)
    write_feature_scaling(joinpath(output_dir, "clustering_feature_scaling.csv"), feature_names, centers, scales)

    rng = MersenneTwister(seed)
    fit_indices = choose_fit_indices(rng, nrow(df), fit_rows)
    X_fit = Z[fit_indices, :]
    summary_rows = NamedTuple[]
    all_detail_rows = NamedTuple[]
    best_detail_rows = NamedTuple[]

    for k in ks
        @info "Fitting k-means" k fit_rows=length(fit_indices) restarts max_iter
        model = fit_kmeans(rng, X_fit, k; max_iter=max_iter, restarts=restarts)
        @info "Assigning full table to k-means centers" k rows=nrow(df)
        assignments, distances, full_sse = assign_all(Z, model.centers)
        detail_rows = cluster_detail_rows(assignments, truth, distances, k, artifact_prevalence)
        best = best_artifact_cluster_selection(assignments, truth, detail_rows)
        ari = adjusted_rand_index(assignments, truth, k)
        nmi = normalized_mutual_information(assignments, truth, k)
        max_rate = maximum(row.artifact_rate for row in detail_rows)
        min_rate = minimum(row.artifact_rate for row in detail_rows)
        push!(summary_rows, (
            k=k,
            rows=nrow(df),
            artifact_prevalence=artifact_prevalence,
            fit_rows=length(fit_indices),
            restarts=restarts,
            max_iter=max_iter,
            sample_sse=model.sse,
            full_sse=full_sse,
            adjusted_rand_index=ari,
            normalized_mutual_information=nmi,
            min_cluster_artifact_rate=min_rate,
            max_cluster_artifact_rate=max_rate,
            best_selected_cluster_count=best.selected_cluster_count,
            best_selected_clusters=best.selected_clusters,
            best_selected_artifact_rate_minimum=best.selected_artifact_rate_minimum,
            best_precision=best.precision,
            best_recall=best.recall,
            best_specificity=best.specificity,
            best_f1=best.f1,
            best_balanced_accuracy=best.balanced_accuracy,
            best_predicted_artifact_fraction=best.predicted_artifact_fraction,
        ))
        for row in detail_rows
            push!(all_detail_rows, merge((k=k,), row))
        end
        for row in detail_rows[1:min(5, length(detail_rows))]
            push!(best_detail_rows, merge((k=k,), row))
        end
    end

    summary_df = DataFrame(summary_rows)
    sort!(summary_df, [:best_f1, :normalized_mutual_information]; rev=true)
    details_df = DataFrame(all_detail_rows)
    enriched_df = DataFrame(best_detail_rows)
    CSV.write(joinpath(output_dir, "kmeans_artifact_capture_summary.csv"), summary_df)
    CSV.write(joinpath(output_dir, "kmeans_cluster_artifact_details.csv"), details_df)
    CSV.write(joinpath(output_dir, "kmeans_top_artifact_enriched_clusters.csv"), enriched_df)
    write_readme(
        joinpath(output_dir, "artifact_clustering_readme.md");
        input_csv,
        output_dir,
        feature_names,
        summary_df,
        best_details_df=enriched_df,
        fit_rows=length(fit_indices),
        restarts,
        max_iter,
    )

    println("Wrote artifact clustering analysis to $(output_dir)")
    show(stdout, summary_df; allrows=true, allcols=true)
    println()
end

main()
