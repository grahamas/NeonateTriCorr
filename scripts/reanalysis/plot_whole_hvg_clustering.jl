using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics
using TSne

include(scriptsdir("include_src.jl"))

const WHOLE_HVG_CLUSTER_ID_COLUMNS = Set([
    :patient,
    :channel,
    :channel_label,
    :bin_index,
    :bin_start,
    :bin_stop,
    :sample_start_index,
    :sample_stop_index,
    :sample_count,
    :artifact_truth,
    :seizure_truth,
    :clean_seizure_truth,
    :seizure_event_id
])

function newest_whole_hvg_bin_scores(; root=datadir("exp_pro", "whole_hvg_anomaly"))
    isdir(root) || return nothing
    candidates = String[]
    for (walk_root, _, files) in walkdir(root)
        if "whole_hvg_bin_scores.csv" in files
            push!(candidates, joinpath(walk_root, "whole_hvg_bin_scores.csv"))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=path -> stat(path).mtime)
    candidates[end]
end

function numeric_hvg_feature_names(df::AbstractDataFrame)
    names_out = Symbol[]
    for name in Symbol.(names(df))
        name in WHOLE_HVG_CLUSTER_ID_COLUMNS && continue
        startswith(String(name), "hvg_") || continue
        eltype(df[!, name]) <: Union{Missing,Real} || continue
        push!(names_out, name)
    end
    sort(names_out; by=String)
end

function whole_hvg_clustering_feature_matrix(df::AbstractDataFrame, feature_names)
    X = Matrix{Float64}(undef, nrow(df), length(feature_names))
    for (feature_idx, feature_name) in enumerate(feature_names)
        values = df[!, feature_name]
        for row_idx in 1:nrow(df)
            value = values[row_idx]
            X[row_idx, feature_idx] = artifact_feature_value(value)
        end
    end
    X
end

function whole_hvg_standardized_feature_matrix(X::AbstractMatrix, feature_names; min_scale=1e-8, z_clip=20.0)
    keep_idxs = Int[]
    centers = Float64[]
    scales = Float64[]
    kept_features = Symbol[]
    skipped_rows = NamedTuple[]

    for feature_idx in axes(X, 2)
        finite_values = Float64[x for x in view(X, :, feature_idx) if isfinite(x)]
        if isempty(finite_values)
            push!(skipped_rows, (
                feature=String(feature_names[feature_idx]),
                reason="no finite values"
            ))
            continue
        end

        center, scale = robust_center_and_scale(finite_values)
        isfinite(center) || (center = 0.0)
        if !isfinite(scale) || scale <= min_scale
            push!(skipped_rows, (
                feature=String(feature_names[feature_idx]),
                reason="near-constant robust scale ($(scale))"
            ))
            continue
        end

        push!(keep_idxs, feature_idx)
        push!(centers, center)
        push!(scales, scale)
        push!(kept_features, feature_names[feature_idx])
    end

    isempty(keep_idxs) && throw(ArgumentError("No finite, nonconstant HVG feature columns available for clustering."))

    Z = Matrix{Float64}(undef, size(X, 1), length(keep_idxs))
    for (out_idx, feature_idx) in enumerate(keep_idxs)
        center = centers[out_idx]
        scale = scales[out_idx]
        for row_idx in axes(X, 1)
            value = X[row_idx, feature_idx]
            value = isfinite(value) ? value : center
            standardized = (value - center) / scale
            if !isfinite(standardized)
                standardized = 0.0
            end
            Z[row_idx, out_idx] = clamp(standardized, -z_clip, z_clip)
        end
    end

    (
        Z=Z,
        feature_names=kept_features,
        centers=centers,
        scales=scales,
        skipped_df=DataFrame(skipped_rows)
    )
end

function pca_embedding(Z::AbstractMatrix; outdim=2)
    outdim >= 1 || throw(ArgumentError("outdim must be positive"))
    X = Matrix{Float64}(Z)
    X .-= mean(X; dims=1)
    decomposition = svd(X; full=false)
    dims = min(outdim, length(decomposition.S))
    coords = decomposition.U[:, 1:dims] * Diagonal(decomposition.S[1:dims])
    if dims < outdim
        coords = hcat(coords, zeros(Float64, size(coords, 1), outdim - dims))
    end
    coords
end

function whole_hvg_embedding_sample_indices(
        df::AbstractDataFrame;
        max_rows,
        seed,
        strata_columns=[:artifact_truth, :seizure_truth]
    )
    n_rows = nrow(df)
    n_rows <= max_rows && return collect(1:n_rows)

    rng = MersenneTwister(seed)
    strata = Dict{Tuple,Vector{Int}}()
    for row_idx in 1:n_rows
        key = Tuple(df[row_idx, column] for column in strata_columns)
        push!(get!(strata, key, Int[]), row_idx)
    end

    selected = Int[]
    sorted_strata = sort(collect(strata); by=pair -> string(pair[1]))
    remaining_budget = max_rows
    remaining_strata = length(sorted_strata)
    for (_, idxs) in sorted_strata
        target = max(1, floor(Int, max_rows * length(idxs) / n_rows))
        target = min(target, length(idxs), remaining_budget - remaining_strata + 1)
        if target > 0
            shuffled = shuffle(rng, idxs)
            append!(selected, shuffled[1:target])
            remaining_budget -= target
        end
        remaining_strata -= 1
    end

    if length(selected) < max_rows
        selected_set = Set(selected)
        remaining = [idx for idx in 1:n_rows if idx ∉ selected_set]
        shuffled = shuffle(rng, remaining)
        append!(selected, shuffled[1:min(max_rows - length(selected), length(shuffled))])
    elseif length(selected) > max_rows
        selected = shuffle(rng, selected)[1:max_rows]
    end

    sort!(unique(selected))
end

function tsne_embedding(
        Z::AbstractMatrix;
        outdim=2,
        reduce_dims=min(30, size(Z, 2)),
        iterations=1000,
        perplexity=30.0,
        seed=1,
        pca_init=true
    )
    Random.seed!(seed)
    tsne(Matrix{Float64}(Z), outdim, reduce_dims, iterations, perplexity;
        pca_init=pca_init,
        progress=false,
        verbose=false
    )
end

function whole_hvg_distribution_values(df::AbstractDataFrame, feature_name::Symbol)
    values = df[!, feature_name]
    finite_values = Float64[]
    invalid_count = 0
    for value in values
        feature_value = artifact_feature_value(value)
        if isfinite(feature_value)
            push!(finite_values, feature_value)
        else
            invalid_count += 1
        end
    end
    finite_values, invalid_count
end

function whole_hvg_feature_distribution_rows(df::AbstractDataFrame, feature_names, centers, scales)
    rows = NamedTuple[]
    for (idx, feature_name) in enumerate(feature_names)
        finite_values, invalid_count = whole_hvg_distribution_values(df, feature_name)
        if isempty(finite_values)
            push!(rows, (
                feature=String(feature_name),
                n_finite=0,
                n_invalid=invalid_count,
                min=NaN,
                p05=NaN,
                q25=NaN,
                median=NaN,
                mean=NaN,
                q75=NaN,
                p95=NaN,
                max=NaN,
                std=NaN,
                robust_center=centers[idx],
                robust_scale=scales[idx]
            ))
        else
            push!(rows, (
                feature=String(feature_name),
                n_finite=length(finite_values),
                n_invalid=invalid_count,
                min=minimum(finite_values),
                p05=quantile(finite_values, 0.05),
                q25=quantile(finite_values, 0.25),
                median=median(finite_values),
                mean=mean(finite_values),
                q75=quantile(finite_values, 0.75),
                p95=quantile(finite_values, 0.95),
                max=maximum(finite_values),
                std=length(finite_values) > 1 ? std(finite_values) : 0.0,
                robust_center=centers[idx],
                robust_scale=scales[idx]
            ))
        end
    end
    DataFrame(rows)
end

function plain_float(value)
    if value isa Real && isfinite(value)
        return string(round(Float64(value); sigdigits=6))
    end
    string(value)
end

function write_whole_hvg_feature_distribution_report(
        path,
        distribution_df::AbstractDataFrame,
        skipped_df::AbstractDataFrame;
        method,
        sample_df=DataFrame()
    )
    open(path, "w") do io
        println(io, "Whole-HVG $(uppercase(String(method))) Feature Report")
        println(io, "Generated: $(Dates.now())")
        println(io)
        if nrow(sample_df) > 0
            sample = first(eachrow(sample_df))
            println(io, "Embedding rows: $(sample.embedded_rows) / $(sample.source_rows)")
            if method == :tsne
                println(io, "t-SNE max rows: $(sample.tsne_max_rows)")
                println(io, "t-SNE seed: $(sample.tsne_seed)")
                println(io, "t-SNE iterations: $(sample.tsne_iterations)")
                println(io, "t-SNE perplexity: $(sample.tsne_perplexity)")
            end
            println(io)
        end
        println(io, "Features used: $(nrow(distribution_df))")
        println(io, "Features skipped: $(nrow(skipped_df))")
        println(io)

        if nrow(skipped_df) > 0
            println(io, "Skipped features")
            for row in eachrow(skipped_df)
                println(io, "- $(row.feature): $(row.reason)")
            end
            println(io)
        end

        println(io, "Used feature distributions")
        for row in eachrow(distribution_df)
            println(io, "- $(row.feature)")
            println(io,
                "  n=$(row.n_finite), invalid=$(row.n_invalid), " *
                "min=$(plain_float(row.min)), p05=$(plain_float(row.p05)), " *
                "q25=$(plain_float(row.q25)), median=$(plain_float(row.median)), " *
                "mean=$(plain_float(row.mean)), q75=$(plain_float(row.q75)), " *
                "p95=$(plain_float(row.p95)), max=$(plain_float(row.max)), " *
                "std=$(plain_float(row.std)), robust_center=$(plain_float(row.robust_center)), " *
                "robust_scale=$(plain_float(row.robust_scale))"
            )
        end
    end
end

function whole_hvg_embedding_df(
        df::AbstractDataFrame;
        method=:pca,
        tsne_max_rows=5000,
        tsne_seed=1,
        tsne_iterations=1000,
        tsne_perplexity=30.0
    )
    feature_names = numeric_hvg_feature_names(df)
    raw_X = whole_hvg_clustering_feature_matrix(df, feature_names)
    standardized = whole_hvg_standardized_feature_matrix(raw_X, feature_names)

    sample_idxs = if method == :tsne
        whole_hvg_embedding_sample_indices(df;
            max_rows=tsne_max_rows,
            seed=tsne_seed
        )
    else
        collect(1:nrow(df))
    end

    Z_embedding = standardized.Z[sample_idxs, :]
    coords = if method == :pca
        pca_embedding(Z_embedding; outdim=2)
    elseif method == :tsne
        tsne_embedding(Z_embedding;
            outdim=2,
            iterations=tsne_iterations,
            perplexity=tsne_perplexity,
            seed=tsne_seed
        )
    else
        throw(ArgumentError("Unsupported whole-HVG clustering method: $(method)"))
    end

    embedding_df = select(
        DataFrame(df[sample_idxs, :]),
        collect(WHOLE_HVG_CLUSTER_ID_COLUMNS ∩ Set(Symbol.(names(df))))...
    )
    embedding_df.embedding_x = coords[:, 1]
    embedding_df.embedding_y = coords[:, 2]
    embedding_df.embedding_method = fill(String(method), nrow(embedding_df))
    embedding_df.source_row_index = sample_idxs

    (
        embedding_df=embedding_df,
        feature_names=standardized.feature_names,
        centers=standardized.centers,
        scales=standardized.scales,
        skipped_df=standardized.skipped_df,
        sample_df=DataFrame(
            method=String(method),
            source_rows=nrow(df),
            embedded_rows=length(sample_idxs),
            tsne_max_rows=method == :tsne ? tsne_max_rows : missing,
            tsne_seed=method == :tsne ? tsne_seed : missing,
            tsne_iterations=method == :tsne ? tsne_iterations : missing,
            tsne_perplexity=method == :tsne ? tsne_perplexity : missing
        ),
        distribution_df=whole_hvg_feature_distribution_rows(
            df,
            standardized.feature_names,
            standardized.centers,
            standardized.scales
        )
    )
end

function label_colors(labels; false_color=(:gray60, 0.35), true_color=(:red, 0.75))
    [Bool(label) ? true_color : false_color for label in labels]
end

function plot_whole_hvg_embedding(
        embedding_df::AbstractDataFrame;
        color_column,
        title,
        true_label,
        false_label="other",
        true_color=:red,
        resolution=(1100, 850)
    )
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1];
        xlabel="embedding 1",
        ylabel="embedding 2",
        title=title
    )

    colors = label_colors(embedding_df[!, color_column];
        true_color=(true_color, 0.8),
        false_color=(:gray70, 0.25)
    )
    scatter!(ax, embedding_df.embedding_x, embedding_df.embedding_y;
        color=colors,
        markersize=4,
        strokewidth=0
    )
    scatter!(ax, [NaN], [NaN]; color=(:gray70, 0.55), markersize=10, label=false_label)
    scatter!(ax, [NaN], [NaN]; color=true_color, markersize=10, label=true_label)
    axislegend(ax; position=:rt)
    fig
end

function save_whole_hvg_clustering_outputs(
        bin_df::AbstractDataFrame;
        output_dir=datadir("exp_pro", "whole_hvg_clustering", string(Dates.now())),
        plot_dir=plotsdir("whole_hvg_clustering_$(Dates.now())"),
        method=:pca,
        tsne_max_rows=5000,
        tsne_seed=1,
        tsne_iterations=1000,
        tsne_perplexity=30.0
    )
    mkpath(output_dir)
    mkpath(plot_dir)

    result = whole_hvg_embedding_df(bin_df;
        method=method,
        tsne_max_rows=tsne_max_rows,
        tsne_seed=tsne_seed,
        tsne_iterations=tsne_iterations,
        tsne_perplexity=tsne_perplexity
    )
    embedding_df = result.embedding_df
    CSV.write(joinpath(output_dir, "whole_hvg_embedding.csv"), embedding_df)
    CSV.write(joinpath(output_dir, "whole_hvg_embedding_sample.csv"), result.sample_df)
    CSV.write(
        joinpath(output_dir, "whole_hvg_embedding_features.csv"),
        DataFrame(feature=String.(result.feature_names), center=result.centers, scale=result.scales)
    )
    CSV.write(joinpath(output_dir, "whole_hvg_embedding_skipped_features.csv"), result.skipped_df)
    CSV.write(joinpath(output_dir, "whole_hvg_embedding_feature_distributions.csv"), result.distribution_df)
    write_whole_hvg_feature_distribution_report(
        joinpath(output_dir, "whole_hvg_embedding_feature_report.txt"),
        result.distribution_df,
        result.skipped_df;
        method=method,
        sample_df=result.sample_df
    )

    artifact_fig = plot_whole_hvg_embedding(embedding_df;
        color_column=:artifact_truth,
        title="Whole-recording HVG $(uppercase(String(method))) embedding by artifact truth",
        true_label="artifact",
        true_color=:red
    )
    save(joinpath(plot_dir, "whole_hvg_embedding_artifact_truth.png"), artifact_fig)

    seizure_fig = plot_whole_hvg_embedding(embedding_df;
        color_column=:seizure_truth,
        title="Whole-recording HVG $(uppercase(String(method))) embedding by seizure truth",
        true_label="seizure",
        true_color=:dodgerblue
    )
    save(joinpath(plot_dir, "whole_hvg_embedding_seizure_truth.png"), seizure_fig)

    if :clean_seizure_truth in Symbol.(names(embedding_df))
        clean_seizure_fig = plot_whole_hvg_embedding(embedding_df;
            color_column=:clean_seizure_truth,
            title="Whole-recording HVG $(uppercase(String(method))) embedding by clean seizure truth",
            true_label="clean seizure",
            true_color=:purple
        )
        save(joinpath(plot_dir, "whole_hvg_embedding_clean_seizure_truth.png"), clean_seizure_fig)
    end

    (
        embedding_df=embedding_df,
        feature_names=result.feature_names,
        distribution_df=result.distribution_df,
        skipped_df=result.skipped_df,
        output_dir=output_dir,
        plot_dir=plot_dir
    )
end

function whole_hvg_cluster_methods()
    if isdefined(Main, :WHOLE_HVG_CLUSTER_METHODS)
        methods = Symbol.(collect(WHOLE_HVG_CLUSTER_METHODS))
    elseif isdefined(Main, :WHOLE_HVG_CLUSTER_METHOD)
        methods = [Symbol(WHOLE_HVG_CLUSTER_METHOD)]
    else
        methods = [:pca, :tsne]
    end
    isempty(methods) && throw(ArgumentError("At least one whole-HVG clustering method is required."))
    unique(methods)
end

function whole_hvg_method_output_dir(base_dir, method, methods)
    length(methods) == 1 ? base_dir : joinpath(base_dir, String(method))
end

methods = whole_hvg_cluster_methods()
tsne_max_rows = isdefined(Main, :WHOLE_HVG_TSNE_MAX_ROWS) ? WHOLE_HVG_TSNE_MAX_ROWS : 5000
tsne_seed = isdefined(Main, :WHOLE_HVG_TSNE_SEED) ? WHOLE_HVG_TSNE_SEED : 1
tsne_iterations = isdefined(Main, :WHOLE_HVG_TSNE_ITERATIONS) ? WHOLE_HVG_TSNE_ITERATIONS : 1000
tsne_perplexity = isdefined(Main, :WHOLE_HVG_TSNE_PERPLEXITY) ? WHOLE_HVG_TSNE_PERPLEXITY : 30.0
input_csv = isdefined(Main, :WHOLE_HVG_CLUSTER_INPUT_CSV) ?
    WHOLE_HVG_CLUSTER_INPUT_CSV :
    newest_whole_hvg_bin_scores()
compute_if_missing = isdefined(Main, :WHOLE_HVG_CLUSTER_COMPUTE_IF_MISSING) ?
    WHOLE_HVG_CLUSTER_COMPUTE_IF_MISSING :
    true

bin_df = if !isnothing(input_csv)
    @info "Loading whole-HVG bin scores for clustering" input_csv
    CSV.read(input_csv, DataFrame)
elseif compute_if_missing
    patients = isdefined(Main, :WHOLE_HVG_PATIENTS) ?
        WHOLE_HVG_PATIENTS :
        artifact_labeled_patients()
    @info "No whole-HVG bin score CSV found; computing features first" patients
    run_whole_hvg_anomaly_survey(;
        patients=patients,
        bin_s=isdefined(Main, :WHOLE_HVG_BIN_S) ? WHOLE_HVG_BIN_S : 1,
        artifact_grades=isdefined(Main, :WHOLE_HVG_ARTIFACT_GRADES) ? WHOLE_HVG_ARTIFACT_GRADES : [1, 2],
        min_reviewers_per_seizure=isdefined(Main, :WHOLE_HVG_MIN_REVIEWERS_PER_SEIZURE) ?
            WHOLE_HVG_MIN_REVIEWERS_PER_SEIZURE :
            3,
        save_node_scores=isdefined(Main, :WHOLE_HVG_SAVE_NODE_SCORES) ?
            WHOLE_HVG_SAVE_NODE_SCORES :
            false,
        include_edge_span_features=isdefined(Main, :WHOLE_HVG_INCLUDE_EDGE_SPAN_FEATURES) ?
            WHOLE_HVG_INCLUDE_EDGE_SPAN_FEATURES :
            true,
        include_centrality_features=isdefined(Main, :WHOLE_HVG_INCLUDE_CENTRALITY_FEATURES) ?
            WHOLE_HVG_INCLUDE_CENTRALITY_FEATURES :
            true,
        include_within_bin_graph_features=isdefined(Main, :WHOLE_HVG_INCLUDE_WITHIN_BIN_GRAPH_FEATURES) ?
            WHOLE_HVG_INCLUDE_WITHIN_BIN_GRAPH_FEATURES :
            true,
        include_within_bin_topology_features=isdefined(Main, :WHOLE_HVG_INCLUDE_WITHIN_BIN_TOPOLOGY_FEATURES) ?
            WHOLE_HVG_INCLUDE_WITHIN_BIN_TOPOLOGY_FEATURES :
            true,
        centrality_landmark_count=isdefined(Main, :WHOLE_HVG_CENTRALITY_LANDMARK_COUNT) ?
            WHOLE_HVG_CENTRALITY_LANDMARK_COUNT :
            64,
        centrality_max_nodes=isdefined(Main, :WHOLE_HVG_CENTRALITY_MAX_NODES) ?
            WHOLE_HVG_CENTRALITY_MAX_NODES :
            WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES,
        save_outputs=isdefined(Main, :WHOLE_HVG_SAVE_OUTPUTS) ? WHOLE_HVG_SAVE_OUTPUTS : true
    ).bin_df
else
    throw(ArgumentError("No whole-HVG bin score CSV found and WHOLE_HVG_CLUSTER_COMPUTE_IF_MISSING=false."))
end

output_root = isdefined(Main, :WHOLE_HVG_CLUSTER_OUTPUT_ROOT) ?
    WHOLE_HVG_CLUSTER_OUTPUT_ROOT :
    datadir("exp_pro", "whole_hvg_clustering", string(Dates.now()))
plot_root = isdefined(Main, :WHOLE_HVG_CLUSTER_PLOT_ROOT) ?
    WHOLE_HVG_CLUSTER_PLOT_ROOT :
    plotsdir("whole_hvg_clustering_$(Dates.now())")

method_results = Dict{Symbol,Any}()
for method in methods
    method_output_dir = whole_hvg_method_output_dir(output_root, method, methods)
    method_plot_dir = whole_hvg_method_output_dir(plot_root, method, methods)
    @info "Saving whole-HVG clustering outputs" method output_dir=method_output_dir plot_dir=method_plot_dir
    method_results[method] = save_whole_hvg_clustering_outputs(bin_df;
        output_dir=method_output_dir,
        plot_dir=method_plot_dir,
        method=method,
        tsne_max_rows=tsne_max_rows,
        tsne_seed=tsne_seed,
        tsne_iterations=tsne_iterations,
        tsne_perplexity=tsne_perplexity
    )
end

whole_hvg_clustering_results = length(methods) == 1 ?
    method_results[first(methods)] :
    method_results

whole_hvg_clustering_results
