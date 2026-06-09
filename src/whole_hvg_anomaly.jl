using CSV
using DataFrames
using Dates
using JLD2
using Statistics

import Graphs
import VisibilityGraphs:
    edge_list,
    horizontal_visibility_graph,
    hvg_degrees,
    hvg_random_degree_probability,
    to_graphs

const WHOLE_HVG_BASE_FEATURE_NAMES = [
    :hvg_degree,
    :hvg_indegree,
    :hvg_outdegree,
    :hvg_degree_imbalance,
    :hvg_random_law_surprisal
]

const WHOLE_HVG_EDGE_SPAN_FEATURE_NAMES = [
    :hvg_edge_span_mean,
    :hvg_edge_span_max
]

const WHOLE_HVG_CENTRALITY_FEATURE_NAMES = [
    :hvg_pagerank,
    :hvg_eigenvector_centrality,
    :hvg_sampled_betweenness,
    :hvg_landmark_closeness
]

const WHOLE_HVG_WITHIN_BIN_GRAPH_FEATURE_NAMES = [
    :hvg_bin_vertex_count,
    :hvg_bin_edge_count,
    :hvg_bin_edge_density,
    :hvg_bin_mean_degree,
    :hvg_bin_degree_std,
    :hvg_bin_degree_max,
    :hvg_bin_degree_q95,
    :hvg_bin_degree_entropy,
    :hvg_bin_random_law_surprisal_mean,
    :hvg_bin_edge_span_mean,
    :hvg_bin_edge_span_max,
    :hvg_bin_edge_span_q95
]

const WHOLE_HVG_WITHIN_BIN_TOPOLOGY_FEATURE_NAMES = [
    :hvg_bin_global_clustering_coefficient,
    :hvg_bin_degree_assortativity
]

const WHOLE_HVG_AGGREGATORS = (:mean, :std, :max, :q95)
const WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES = 20_000

function whole_hvg_feature_names(;
        include_edge_span_features=false,
        include_centrality_features=false
    )
    feature_names = Symbol[WHOLE_HVG_BASE_FEATURE_NAMES...]
    include_edge_span_features && append!(feature_names, WHOLE_HVG_EDGE_SPAN_FEATURE_NAMES)
    include_centrality_features && append!(feature_names, WHOLE_HVG_CENTRALITY_FEATURE_NAMES)
    feature_names
end

function whole_hvg_metric_notes(;
        include_edge_span_features=false,
        include_centrality_features=false,
        include_within_bin_graph_features=false,
        include_within_bin_topology_features=true
    )
    within_bin_graph_status = include_within_bin_graph_features ?
        "included_optional" :
        "optional_skipped"
    within_bin_topology_status =
        include_within_bin_graph_features && include_within_bin_topology_features ?
            "included_optional" :
            "optional_skipped"

    rows = [
        (
            metric="hvg_degree",
            status="included",
            reason="Computed with VisibilityGraphs.hvg_degrees over the full channel recording."
        ),
        (
            metric="hvg_indegree",
            status="included",
            reason="Computed with VisibilityGraphs.hvg_degrees direction=:in over the full channel recording."
        ),
        (
            metric="hvg_outdegree",
            status="included",
            reason="Computed with VisibilityGraphs.hvg_degrees direction=:out over the full channel recording."
        ),
        (
            metric="hvg_degree_imbalance",
            status="included",
            reason="Outdegree minus indegree; captures local time-asymmetric visibility structure."
        ),
        (
            metric="hvg_random_law_surprisal",
            status="included",
            reason="Negative log probability under the analytic i.i.d. HVG degree law."
        ),
        (
            metric="hvg_edge_span",
            status=include_edge_span_features ? "included_optional" : "optional_skipped",
            reason="Requires retaining full HVG edge lists per channel; disabled by default to limit memory."
        ),
        (
            metric="hvg_pagerank",
            status=include_centrality_features ? "included_optional" : "optional_skipped",
            reason="Undirected PageRank centrality over the full recording; long recordings use deterministic node sampling to bound runtime."
        ),
        (
            metric="hvg_eigenvector_centrality",
            status=include_centrality_features ? "included_optional" : "optional_skipped",
            reason="Undirected eigenvector centrality over the full recording; long recordings use deterministic node sampling to bound runtime."
        ),
        (
            metric="hvg_sampled_betweenness",
            status=include_centrality_features ? "included_optional" : "optional_skipped",
            reason="Betweenness centrality sampled from deterministic evenly spaced landmark nodes; long recordings use deterministic node sampling to bound runtime."
        ),
        (
            metric="hvg_landmark_closeness",
            status=include_centrality_features ? "included_optional" : "optional_skipped",
            reason="Approximate closeness to deterministic evenly spaced landmark nodes; long recordings use deterministic node sampling to bound runtime."
        ),
        (
            metric="hvg_bin_graph_metrics",
            status=within_bin_graph_status,
            reason="Graph-wide size, density, degree, random-law, and edge-span metrics computed on each bin's own HVG."
        ),
        (
            metric="hvg_bin_graph_topology_metrics",
            status=within_bin_topology_status,
            reason="Global clustering coefficient and degree assortativity computed on each bin's own undirected HVG."
        ),
        (
            metric="connectedness_components",
            status="skipped",
            reason="Valid HVGs are structurally connected through adjacent-sample visibility edges."
        ),
        (
            metric="betweenness_closeness_diameter_radius_eccentricity",
            status="skipped",
            reason="Exact shortest-path metrics are too expensive for sample-level full-recording graphs."
        ),
        (
            metric="community_clique_metrics",
            status="skipped",
            reason="Expensive and semantically unstable on long per-channel time-series visibility graphs."
        ),
        (
            metric="nvg_metrics",
            status="deferred",
            reason="Natural visibility graphs are useful for comparison but heavier than HVGs for whole recordings."
        )
    ]
    DataFrame(rows)
end

function whole_hvg_signal_values(signal::AbstractVector)
    values = Vector{Float64}(undef, length(signal))
    for idx in eachindex(signal)
        value = signal[idx]
        if ismissing(value)
            return nothing
        end
        value = Float64(value)
        if !isfinite(value)
            return nothing
        end
        values[idx] = value
    end
    values
end

function whole_hvg_random_law_surprisal(degrees)
    [-(log(max(hvg_random_degree_probability(Int(degree)), eps(Float64)))) for degree in degrees]
end

function whole_hvg_nan_features(feature_names, n_nodes)
    Dict{Symbol,Vector{Float64}}(
        feature_name => fill(NaN, n_nodes) for feature_name in feature_names
    )
end

function whole_hvg_robust_zscores(values::AbstractVector)
    zscores = fill(NaN, length(values))
    finite_values = Float64[Float64(value) for value in values if isfinite(value)]
    isempty(finite_values) && return zscores

    center, scale = robust_center_and_scale(finite_values)
    if !isfinite(center)
        center = 0.0
    end
    if !isfinite(scale) || scale <= eps(Float64)
        scale = 1.0
    end

    for idx in eachindex(values)
        value = Float64(values[idx])
        zscores[idx] = isfinite(value) ? (value - center) / scale : NaN
    end
    zscores
end

function whole_hvg_edge_span_features(graph)
    n_nodes = length(graph)
    span_sum = zeros(Float64, n_nodes)
    span_count = zeros(Int, n_nodes)
    span_max = zeros(Float64, n_nodes)

    for (src, dst) in edge_list(graph)
        span = Float64(dst - src)
        span_sum[src] += span
        span_sum[dst] += span
        span_count[src] += 1
        span_count[dst] += 1
        span_max[src] = max(span_max[src], span)
        span_max[dst] = max(span_max[dst], span)
    end

    span_mean = Vector{Float64}(undef, n_nodes)
    for idx in 1:n_nodes
        span_mean[idx] = span_count[idx] == 0 ? NaN : span_sum[idx] / span_count[idx]
        if span_count[idx] == 0
            span_max[idx] = NaN
        end
    end

    Dict{Symbol,Vector{Float64}}(
        :hvg_edge_span_mean => span_mean,
        :hvg_edge_span_max => span_max
    )
end

function whole_hvg_landmark_indices(n_nodes::Integer, landmark_count::Integer)
    n_nodes >= 0 || throw(ArgumentError("n_nodes must be nonnegative"))
    landmark_count > 0 || throw(ArgumentError("centrality_landmark_count must be positive"))
    n_nodes == 0 && return Int[]

    count = min(n_nodes, landmark_count)
    count == n_nodes && return collect(1:n_nodes)
    unique(clamp.(round.(Int, range(1, n_nodes; length=count)), 1, n_nodes))
end

function whole_hvg_centrality_sample_indices(n_nodes::Integer, centrality_max_nodes::Integer)
    centrality_max_nodes > 1 ||
        throw(ArgumentError("centrality_max_nodes must be greater than 1"))
    n_nodes <= centrality_max_nodes && return collect(1:n_nodes)
    whole_hvg_landmark_indices(n_nodes, centrality_max_nodes)
end

function whole_hvg_expand_sampled_values(sample_indices, sample_values, n_nodes::Integer)
    if length(sample_indices) == n_nodes
        return Float64.(sample_values)
    end

    expanded = Vector{Float64}(undef, n_nodes)
    nearest_sample_idx = 1
    for node_idx in 1:n_nodes
        while nearest_sample_idx < length(sample_indices) &&
                abs(sample_indices[nearest_sample_idx + 1] - node_idx) <=
                abs(sample_indices[nearest_sample_idx] - node_idx)
            nearest_sample_idx += 1
        end
        expanded[node_idx] = sample_values[nearest_sample_idx]
    end
    expanded
end

function whole_hvg_expand_centrality_features(raw_features, sample_indices, n_nodes::Integer)
    Dict{Symbol,Vector{Float64}}(
        feature_name => whole_hvg_expand_sampled_values(sample_indices, values, n_nodes)
        for (feature_name, values) in raw_features
    )
end

function whole_hvg_landmark_closeness(graph, landmarks::AbstractVector{<:Integer})
    n_nodes = Graphs.nv(graph)
    isempty(landmarks) && return fill(NaN, n_nodes)

    distance_sums = zeros(Float64, n_nodes)
    reachable_counts = zeros(Int, n_nodes)
    landmark_flags = falses(n_nodes)
    for landmark in landmarks
        landmark_flags[landmark] = true
        distances = Graphs.gdistances(graph, landmark)
        unreachable = typemax(eltype(distances))
        for node_idx in 1:n_nodes
            distance = distances[node_idx]
            if distance != unreachable && distance > 0
                distance_sums[node_idx] += Float64(distance)
                reachable_counts[node_idx] += 1
            end
        end
    end

    scores = zeros(Float64, n_nodes)
    for node_idx in 1:n_nodes
        distance_sums[node_idx] > 0 || continue
        possible_count = length(landmarks) - (landmark_flags[node_idx] ? 1 : 0)
        possible_count = max(possible_count, 1)
        scores[node_idx] = reachable_counts[node_idx] / distance_sums[node_idx]
        scores[node_idx] *= reachable_counts[node_idx] / possible_count
    end
    scores
end

function whole_hvg_sampled_betweenness(graph, landmarks::AbstractVector{<:Integer})
    n_nodes = Graphs.nv(graph)
    scores = zeros(Float64, n_nodes)
    (n_nodes <= 2 || isempty(landmarks)) && return scores

    distances = fill(-1, n_nodes)
    path_counts = zeros(Float64, n_nodes)
    dependencies = zeros(Float64, n_nodes)
    stack = Vector{Int}(undef, n_nodes)
    queue = Vector{Int}(undef, n_nodes)

    for source in landmarks
        fill!(distances, -1)
        fill!(path_counts, 0.0)
        fill!(dependencies, 0.0)

        queue_head = 1
        queue_tail = 1
        stack_len = 0
        queue[queue_tail] = source
        distances[source] = 0
        path_counts[source] = 1.0

        while queue_head <= queue_tail
            vertex = queue[queue_head]
            queue_head += 1
            stack_len += 1
            stack[stack_len] = vertex

            for neighbor in Graphs.neighbors(graph, vertex)
                if distances[neighbor] < 0
                    queue_tail += 1
                    queue[queue_tail] = neighbor
                    distances[neighbor] = distances[vertex] + 1
                end
                if distances[neighbor] == distances[vertex] + 1
                    path_counts[neighbor] += path_counts[vertex]
                end
            end
        end

        for stack_idx in stack_len:-1:1
            vertex = stack[stack_idx]
            path_counts[vertex] == 0 && continue
            coefficient = (1.0 + dependencies[vertex]) / path_counts[vertex]
            for predecessor in Graphs.neighbors(graph, vertex)
                if distances[predecessor] == distances[vertex] - 1
                    dependencies[predecessor] += path_counts[predecessor] * coefficient
                end
            end
            vertex == source || (scores[vertex] += dependencies[vertex])
        end
    end

    scores .*= (1.0 / ((n_nodes - 1) * (n_nodes - 2))) * n_nodes / length(landmarks)
    scores
end

function whole_hvg_safe_centrality_values(compute::Function, feature_name::Symbol, n_nodes::Integer)
    try
        values = Float64.(compute())
        if length(values) != n_nodes
            throw(ArgumentError("expected $(n_nodes) values, got $(length(values))"))
        end
        values
    catch err
        @warn "Failed to compute whole-HVG centrality feature" feature=feature_name exception=(err, catch_backtrace())
        fill(NaN, n_nodes)
    end
end

function whole_hvg_centrality_features(
        values,
        graph;
        tie_policy=:strict,
        centrality_landmark_count=64,
        centrality_max_nodes=WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES
    )
    n_full_nodes = length(values)
    sample_indices = whole_hvg_centrality_sample_indices(
        n_full_nodes,
        Int(centrality_max_nodes)
    )
    centrality_graph = length(sample_indices) == n_full_nodes ?
        graph :
        horizontal_visibility_graph(values[sample_indices]; tie_policy=tie_policy)

    hvg_graph = to_graphs(centrality_graph; directed=false)
    n_nodes = Graphs.nv(hvg_graph)
    landmarks = whole_hvg_landmark_indices(n_nodes, Int(centrality_landmark_count))

    raw_features = Dict{Symbol,Vector{Float64}}(
        :hvg_pagerank => whole_hvg_safe_centrality_values(:hvg_pagerank, n_nodes) do
            Graphs.pagerank(hvg_graph, 0.85, 100, 1.0e-6)
        end,
        :hvg_eigenvector_centrality => whole_hvg_safe_centrality_values(:hvg_eigenvector_centrality, n_nodes) do
            Graphs.eigenvector_centrality(hvg_graph)
        end,
        :hvg_sampled_betweenness => whole_hvg_safe_centrality_values(:hvg_sampled_betweenness, n_nodes) do
            whole_hvg_sampled_betweenness(hvg_graph, landmarks)
        end,
        :hvg_landmark_closeness => whole_hvg_safe_centrality_values(:hvg_landmark_closeness, n_nodes) do
            whole_hvg_landmark_closeness(hvg_graph, landmarks)
        end
    )
    whole_hvg_expand_centrality_features(raw_features, sample_indices, n_full_nodes)
end

function whole_hvg_node_features(
        signal::AbstractVector;
        tie_policy=:strict,
        include_edge_span_features=false,
        include_centrality_features=false,
        centrality_landmark_count=64,
        centrality_max_nodes=WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES
    )
    values = whole_hvg_signal_values(signal)
    n_nodes = length(signal)
    feature_values = Dict{Symbol,Vector{Float64}}()

    if isnothing(values) || n_nodes < 2
        merge!(
            feature_values,
            whole_hvg_nan_features(
                whole_hvg_feature_names(;
                    include_edge_span_features=include_edge_span_features,
                    include_centrality_features=include_centrality_features
                ),
                n_nodes
            )
        )
    else
        degree_values = Float64.(hvg_degrees(values; direction=:both, tie_policy=tie_policy))
        indegree_values = Float64.(hvg_degrees(values; direction=:in, tie_policy=tie_policy))
        outdegree_values = Float64.(hvg_degrees(values; direction=:out, tie_policy=tie_policy))

        feature_values[:hvg_degree] = degree_values
        feature_values[:hvg_indegree] = indegree_values
        feature_values[:hvg_outdegree] = outdegree_values
        feature_values[:hvg_degree_imbalance] = outdegree_values .- indegree_values
        feature_values[:hvg_random_law_surprisal] = whole_hvg_random_law_surprisal(degree_values)

        if include_edge_span_features || include_centrality_features
            graph = horizontal_visibility_graph(values; tie_policy=tie_policy)
            if include_edge_span_features
                merge!(feature_values, whole_hvg_edge_span_features(graph))
            end
            if include_centrality_features
                merge!(
                    feature_values,
                    whole_hvg_centrality_features(values, graph;
                        tie_policy=tie_policy,
                        centrality_landmark_count=centrality_landmark_count,
                        centrality_max_nodes=centrality_max_nodes
                    )
                )
            end
        end
    end

    robust_zscores = Dict{Symbol,Vector{Float64}}()
    abs_robust_zscores = Dict{Symbol,Vector{Float64}}()
    for (feature_name, values) in feature_values
        zscores = whole_hvg_robust_zscores(values)
        robust_zscores[feature_name] = zscores
        abs_robust_zscores[feature_name] = abs.(zscores)
    end

    node_anomaly = whole_hvg_node_anomaly(abs_robust_zscores, n_nodes)
    (
        feature_values=feature_values,
        robust_zscores=robust_zscores,
        abs_robust_zscores=abs_robust_zscores,
        node_anomaly=node_anomaly
    )
end

function whole_hvg_node_anomaly(abs_robust_zscores, n_nodes)
    node_anomaly = fill(NaN, n_nodes)
    for node_idx in 1:n_nodes
        score = -Inf
        for feature_values in Base.values(abs_robust_zscores)
            value = feature_values[node_idx]
            if isfinite(value)
                score = max(score, value)
            end
        end
        node_anomaly[node_idx] = isfinite(score) ? score : NaN
    end
    node_anomaly
end

function whole_hvg_feature_values_for_range(values::AbstractVector, sample_range)
    isempty(sample_range) && return Float64[]
    start_idx = max(first(sample_range), firstindex(values))
    stop_idx = min(last(sample_range), lastindex(values))
    stop_idx < start_idx && return Float64[]
    Float64[Float64(value) for value in values[start_idx:stop_idx] if isfinite(value)]
end

function whole_hvg_sample_count_for_range(values::AbstractVector, sample_range)
    isempty(sample_range) && return 0
    start_idx = max(first(sample_range), firstindex(values))
    stop_idx = min(last(sample_range), lastindex(values))
    max(0, stop_idx - start_idx + 1)
end

function whole_hvg_sample_index_bounds(values::AbstractVector, sample_range)
    isempty(sample_range) && return (0, -1)
    start_idx = max(first(sample_range), firstindex(values))
    stop_idx = min(last(sample_range), lastindex(values))
    stop_idx < start_idx ? (0, -1) : (start_idx, stop_idx)
end

function whole_hvg_bin_signal_values(signal::AbstractVector, sample_range)
    isempty(sample_range) && return Float64[]
    start_idx = max(first(sample_range), firstindex(signal))
    stop_idx = min(last(sample_range), lastindex(signal))
    stop_idx < start_idx && return Float64[]

    values = Vector{Float64}(undef, stop_idx - start_idx + 1)
    out_idx = 1
    for idx in start_idx:stop_idx
        value = signal[idx]
        if ismissing(value)
            return nothing
        end
        value = Float64(value)
        if !isfinite(value)
            return nothing
        end
        values[out_idx] = value
        out_idx += 1
    end
    values
end

function whole_hvg_within_bin_graph_feature_names(; include_within_bin_topology_features=true)
    feature_names = Symbol[WHOLE_HVG_WITHIN_BIN_GRAPH_FEATURE_NAMES...]
    if include_within_bin_topology_features
        append!(feature_names, WHOLE_HVG_WITHIN_BIN_TOPOLOGY_FEATURE_NAMES)
    end
    feature_names
end

function whole_hvg_nan_within_bin_graph_pairs(; include_within_bin_topology_features=true)
    Pair{Symbol,Any}[
        feature_name => NaN for feature_name in whole_hvg_within_bin_graph_feature_names(;
            include_within_bin_topology_features=include_within_bin_topology_features
        )
    ]
end

function whole_hvg_degree_entropy(degrees::AbstractVector{<:Integer})
    n_nodes = length(degrees)
    n_nodes == 0 && return NaN

    counts = Dict{Int,Int}()
    for degree in degrees
        counts[degree] = get(counts, degree, 0) + 1
    end

    entropy = 0.0
    for count in Base.values(counts)
        probability = count / n_nodes
        entropy -= probability * log(probability)
    end
    entropy
end

function whole_hvg_safe_within_bin_graph_value(compute::Function, feature_name::Symbol)
    try
        value = Float64(compute())
        isfinite(value) ? value : NaN
    catch err
        @warn "Failed to compute whole-HVG within-bin graph feature" feature=feature_name exception=(err, catch_backtrace())
        NaN
    end
end

function whole_hvg_within_bin_graph_feature_pairs(
        signal::AbstractVector,
        sample_range;
        tie_policy=:strict,
        include_within_bin_topology_features=true
    )
    values = whole_hvg_bin_signal_values(signal, sample_range)
    if isnothing(values)
        return whole_hvg_nan_within_bin_graph_pairs(;
            include_within_bin_topology_features=include_within_bin_topology_features
        )
    end

    n_nodes = length(values)
    graph = nothing
    if n_nodes >= 2
        try
            graph = horizontal_visibility_graph(values; tie_policy=tie_policy)
        catch err
            @warn "Failed to construct whole-HVG within-bin graph" exception=(err, catch_backtrace())
            return whole_hvg_nan_within_bin_graph_pairs(;
                include_within_bin_topology_features=include_within_bin_topology_features
            )
        end
    end

    degree_counts = zeros(Int, n_nodes)
    edge_spans = Float64[]
    if !isnothing(graph)
        for (src, dst) in edge_list(graph)
            degree_counts[src] += 1
            degree_counts[dst] += 1
            push!(edge_spans, Float64(abs(dst - src)))
        end
    end

    degree_values = Float64.(degree_counts)
    edge_count = length(edge_spans)
    edge_density = n_nodes >= 2 ? (2.0 * edge_count) / (n_nodes * (n_nodes - 1)) : NaN

    pairs = Pair{Symbol,Any}[
        :hvg_bin_vertex_count => Float64(n_nodes),
        :hvg_bin_edge_count => Float64(edge_count),
        :hvg_bin_edge_density => edge_density,
        :hvg_bin_mean_degree => n_nodes > 0 ? mean(degree_values) : NaN,
        :hvg_bin_degree_std => n_nodes > 1 ? std(degree_values) : (n_nodes == 1 ? 0.0 : NaN),
        :hvg_bin_degree_max => n_nodes > 0 ? maximum(degree_values) : NaN,
        :hvg_bin_degree_q95 => n_nodes > 0 ? quantile(degree_values, 0.95) : NaN,
        :hvg_bin_degree_entropy => whole_hvg_degree_entropy(degree_counts),
        :hvg_bin_random_law_surprisal_mean =>
            n_nodes > 0 ? mean(whole_hvg_random_law_surprisal(degree_values)) : NaN,
        :hvg_bin_edge_span_mean => edge_count > 0 ? mean(edge_spans) : NaN,
        :hvg_bin_edge_span_max => edge_count > 0 ? maximum(edge_spans) : NaN,
        :hvg_bin_edge_span_q95 => edge_count > 0 ? quantile(edge_spans, 0.95) : NaN
    ]

    if include_within_bin_topology_features
        hvg_graph = if isnothing(graph)
            Graphs.SimpleGraph(n_nodes)
        else
            to_graphs(graph; directed=false)
        end
        append!(
            pairs,
            Pair{Symbol,Any}[
                :hvg_bin_global_clustering_coefficient => n_nodes >= 3 ?
                    whole_hvg_safe_within_bin_graph_value(
                        :hvg_bin_global_clustering_coefficient
                    ) do
                        Graphs.global_clustering_coefficient(hvg_graph)
                    end :
                    NaN,
                :hvg_bin_degree_assortativity => n_nodes >= 2 ?
                    whole_hvg_safe_within_bin_graph_value(
                        :hvg_bin_degree_assortativity
                    ) do
                        Graphs.assortativity(hvg_graph)
                    end :
                    NaN
            ]
        )
    end

    pairs
end

function whole_hvg_aggregate_pairs(prefix::Symbol, values)
    finite_values = Float64[Float64(value) for value in values if isfinite(value)]
    if isempty(finite_values)
        return Pair{Symbol,Any}[
            Symbol(prefix, :_mean) => NaN,
            Symbol(prefix, :_std) => NaN,
            Symbol(prefix, :_max) => NaN,
            Symbol(prefix, :_q95) => NaN
        ]
    end

    Pair{Symbol,Any}[
        Symbol(prefix, :_mean) => mean(finite_values),
        Symbol(prefix, :_std) => length(finite_values) > 1 ? std(finite_values) : 0.0,
        Symbol(prefix, :_max) => maximum(finite_values),
        Symbol(prefix, :_q95) => quantile(finite_values, 0.95)
    ]
end

function whole_hvg_channel_label(eeg, channel)
    if hasproperty(eeg, :labels)
        labels = getproperty(eeg, :labels)
        if channel in eachindex(labels)
            return String(labels[channel])
        end
    end
    "channel_$(channel)"
end

function whole_hvg_recording_duration(eeg, recording_duration)
    !isnothing(recording_duration) && return Float64(recording_duration)
    hasproperty(eeg, :duration) ||
        throw(ArgumentError("recording_duration must be provided for EEG values without a duration field"))
    Float64(getproperty(eeg, :duration))
end

function whole_hvg_default_seizure_bounds(eeg)
    hasproperty(eeg, :seizure_annotations) ? getproperty(eeg, :seizure_annotations) : Tuple{Float64,Float64}[]
end

function whole_hvg_bin_feature_pairs(node_features, sample_range)
    pairs = Pair{Symbol,Any}[]

    for feature_name in sort(collect(keys(node_features.feature_values)); by=String)
        raw_values = whole_hvg_feature_values_for_range(
            node_features.feature_values[feature_name],
            sample_range
        )
        append!(pairs, whole_hvg_aggregate_pairs(feature_name, raw_values))

        abs_z_values = whole_hvg_feature_values_for_range(
            node_features.abs_robust_zscores[feature_name],
            sample_range
        )
        append!(pairs, whole_hvg_aggregate_pairs(Symbol(feature_name, :_abs_robust_z), abs_z_values))
    end

    node_anomaly_values = whole_hvg_feature_values_for_range(node_features.node_anomaly, sample_range)
    append!(pairs, whole_hvg_aggregate_pairs(:hvg_node_anomaly, node_anomaly_values))
    pairs
end

function whole_hvg_node_payload(patient, channel, channel_label, node_features)
    Dict{String,Any}(
        "patient" => patient,
        "channel" => channel,
        "channel_label" => channel_label,
        "feature_values" => node_features.feature_values,
        "robust_zscores" => node_features.robust_zscores,
        "abs_robust_zscores" => node_features.abs_robust_zscores,
        "node_anomaly" => node_features.node_anomaly
    )
end

function whole_hvg_bin_score_rows(
        eeg::AbstractProcessedEEG;
        patient,
        artifact_bounds=Tuple{Float64,Float64}[],
        seizure_bounds=whole_hvg_default_seizure_bounds(eeg),
        bin_s=15,
        recording_duration=nothing,
        tie_policy=:strict,
        include_edge_span_features=false,
        include_centrality_features=false,
        include_within_bin_graph_features=false,
        include_within_bin_topology_features=true,
        centrality_landmark_count=64,
        centrality_max_nodes=WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES,
        save_node_scores=false
    )
    recording_duration = whole_hvg_recording_duration(eeg, recording_duration)
    bin_starts, bin_stops = artifact_bin_bounds(recording_duration, bin_s)
    sample_ranges = artifact_sample_ranges(eeg, bin_starts, bin_stops)
    artifact_truth = bounds_to_bin_labels(artifact_bounds, bin_starts, bin_stops)
    seizure_truth = bounds_to_bin_labels(seizure_bounds, bin_starts, bin_stops)
    seizure_event_ids = seizure_event_ids_for_bins(seizure_bounds, bin_starts, bin_stops)

    rows = NamedTuple[]
    node_scores = Dict{String,Any}()
    signals = getproperty(eeg, :signals)

    for channel in axes(signals, 1)
        channel_label = whole_hvg_channel_label(eeg, channel)
        signal = vec(signals[channel, :])
        node_features = whole_hvg_node_features(signal;
            tie_policy=tie_policy,
            include_edge_span_features=include_edge_span_features,
            include_centrality_features=include_centrality_features,
            centrality_landmark_count=centrality_landmark_count,
            centrality_max_nodes=centrality_max_nodes
        )

        if save_node_scores
            key = "patient_$(patient)_channel_$(channel)"
            node_scores[key] = whole_hvg_node_payload(patient, channel, channel_label, node_features)
        end

        for bin_idx in eachindex(sample_ranges)
            sample_start_idx, sample_stop_idx =
                whole_hvg_sample_index_bounds(node_features.node_anomaly, sample_ranges[bin_idx])
            pairs = Pair{Symbol,Any}[
                :patient => patient,
                :channel => channel,
                :channel_label => channel_label,
                :bin_index => bin_idx,
                :bin_start => bin_starts[bin_idx],
                :bin_stop => bin_stops[bin_idx],
                :sample_start_index => sample_start_idx,
                :sample_stop_index => sample_stop_idx,
                :sample_count => whole_hvg_sample_count_for_range(
                    node_features.node_anomaly,
                    sample_ranges[bin_idx]
                ),
                :artifact_truth => artifact_truth[bin_idx],
                :seizure_truth => seizure_truth[bin_idx],
                :clean_seizure_truth => seizure_truth[bin_idx] && !artifact_truth[bin_idx],
                :seizure_event_id => seizure_event_ids[bin_idx]
            ]
            append!(pairs, whole_hvg_bin_feature_pairs(node_features, sample_ranges[bin_idx]))
            if include_within_bin_graph_features
                append!(
                    pairs,
                    whole_hvg_within_bin_graph_feature_pairs(
                        signal,
                        sample_ranges[bin_idx];
                        tie_policy=tie_policy,
                        include_within_bin_topology_features=include_within_bin_topology_features
                    )
                )
            end
            push!(rows, (; pairs...))
        end
    end

    (
        bin_df=DataFrame(rows),
        node_scores=node_scores
    )
end

function whole_hvg_bin_score_rows_for_patient(
        patient_num;
        bin_s=15,
        artifact_grades=[1, 2],
        min_reviewers_per_seizure=3,
        artifact_csv_path=scriptsdir("helsinki_artifacts.csv"),
        tie_policy=:strict,
        include_edge_span_features=false,
        include_centrality_features=false,
        include_within_bin_graph_features=false,
        include_within_bin_topology_features=true,
        centrality_landmark_count=64,
        centrality_max_nodes=WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES,
        save_node_scores=false
    )
    eeg = load_helsinki_eeg(patient_num;
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        excluded_artifact_grades=Int[]
    )
    artifact_bounds = load_helsinki_artifact_ground_truth(patient_num, artifact_grades;
        csv_path=artifact_csv_path,
        recording_duration=eeg.duration
    )

    whole_hvg_bin_score_rows(eeg;
        patient=patient_num,
        artifact_bounds=artifact_bounds,
        seizure_bounds=eeg.seizure_annotations,
        bin_s=bin_s,
        tie_policy=tie_policy,
        include_edge_span_features=include_edge_span_features,
        include_centrality_features=include_centrality_features,
        include_within_bin_graph_features=include_within_bin_graph_features,
        include_within_bin_topology_features=include_within_bin_topology_features,
        centrality_landmark_count=centrality_landmark_count,
        centrality_max_nodes=centrality_max_nodes,
        save_node_scores=save_node_scores
    )
end

function run_whole_hvg_anomaly_survey(;
        patients=artifact_labeled_patients(),
        bin_s=15,
        artifact_grades=[1, 2],
        min_reviewers_per_seizure=3,
        artifact_csv_path=scriptsdir("helsinki_artifacts.csv"),
        tie_policy=:strict,
        include_edge_span_features=false,
        include_centrality_features=false,
        include_within_bin_graph_features=false,
        include_within_bin_topology_features=true,
        centrality_landmark_count=64,
        centrality_max_nodes=WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES,
        save_node_scores=false,
        output_root=datadir("exp_pro", "whole_hvg_anomaly"),
        save_outputs=true
    )
    bin_dfs = DataFrame[]
    node_scores = Dict{String,Any}()
    event_offset = 0

    for patient in patients
        @info "Scoring whole-recording HVG anomaly features" patient
        patient_results = whole_hvg_bin_score_rows_for_patient(patient;
            bin_s=bin_s,
            artifact_grades=artifact_grades,
            min_reviewers_per_seizure=min_reviewers_per_seizure,
            artifact_csv_path=artifact_csv_path,
            tie_policy=tie_policy,
            include_edge_span_features=include_edge_span_features,
            include_centrality_features=include_centrality_features,
            include_within_bin_graph_features=include_within_bin_graph_features,
            include_within_bin_topology_features=include_within_bin_topology_features,
            centrality_landmark_count=centrality_landmark_count,
            centrality_max_nodes=centrality_max_nodes,
            save_node_scores=save_node_scores
        )

        patient_df = patient_results.bin_df
        local_event_count = nrow(patient_df) == 0 ? 0 : maximum(patient_df.seizure_event_id)
        patient_df.seizure_event_id = offset_event_ids(patient_df, event_offset)
        event_offset += local_event_count
        push!(bin_dfs, patient_df)

        if save_node_scores
            merge!(node_scores, patient_results.node_scores)
        end
    end

    bin_df = isempty(bin_dfs) ? DataFrame() : reduce(vcat, bin_dfs)
    metric_notes_df = whole_hvg_metric_notes(
        include_edge_span_features=include_edge_span_features,
        include_centrality_features=include_centrality_features,
        include_within_bin_graph_features=include_within_bin_graph_features,
        include_within_bin_topology_features=include_within_bin_topology_features
    )

    if save_outputs
        session_id = Dates.now()
        output_dir = joinpath(output_root, string(session_id))
        mkpath(output_dir)
        CSV.write(joinpath(output_dir, "whole_hvg_bin_scores.csv"), bin_df)
        CSV.write(joinpath(output_dir, "whole_hvg_metric_notes.csv"), metric_notes_df)
        if save_node_scores
            jldsave(joinpath(output_dir, "whole_hvg_node_scores.jld2"); node_scores=node_scores)
        end
        @info "Saved whole-recording HVG anomaly survey" output_dir
        return (
            bin_df=bin_df,
            metric_notes_df=metric_notes_df,
            node_scores=save_node_scores ? node_scores : nothing,
            output_dir=output_dir
        )
    else
        return (
            bin_df=bin_df,
            metric_notes_df=metric_notes_df,
            node_scores=save_node_scores ? node_scores : nothing,
            output_dir=nothing
        )
    end
end
