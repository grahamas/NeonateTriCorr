using Test

using DrWatson
@quickactivate "NeonateTriCorr"

include(srcdir("types.jl"))
include(srcdir("detect_artifacts.jl"))
include(srcdir("whole_hvg_anomaly.jl"))

struct ArtifactScoreTestEEG <: AbstractProcessedEEG
    signals::Matrix{Float64}
    sample_rate::Int
end

struct WholeHVGTestEEG <: AbstractProcessedEEG
    signals::Matrix{Float64}
    sample_rate::Int
    duration::Float64
    labels::Vector{String}
    seizure_annotations::Vector{Tuple{Float64,Float64}}
end

function repeated_window_eeg(; sample_rate=16)
    samples_per_bin = 32
    base = [sin(2 * pi * i / samples_per_bin) for i in 1:samples_per_bin]
    artifact = copy(base)
    artifact[div(samples_per_bin, 2) + 1] = 12.0
    signal = vcat(base, base, base, artifact)
    ArtifactScoreTestEEG(reshape(signal, 1, :), sample_rate)
end

function synthetic_artifact_score_df()
    rows = NamedTuple[]
    for patient in 1:3
        for bin_idx in 1:6
            artifact_truth = bin_idx >= 4
            for method in ["method_a", "method_b"]
                score = if method == "method_a"
                    artifact_truth ? 4.0 + 0.1 * bin_idx : -4.0 + 0.1 * bin_idx
                else
                    artifact_truth ? 2.0 + 0.2 * patient : -2.0 - 0.1 * patient
                end
                push!(rows, (
                    patient=patient,
                    method=method,
                    method_description="synthetic",
                    bin_index=bin_idx,
                    bin_start=15.0 * (bin_idx - 1),
                    bin_stop=15.0 * bin_idx,
                    artifact_truth=artifact_truth,
                    seizure_truth=false,
                    clean_seizure_truth=false,
                    seizure_event_id=0,
                    score=score
                ))
            end
        end
    end
    DataFrame(rows)
end

@testset "Visibility graph artifact scores" begin
    eeg = repeated_window_eeg()
    sample_ranges = [1:32, 33:64, 65:96, 97:128]

    hvg_scores = score_hvg_degree_anomaly(eeg, sample_ranges; sample_stride=1)
    nvg_scores = score_nvg_degree_anomaly(eeg, sample_ranges; sample_stride=1)
    nvg_stride8_scores = score_nvg_degree_anomaly(eeg, sample_ranges; sample_stride=8)

    @test length(hvg_scores) == length(sample_ranges)
    @test length(nvg_scores) == length(sample_ranges)
    @test length(nvg_stride8_scores) == length(sample_ranges)
    @test all(isfinite, hvg_scores)
    @test all(isfinite, nvg_scores)
    @test all(isfinite, nvg_stride8_scores)
    @test hvg_scores[end] > maximum(hvg_scores[begin:end-1])
    @test nvg_scores[end] > maximum(nvg_scores[begin:end-1])

    short_eeg = ArtifactScoreTestEEG(reshape([1.0], 1, :), 1)
    @test isnan(only(score_hvg_degree_anomaly(short_eeg, [1:1]; sample_stride=1)))
    @test_throws ArgumentError score_hvg_degree_anomaly(eeg, sample_ranges; sample_stride=0)
end

@testset "Visibility graph artifact methods" begin
    names = getfield.(artifact_detection_methods(visibility_graph_sample_stride=2), :name)
    @test "hvg_degree_anomaly_stride2" in names
    @test "nvg_degree_anomaly_stride2" in names

    decoupled_names = getfield.(artifact_detection_methods(
        hvg_sample_strides=(1, 4),
        nvg_sample_strides=(8,)
    ), :name)
    @test "hvg_degree_anomaly_stride1" in decoupled_names
    @test "hvg_degree_anomaly_stride4" in decoupled_names
    @test "nvg_degree_anomaly_stride8" in decoupled_names
    @test "nvg_degree_anomaly_stride4" ∉ decoupled_names

    @test_throws ArgumentError artifact_detection_methods(visibility_graph_sample_stride=0)
    @test_throws ArgumentError artifact_detection_methods(hvg_sample_strides=(1, 0))
end

@testset "Artifact score context methods" begin
    scores = [1.0, NaN, 3.0, 2.0]
    @test contextualize_artifact_scores(scores; reducer=:max, radius=1) == [1.0, 3.0, 3.0, 3.0]
    @test contextualize_artifact_scores(scores; reducer=:mean, radius=1) == [1.0, 2.0, 2.5, 2.5]

    base_method = ArtifactDetectionMethod("toy", "Synthetic score.", (eeg, sample_ranges) -> scores)
    context_names = getfield.(contextual_artifact_detection_methods([base_method]), :name)
    @test "toy_ctxmax_w1" in context_names
    @test "toy_ctxmax_w2" in context_names
    @test "toy_ctxmean_w1" in context_names
    @test_throws ArgumentError contextualize_artifact_scores(scores; reducer=:median, radius=1)
    @test_throws ArgumentError contextualize_artifact_scores(scores; reducer=:max, radius=-1)
end

@testset "Artifact logistic classifier" begin
    score_df = synthetic_artifact_score_df()
    feature_df = artifact_score_feature_table(score_df)
    @test nrow(feature_df) == 18
    @test Set(artifact_feature_names(feature_df)) == Set([:method_a, :method_b])

    results = leave_one_patient_out_artifact_classifier(score_df; n_thresholds=20)
    @test nrow(results.summary_df) == 1
    @test nrow(results.patient_df) == 3
    @test nrow(results.score_df) == 18
    @test all(isfinite, results.score_df.score)
    @test all(0 .<= results.score_df.score .<= 1)
    @test only(results.summary_df.f1) >= 0.8
    @test isfinite(only(results.summary_df.auprc))
end

@testset "Whole-recording HVG node features" begin
    node_features = whole_hvg_node_features([2.0, 1.0, 3.0])

    @test node_features.feature_values[:hvg_degree] == [2.0, 2.0, 2.0]
    @test node_features.feature_values[:hvg_outdegree] == [2.0, 1.0, 0.0]
    @test node_features.feature_values[:hvg_indegree] == [0.0, 1.0, 2.0]
    @test node_features.feature_values[:hvg_degree_imbalance] == [2.0, 0.0, -2.0]
    @test all(isfinite, node_features.feature_values[:hvg_random_law_surprisal])
    @test all(isfinite, node_features.node_anomaly)

    centrality_features = whole_hvg_node_features(
        [2.0, 1.0, 3.0, 2.0, 4.0];
        include_edge_span_features=true,
        include_centrality_features=true,
        centrality_landmark_count=3,
        centrality_max_nodes=3
    )
    for feature_name in (
            :hvg_edge_span_mean,
            :hvg_edge_span_max,
            :hvg_pagerank,
            :hvg_eigenvector_centrality,
            :hvg_sampled_betweenness,
            :hvg_landmark_closeness
        )
        @test haskey(centrality_features.feature_values, feature_name)
        @test length(centrality_features.feature_values[feature_name]) == 5
        @test all(isfinite, centrality_features.feature_values[feature_name])
    end
end

@testset "Within-bin whole-HVG graph features" begin
    pairs = Dict(whole_hvg_within_bin_graph_feature_pairs([2.0, 1.0, 3.0], 1:3))

    @test pairs[:hvg_bin_vertex_count] == 3.0
    @test pairs[:hvg_bin_edge_count] == 3.0
    @test pairs[:hvg_bin_edge_density] == 1.0
    @test pairs[:hvg_bin_mean_degree] == 2.0
    @test pairs[:hvg_bin_degree_std] == 0.0
    @test pairs[:hvg_bin_degree_max] == 2.0
    @test pairs[:hvg_bin_degree_q95] == 2.0
    @test pairs[:hvg_bin_degree_entropy] == 0.0
    @test isfinite(pairs[:hvg_bin_random_law_surprisal_mean])
    @test pairs[:hvg_bin_edge_span_mean] ≈ 4 / 3
    @test pairs[:hvg_bin_edge_span_max] == 2.0
    @test isfinite(pairs[:hvg_bin_edge_span_q95])
    @test pairs[:hvg_bin_global_clustering_coefficient] == 1.0
    @test isnan(pairs[:hvg_bin_degree_assortativity])

    one_node_pairs = Dict(whole_hvg_within_bin_graph_feature_pairs([2.0], 1:1))
    @test one_node_pairs[:hvg_bin_vertex_count] == 1.0
    @test one_node_pairs[:hvg_bin_edge_count] == 0.0
    @test isnan(one_node_pairs[:hvg_bin_edge_density])
end

@testset "Whole-recording HVG bin aggregation and labels" begin
    eeg = WholeHVGTestEEG(
        [
            2.0 1.0 3.0 2.0 1.0 3.0
            1.0 2.0 1.0 3.0 1.0 2.0
        ],
        1,
        6.0,
        ["C1", "C2"],
        [(4.0, 5.0)]
    )

    result = whole_hvg_bin_score_rows(eeg;
        patient=99,
        artifact_bounds=[(1.0, 2.0)],
        bin_s=3,
        recording_duration=6.0,
        include_edge_span_features=true,
        include_centrality_features=true,
        include_within_bin_graph_features=true,
        centrality_landmark_count=3
    )
    df = result.bin_df

    @test nrow(df) == 4
    @test Set(df.channel_label) == Set(["C1", "C2"])
    @test all(df.sample_count .== 3)
    @test :hvg_degree_mean in propertynames(df)
    @test :hvg_degree_abs_robust_z_q95 in propertynames(df)
    @test :hvg_edge_span_mean_mean in propertynames(df)
    @test :hvg_edge_span_max_abs_robust_z_q95 in propertynames(df)
    @test :hvg_pagerank_mean in propertynames(df)
    @test :hvg_pagerank_abs_robust_z_q95 in propertynames(df)
    @test :hvg_sampled_betweenness_q95 in propertynames(df)
    @test :hvg_landmark_closeness_abs_robust_z_max in propertynames(df)
    @test :hvg_node_anomaly_max in propertynames(df)
    @test :hvg_bin_vertex_count in propertynames(df)
    @test :hvg_bin_edge_density in propertynames(df)
    @test :hvg_bin_degree_entropy in propertynames(df)
    @test :hvg_bin_edge_span_q95 in propertynames(df)
    @test :hvg_bin_global_clustering_coefficient in propertynames(df)
    @test :hvg_bin_degree_assortativity in propertynames(df)
    @test all(df.hvg_bin_vertex_count .== 3.0)
    @test all(isfinite, df.hvg_bin_edge_density)

    first_bin = df[df.bin_index .== 1, :]
    second_bin = df[df.bin_index .== 2, :]
    @test all(first_bin.artifact_truth)
    @test !any(first_bin.seizure_truth)
    @test !any(first_bin.clean_seizure_truth)
    @test !any(second_bin.artifact_truth)
    @test all(second_bin.seizure_truth)
    @test all(second_bin.clean_seizure_truth)
    @test all(second_bin.seizure_event_id .== 1)
end

@testset "Whole-recording HVG metric notes" begin
    notes = whole_hvg_metric_notes()
    @test only(notes[notes.metric .== "connectedness_components", :status]) == "skipped"
    @test only(notes[notes.metric .== "nvg_metrics", :status]) == "deferred"
    @test only(notes[notes.metric .== "hvg_edge_span", :status]) == "optional_skipped"
    @test only(notes[notes.metric .== "hvg_pagerank", :status]) == "optional_skipped"
    @test only(notes[notes.metric .== "hvg_sampled_betweenness", :status]) == "optional_skipped"
    @test only(notes[notes.metric .== "hvg_bin_graph_metrics", :status]) == "optional_skipped"
    @test only(notes[notes.metric .== "hvg_bin_graph_topology_metrics", :status]) == "optional_skipped"

    edge_notes = whole_hvg_metric_notes(
        include_edge_span_features=true,
        include_centrality_features=true,
        include_within_bin_graph_features=true
    )
    @test only(edge_notes[edge_notes.metric .== "hvg_edge_span", :status]) == "included_optional"
    @test only(edge_notes[edge_notes.metric .== "hvg_pagerank", :status]) == "included_optional"
    @test only(edge_notes[edge_notes.metric .== "hvg_landmark_closeness", :status]) == "included_optional"
    @test only(edge_notes[edge_notes.metric .== "hvg_bin_graph_metrics", :status]) == "included_optional"
    @test only(edge_notes[edge_notes.metric .== "hvg_bin_graph_topology_metrics", :status]) == "included_optional"

    cheap_bin_notes = whole_hvg_metric_notes(
        include_within_bin_graph_features=true,
        include_within_bin_topology_features=false
    )
    @test only(cheap_bin_notes[cheap_bin_notes.metric .== "hvg_bin_graph_metrics", :status]) == "included_optional"
    @test only(cheap_bin_notes[cheap_bin_notes.metric .== "hvg_bin_graph_topology_metrics", :status]) == "optional_skipped"
end
