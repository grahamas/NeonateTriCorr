using Test

using DrWatson
@quickactivate "NeonateTriCorr"

include(srcdir("types.jl"))
include(srcdir("detect_artifacts.jl"))

struct ArtifactScoreTestEEG <: AbstractProcessedEEG
    signals::Matrix{Float64}
    sample_rate::Int
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
