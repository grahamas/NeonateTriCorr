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
    artifact[div(samples_per_bin, 2)] = 12.0
    signal = vcat(base, base, base, artifact)
    ArtifactScoreTestEEG(reshape(signal, 1, :), sample_rate)
end

@testset "Visibility graph artifact scores" begin
    eeg = repeated_window_eeg()
    sample_ranges = [1:32, 33:64, 65:96, 97:128]

    hvg_scores = score_hvg_degree_anomaly(eeg, sample_ranges; sample_stride=1)
    nvg_scores = score_nvg_degree_anomaly(eeg, sample_ranges; sample_stride=1)

    @test length(hvg_scores) == length(sample_ranges)
    @test length(nvg_scores) == length(sample_ranges)
    @test all(isfinite, hvg_scores)
    @test all(isfinite, nvg_scores)
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
    @test_throws ArgumentError artifact_detection_methods(visibility_graph_sample_stride=0)
end
