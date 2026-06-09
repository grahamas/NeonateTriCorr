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

const EEG_SIGNAL_CLUSTER_ID_COLUMNS = [
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
]

const EEG_SIGNAL_CLUSTER_ID_COLUMN_SET = Set(EEG_SIGNAL_CLUSTER_ID_COLUMNS)

function eeg_signal_feature_name(idx, n_features)
    Symbol("eeg_signal_block_$(lpad(string(idx), ndigits(n_features), '0'))")
end

function eeg_signal_feature_names(n_features)
    [eeg_signal_feature_name(idx, n_features) for idx in 1:n_features]
end

function eeg_signal_channel_label(eeg::AbstractProcessedEEG, channel)
    if hasproperty(eeg, :labels) && channel <= length(eeg.labels)
        return eeg.labels[channel]
    end
    "channel_$(channel)"
end

function eeg_signal_sample_index_bounds(signal_length, sample_range)
    if isempty(sample_range)
        return missing, missing
    end
    first_idx = max(1, first(sample_range))
    last_idx = min(signal_length, last(sample_range))
    if first_idx > last_idx
        return missing, missing
    end
    first_idx, last_idx
end

function eeg_signal_sample_count_for_range(signal_length, sample_range)
    first_idx, last_idx = eeg_signal_sample_index_bounds(signal_length, sample_range)
    if ismissing(first_idx) || ismissing(last_idx)
        return 0
    end
    last_idx - first_idx + 1
end

function eeg_signal_block_mean_vector(
        segment::AbstractVector,
        n_features;
        center_bin=false,
        scale_bin=false
    )
    n_features >= 1 || throw(ArgumentError("n_features must be positive"))
    if isempty(segment)
        return fill(NaN, n_features)
    end

    values = Float64.(segment)
    if center_bin || scale_bin
        center = center_bin ? median(values) : 0.0
        scale = 1.0
        if scale_bin
            _, scale = robust_center_and_scale(values)
        end
        values .= (values .- center) ./ scale
    end

    output = Vector{Float64}(undef, n_features)
    n = length(values)
    for feature_idx in 1:n_features
        lo = floor(Int, (feature_idx - 1) * n / n_features) + 1
        hi = floor(Int, feature_idx * n / n_features)
        if hi < lo
            nearest_idx = clamp(round(Int, (feature_idx - 0.5) * n / n_features), 1, n)
            output[feature_idx] = values[nearest_idx]
        else
            output[feature_idx] = mean(@view values[lo:hi])
        end
    end
    output
end

function eeg_signal_bin_rows(
        eeg::AbstractProcessedEEG;
        patient,
        artifact_bounds=Tuple{Float64,Float64}[],
        seizure_bounds=eeg.seizure_annotations,
        bin_s=15,
        n_features=128,
        center_bins=false,
        scale_bins=false,
        recording_duration=nothing
    )
    recording_duration = isnothing(recording_duration) ? eeg.duration : recording_duration
    bin_starts, bin_stops = artifact_bin_bounds(recording_duration, bin_s)
    sample_ranges = artifact_sample_ranges(eeg, bin_starts, bin_stops)
    artifact_truth = bounds_to_bin_labels(artifact_bounds, bin_starts, bin_stops)
    seizure_truth = bounds_to_bin_labels(seizure_bounds, bin_starts, bin_stops)
    seizure_event_ids = seizure_event_ids_for_bins(seizure_bounds, bin_starts, bin_stops)
    feature_names = eeg_signal_feature_names(n_features)

    rows = NamedTuple[]
    signals = getproperty(eeg, :signals)
    signal_length = size(signals, 2)

    for channel in axes(signals, 1)
        channel_label = eeg_signal_channel_label(eeg, channel)
        signal = vec(signals[channel, :])
        for bin_idx in eachindex(sample_ranges)
            sample_start_idx, sample_stop_idx =
                eeg_signal_sample_index_bounds(signal_length, sample_ranges[bin_idx])
            sample_count = eeg_signal_sample_count_for_range(signal_length, sample_ranges[bin_idx])
            segment = sample_count == 0 ?
                Float64[] :
                signal[sample_start_idx:sample_stop_idx]
            feature_values = eeg_signal_block_mean_vector(segment, n_features;
                center_bin=center_bins,
                scale_bin=scale_bins
            )

            pairs = Pair{Symbol,Any}[
                :patient => patient,
                :channel => channel,
                :channel_label => channel_label,
                :bin_index => bin_idx,
                :bin_start => bin_starts[bin_idx],
                :bin_stop => bin_stops[bin_idx],
                :sample_start_index => sample_start_idx,
                :sample_stop_index => sample_stop_idx,
                :sample_count => sample_count,
                :artifact_truth => artifact_truth[bin_idx],
                :seizure_truth => seizure_truth[bin_idx],
                :clean_seizure_truth => seizure_truth[bin_idx] && !artifact_truth[bin_idx],
                :seizure_event_id => seizure_event_ids[bin_idx]
            ]
            for (feature_name, feature_value) in zip(feature_names, feature_values)
                push!(pairs, feature_name => feature_value)
            end
            push!(rows, (; pairs...))
        end
    end

    DataFrame(rows)
end

function eeg_signal_bin_rows_for_patient(
        patient_num;
        bin_s=15,
        artifact_grades=[1, 2],
        min_reviewers_per_seizure=3,
        artifact_csv_path=scriptsdir("helsinki_artifacts.csv"),
        n_features=128,
        center_bins=false,
        scale_bins=false
    )
    eeg = load_helsinki_eeg(patient_num;
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        excluded_artifact_grades=Int[]
    )
    artifact_bounds = load_helsinki_artifact_ground_truth(patient_num, artifact_grades;
        csv_path=artifact_csv_path,
        recording_duration=eeg.duration
    )

    eeg_signal_bin_rows(eeg;
        patient=patient_num,
        artifact_bounds=artifact_bounds,
        seizure_bounds=eeg.seizure_annotations,
        bin_s=bin_s,
        n_features=n_features,
        center_bins=center_bins,
        scale_bins=scale_bins
    )
end

function run_eeg_signal_feature_survey(;
        patients=artifact_labeled_patients(),
        bin_s=15,
        artifact_grades=[1, 2],
        min_reviewers_per_seizure=3,
        artifact_csv_path=scriptsdir("helsinki_artifacts.csv"),
        n_features=128,
        center_bins=false,
        scale_bins=false
    )
    bin_dfs = DataFrame[]
    event_offset = 0

    for patient in patients
        @info "Building EEG signal bin vectors" patient
        patient_df = eeg_signal_bin_rows_for_patient(patient;
            bin_s=bin_s,
            artifact_grades=artifact_grades,
            min_reviewers_per_seizure=min_reviewers_per_seizure,
            artifact_csv_path=artifact_csv_path,
            n_features=n_features,
            center_bins=center_bins,
            scale_bins=scale_bins
        )

        local_event_count = nrow(patient_df) == 0 ? 0 : maximum(patient_df.seizure_event_id)
        patient_df.seizure_event_id = offset_event_ids(patient_df, event_offset)
        event_offset += local_event_count
        push!(bin_dfs, patient_df)
    end

    isempty(bin_dfs) ? DataFrame() : reduce(vcat, bin_dfs)
end

function numeric_eeg_signal_feature_names(df::AbstractDataFrame)
    names_out = Symbol[]
    for name in Symbol.(names(df))
        name in EEG_SIGNAL_CLUSTER_ID_COLUMN_SET && continue
        startswith(String(name), "eeg_signal_") || continue
        eltype(df[!, name]) <: Union{Missing,Real} || continue
        push!(names_out, name)
    end
    sort(names_out; by=String)
end

function eeg_signal_clustering_feature_matrix(df::AbstractDataFrame, feature_names)
    X = Matrix{Float64}(undef, nrow(df), length(feature_names))
    for (feature_idx, feature_name) in enumerate(feature_names)
        values = df[!, feature_name]
        for row_idx in 1:nrow(df)
            X[row_idx, feature_idx] = artifact_feature_value(values[row_idx])
        end
    end
    X
end

function empty_eeg_signal_skipped_feature_df()
    DataFrame(feature=String[], reason=String[])
end

function eeg_signal_standardized_feature_matrix(X::AbstractMatrix, feature_names; min_scale=1e-8, z_clip=20.0)
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

    isempty(keep_idxs) && throw(ArgumentError("No finite, nonconstant EEG signal columns available for clustering."))

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

    skipped_df = isempty(skipped_rows) ? empty_eeg_signal_skipped_feature_df() : DataFrame(skipped_rows)
    (
        Z=Z,
        feature_names=kept_features,
        centers=centers,
        scales=scales,
        skipped_df=skipped_df
    )
end

function eeg_signal_pca_embedding(Z::AbstractMatrix; outdim=2)
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

function eeg_signal_embedding_sample_indices(
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

function eeg_signal_tsne_embedding(
        Z::AbstractMatrix;
        outdim=2,
        reduce_dims=min(30, size(Z, 2)),
        iterations=1000,
        perplexity=30.0,
        seed=1,
        pca_init=true
    )
    Random.seed!(seed)
    effective_perplexity = min(perplexity, max(1.0, (size(Z, 1) - 1) / 3))
    tsne(Matrix{Float64}(Z), outdim, reduce_dims, iterations, effective_perplexity;
        pca_init=pca_init,
        progress=false,
        verbose=false
    )
end

function eeg_signal_distribution_values(df::AbstractDataFrame, feature_name::Symbol)
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

function eeg_signal_feature_distribution_rows(df::AbstractDataFrame, feature_names, centers, scales)
    rows = NamedTuple[]
    for (idx, feature_name) in enumerate(feature_names)
        finite_values, invalid_count = eeg_signal_distribution_values(df, feature_name)
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

function eeg_signal_plain_float(value)
    if value isa Real && isfinite(value)
        return string(round(Float64(value); sigdigits=6))
    end
    string(value)
end

function write_eeg_signal_feature_distribution_report(
        path,
        distribution_df::AbstractDataFrame,
        skipped_df::AbstractDataFrame;
        method,
        n_features,
        bin_s,
        center_bins,
        scale_bins,
        sample_df=DataFrame()
    )
    open(path, "w") do io
        println(io, "EEG Signal $(uppercase(String(method))) Feature Report")
        println(io, "Generated: $(Dates.now())")
        println(io)
        println(io, "Feature vector: one patient-channel-bin row represented by $(n_features) block-mean EEG amplitudes.")
        println(io, "Bin width seconds: $(bin_s)")
        println(io, "Bin centering: $(center_bins)")
        println(io, "Bin scaling: $(scale_bins)")
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
                "min=$(eeg_signal_plain_float(row.min)), p05=$(eeg_signal_plain_float(row.p05)), " *
                "q25=$(eeg_signal_plain_float(row.q25)), median=$(eeg_signal_plain_float(row.median)), " *
                "mean=$(eeg_signal_plain_float(row.mean)), q75=$(eeg_signal_plain_float(row.q75)), " *
                "p95=$(eeg_signal_plain_float(row.p95)), max=$(eeg_signal_plain_float(row.max)), " *
                "std=$(eeg_signal_plain_float(row.std)), robust_center=$(eeg_signal_plain_float(row.robust_center)), " *
                "robust_scale=$(eeg_signal_plain_float(row.robust_scale))"
            )
        end
    end
end

function eeg_signal_embedding_df(
        df::AbstractDataFrame;
        method=:tsne,
        tsne_max_rows=5000,
        tsne_seed=1,
        tsne_iterations=1000,
        tsne_perplexity=30.0,
        n_features=128,
        bin_s=15,
        center_bins=false,
        scale_bins=false
    )
    feature_names = numeric_eeg_signal_feature_names(df)
    raw_X = eeg_signal_clustering_feature_matrix(df, feature_names)
    standardized = eeg_signal_standardized_feature_matrix(raw_X, feature_names)

    sample_idxs = if method == :tsne
        eeg_signal_embedding_sample_indices(df;
            max_rows=tsne_max_rows,
            seed=tsne_seed
        )
    else
        collect(1:nrow(df))
    end

    Z_embedding = standardized.Z[sample_idxs, :]
    coords = if method == :pca
        eeg_signal_pca_embedding(Z_embedding; outdim=2)
    elseif method == :tsne
        eeg_signal_tsne_embedding(Z_embedding;
            outdim=2,
            iterations=tsne_iterations,
            perplexity=tsne_perplexity,
            seed=tsne_seed
        )
    else
        throw(ArgumentError("Unsupported EEG signal clustering method: $(method)"))
    end

    present_id_columns = [column for column in EEG_SIGNAL_CLUSTER_ID_COLUMNS if column in Symbol.(names(df))]
    embedding_df = select(DataFrame(df[sample_idxs, :]), present_id_columns...)
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
            resampled_signal_features=n_features,
            bin_s=bin_s,
            center_bins=center_bins,
            scale_bins=scale_bins,
            tsne_max_rows=method == :tsne ? tsne_max_rows : missing,
            tsne_seed=method == :tsne ? tsne_seed : missing,
            tsne_iterations=method == :tsne ? tsne_iterations : missing,
            tsne_perplexity=method == :tsne ? tsne_perplexity : missing
        ),
        distribution_df=eeg_signal_feature_distribution_rows(
            df,
            standardized.feature_names,
            standardized.centers,
            standardized.scales
        )
    )
end

function eeg_signal_label_colors(labels; false_color=(:gray60, 0.35), true_color=(:red, 0.75))
    [Bool(label) ? true_color : false_color for label in labels]
end

function plot_eeg_signal_embedding(
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

    colors = eeg_signal_label_colors(embedding_df[!, color_column];
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

function save_eeg_signal_clustering_outputs(
        bin_df::AbstractDataFrame;
        output_dir=datadir("exp_pro", "eeg_signal_clustering", string(Dates.now())),
        plot_dir=plotsdir("eeg_signal_clustering_$(Dates.now())"),
        method=:tsne,
        tsne_max_rows=5000,
        tsne_seed=1,
        tsne_iterations=1000,
        tsne_perplexity=30.0,
        n_features=128,
        bin_s=15,
        center_bins=false,
        scale_bins=false
    )
    mkpath(output_dir)
    mkpath(plot_dir)

    result = eeg_signal_embedding_df(bin_df;
        method=method,
        tsne_max_rows=tsne_max_rows,
        tsne_seed=tsne_seed,
        tsne_iterations=tsne_iterations,
        tsne_perplexity=tsne_perplexity,
        n_features=n_features,
        bin_s=bin_s,
        center_bins=center_bins,
        scale_bins=scale_bins
    )
    embedding_df = result.embedding_df
    CSV.write(joinpath(output_dir, "eeg_signal_embedding.csv"), embedding_df)
    CSV.write(joinpath(output_dir, "eeg_signal_embedding_sample.csv"), result.sample_df)
    CSV.write(
        joinpath(output_dir, "eeg_signal_embedding_features.csv"),
        DataFrame(feature=String.(result.feature_names), center=result.centers, scale=result.scales)
    )
    CSV.write(joinpath(output_dir, "eeg_signal_embedding_skipped_features.csv"), result.skipped_df)
    CSV.write(joinpath(output_dir, "eeg_signal_embedding_feature_distributions.csv"), result.distribution_df)
    write_eeg_signal_feature_distribution_report(
        joinpath(output_dir, "eeg_signal_embedding_feature_report.txt"),
        result.distribution_df,
        result.skipped_df;
        method=method,
        n_features=n_features,
        bin_s=bin_s,
        center_bins=center_bins,
        scale_bins=scale_bins,
        sample_df=result.sample_df
    )

    artifact_fig = plot_eeg_signal_embedding(embedding_df;
        color_column=:artifact_truth,
        title="EEG signal $(uppercase(String(method))) embedding by artifact truth",
        true_label="artifact",
        true_color=:red
    )
    save(joinpath(plot_dir, "eeg_signal_embedding_artifact_truth.png"), artifact_fig)

    seizure_fig = plot_eeg_signal_embedding(embedding_df;
        color_column=:seizure_truth,
        title="EEG signal $(uppercase(String(method))) embedding by seizure truth",
        true_label="seizure",
        true_color=:dodgerblue
    )
    save(joinpath(plot_dir, "eeg_signal_embedding_seizure_truth.png"), seizure_fig)

    if :clean_seizure_truth in Symbol.(names(embedding_df))
        clean_seizure_fig = plot_eeg_signal_embedding(embedding_df;
            color_column=:clean_seizure_truth,
            title="EEG signal $(uppercase(String(method))) embedding by clean seizure truth",
            true_label="clean seizure",
            true_color=:purple
        )
        save(joinpath(plot_dir, "eeg_signal_embedding_clean_seizure_truth.png"), clean_seizure_fig)
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

method = isdefined(Main, :EEG_SIGNAL_CLUSTER_METHOD) ? EEG_SIGNAL_CLUSTER_METHOD : :tsne
tsne_max_rows = isdefined(Main, :EEG_SIGNAL_TSNE_MAX_ROWS) ? EEG_SIGNAL_TSNE_MAX_ROWS : 5000
tsne_seed = isdefined(Main, :EEG_SIGNAL_TSNE_SEED) ? EEG_SIGNAL_TSNE_SEED : 1
tsne_iterations = isdefined(Main, :EEG_SIGNAL_TSNE_ITERATIONS) ? EEG_SIGNAL_TSNE_ITERATIONS : 1000
tsne_perplexity = isdefined(Main, :EEG_SIGNAL_TSNE_PERPLEXITY) ? EEG_SIGNAL_TSNE_PERPLEXITY : 30.0
patients = isdefined(Main, :EEG_SIGNAL_CLUSTER_PATIENTS) ?
    EEG_SIGNAL_CLUSTER_PATIENTS :
    artifact_labeled_patients()
bin_s = isdefined(Main, :EEG_SIGNAL_BIN_S) ? EEG_SIGNAL_BIN_S : 15
artifact_grades = isdefined(Main, :EEG_SIGNAL_ARTIFACT_GRADES) ? EEG_SIGNAL_ARTIFACT_GRADES : [1, 2]
min_reviewers_per_seizure = isdefined(Main, :EEG_SIGNAL_MIN_REVIEWERS_PER_SEIZURE) ?
    EEG_SIGNAL_MIN_REVIEWERS_PER_SEIZURE :
    3
n_features = isdefined(Main, :EEG_SIGNAL_RESAMPLED_FEATURES) ? EEG_SIGNAL_RESAMPLED_FEATURES : 128
center_bins = isdefined(Main, :EEG_SIGNAL_CENTER_BINS) ? EEG_SIGNAL_CENTER_BINS : false
scale_bins = isdefined(Main, :EEG_SIGNAL_SCALE_BINS) ? EEG_SIGNAL_SCALE_BINS : false

bin_df = if isdefined(Main, :EEG_SIGNAL_CLUSTER_INPUT_CSV)
    @info "Loading EEG signal bin vectors for clustering" EEG_SIGNAL_CLUSTER_INPUT_CSV
    CSV.read(EEG_SIGNAL_CLUSTER_INPUT_CSV, DataFrame)
else
    run_eeg_signal_feature_survey(;
        patients=patients,
        bin_s=bin_s,
        artifact_grades=artifact_grades,
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        n_features=n_features,
        center_bins=center_bins,
        scale_bins=scale_bins
    )
end

output_root = isdefined(Main, :EEG_SIGNAL_CLUSTER_OUTPUT_ROOT) ?
    EEG_SIGNAL_CLUSTER_OUTPUT_ROOT :
    datadir("exp_pro", "eeg_signal_clustering", string(Dates.now()))
plot_root = isdefined(Main, :EEG_SIGNAL_CLUSTER_PLOT_ROOT) ?
    EEG_SIGNAL_CLUSTER_PLOT_ROOT :
    plotsdir("eeg_signal_clustering_$(Dates.now())")

eeg_signal_clustering_results = save_eeg_signal_clustering_outputs(bin_df;
    output_dir=output_root,
    plot_dir=plot_root,
    method=method,
    tsne_max_rows=tsne_max_rows,
    tsne_seed=tsne_seed,
    tsne_iterations=tsne_iterations,
    tsne_perplexity=tsne_perplexity,
    n_features=n_features,
    bin_s=bin_s,
    center_bins=center_bins,
    scale_bins=scale_bins
)

eeg_signal_clustering_results
