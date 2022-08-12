
function evaluate_detection_posseizure_negepoch(truth_bounds::Vector{<:Tuple}, detections::Vector)
    # One positive per seizure, one negative per non-seizure epoch
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    valid_truth_bounds = [(on,off) for (on,off) ∈ truth_bounds if !all(ismissing.(detections[on:off]))]
    for (on, off) ∈ valid_truth_bounds # should combine adjacent seizures
        detected_seizure = any(skipmissing(detections[on:off]))
        true_positives += detected_seizure
        false_negatives += !detected_seizure
        detections[on:off] .= missing # should exclude adjacent bins
    end
    true_negatives = count(skipmissing(detections) .== false)
    false_positives = count(skipmissing(detections) .== true)
    return (
        gt_positive = length(valid_truth_bounds),
        gt_negative = count(.!ismissing.(detections)), # should combine seconds
        true_positives = true_positives,
        false_positives = false_positives,
        true_negatives = true_negatives,
        false_negatives = false_negatives
    )
end

function add_grace_to_truth_bounds(truth_bounds, time_idxs, grace_s)
    if isempty(truth_bounds)
        return eltype(truth_bounds)[]
    end
    min_idx = Int(minimum(time_idxs)+1)
    max_idx = Int(maximum(time_idxs)+1)
    truth_bounds_with_grace_potentially_overlapping = [(max(min_idx,on-grace_s),min(max_idx,off+grace_s)) for (on,off) ∈ truth_bounds]
    truth_bounds_with_grace = eltype(truth_bounds_with_grace_potentially_overlapping)[]
    prev_bounds = first(truth_bounds_with_grace_potentially_overlapping)
    for bounds ∈ truth_bounds_with_grace_potentially_overlapping[begin+1:end]
        if prev_bounds[end] >= bounds[begin]
            merged_bounds = (prev_bounds[begin], bounds[end])
            prev_bounds = merged_bounds
        else
            push!(truth_bounds_with_grace, prev_bounds)
            prev_bounds = bounds
        end
    end
    push!(truth_bounds_with_grace, prev_bounds)
    return truth_bounds_with_grace
end

function evaluate_detection_posseizure_negalerts(truth_bounds::Vector{<:Tuple}, detections::Vector, detections_times; snippets_duration_s, alert_grace_s=60)
    # One positive per seizure, one negative per alert period
    alert_grace = alert_grace_s ÷ snippets_duration_s
    detections = Vector{Union{Int,Missing}}(detections)
    len = length(detections)
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    valid_truth_bounds = [(on,off) for (on,off) ∈ truth_bounds if !all(ismissing.(detections[on .<= detections_times .< off]))]
    truth_bounds_with_grace = add_grace_to_truth_bounds(valid_truth_bounds, detections_times, alert_grace_s)
    gt_positive = length(truth_bounds_with_grace)
    for (on, off) ∈ truth_bounds_with_grace # should combine adjacent seizures
        detected_seizure = any(skipmissing(detections[on .<= detections_times .< off]) .== 1)
        true_positives += detected_seizure
        false_negatives += !detected_seizure
        detections[on .<= detections_times .< off] .= missing # should exclude adjacent bins
    end
    for i ∈ eachindex(detections)
        detect = detections[i]
        if ismissing(detect)
            continue
        elseif detect == 0
            true_negatives += 1
            detections[min(i+1,len):min(i+alert_grace,len)] .+= 2
        elseif (detect == 1) || (detect == 3)
            false_positives += 1
            true_negatives -= 1
            detections[min(i+1,len):min(i+alert_grace,len)] .-= 10
        end
    end
    return (
        gt_positive = gt_positive,
        gt_negative = true_negatives + false_positives,
        true_positives = true_positives,
        false_positives = false_positives,
        true_negatives = true_negatives,
        false_negatives = false_negatives
    )
end

function apply_rolling_deviation_window(signals, window_fn, window_len)
    n_signals = size(signals, 1)
    window_red = mapreduce(hcat, 1:n_signals) do signal_num
        windowed = [window_fn(signals[signal_num,(i-window_len+1):(i)]) for i ∈ window_len:(size(signals,2))]
        (windowed .- mean(skipmissing(windowed))) ./ std(skipmissing(windowed))
    end
    vcat(zeros(eltype(window_red), window_len-1, n_signals) .+ missing, window_red)
end

function detect_deviation_window(signals, window_fn, window_len, motifs_criterion)
    rolling_deviation = apply_rolling_deviation_window(signals, window_fn, window_len)
    detections = map(motifs_criterion, eachslice(rolling_deviation, dims=1))
    return detections
end

function motif_weights_based_on_pvalues(df, effect_sym, n_nonzero, n_signals)
    if isempty(df)
        return ones(n_signals) / n_signals
    end
    descending_significance = sort(zip(df.p, df[:, effect_sym], 1:n_signals) |> collect)
    weights = zeros(n_signals)
    weights[descending_significance[1:n_nonzero] .|> x -> x[3]] .=(descending_significance[1:n_nonzero] .|> x -> sign(x[2]))
    return weights / norm(weights)
end

function motif0_weight_based_effect(df, effect_sym, n_nonzero, n_signals)
    if isempty(df)
        weights = zeros(n_signals)
        weights[1] = 1
        return weights
    end
    weights = zeros(n_signals)
    weights[1] = sign(df[1, effect_sym])
    return weights
end


function calculate_minimal_threshold_range(arr, arr_times, truth_bounds; alert_grace_s)
    min_θ = Inf
    max_θ = -Inf
    truth_bounds_with_grace = add_grace_to_truth_bounds(truth_bounds, arr_times, alert_grace_s)
    for bounds ∈ truth_bounds_with_grace
        non_artifact = skipmissing(arr[bounds[1] .<= arr_times .< bounds[2]])
        if isempty(non_artifact)
            continue
        end
        bounds_min, bounds_max = extrema(non_artifact)
        min_θ = min(min_θ, bounds_min)
        max_θ = max(max_θ, bounds_max)
    end
    # below min_θ and above max_θ, no seizures are detected
    return (min_θ, max_θ)
end

function calculate_threshold_range(arr, arr_times, truth_bounds; alert_grace_s)
    extrema(skipmissing(arr))
end

function calculate_ROC(full_detection_fn::Function, θs)
    rocnums = map(θs) do θ
        rocnum = full_detection_fn(θ)
        merge((θ=θ,), rocnum)
    end
    DataFrame(rocnums)
end

function calculate_ROC(test_results_df::DataFrame, raw_signals::AbstractArray, signal_times, truth_bounds; signal_sym, window_fn, rolling_window_s, alert_grace_s, snippets_duration_s,  n_signals_used=5, n_θs=100, motif_weights_fn=motif_weights_based_on_pvalues)
    window_len = rolling_window_s ÷ snippets_duration_s
    n_signals = size(raw_signals,1)
    max_significance_weights = motif_weights_fn(test_results_df, signal_sym, n_signals_used, n_signals)
    rolled_signal = apply_rolling_deviation_window(raw_signals, window_fn, window_len) * max_significance_weights
    min_θ, max_θ = calculate_threshold_range(rolled_signal, signal_times,truth_bounds; alert_grace_s=alert_grace_s)
    if !isfinite(min_θ)
        return nothing
    end

    function full_detection_fn(θ)
        evaluate_detection_posseizure_negalerts(truth_bounds, rolled_signal .>= θ, signal_times; snippets_duration_s=snippets_duration_s, alert_grace_s=alert_grace_s)
    end

    calculate_ROC(full_detection_fn, range(min_θ, max_θ, length=n_θs))
end

function plot_μ_and_σ_signals_and_roc(args...; resolution=(1300, 1000), kwargs...)
    fig = Figure(resolution=resolution)
    plot_μ_and_σ_signals_and_roc!(fig, args...; kwargs...)
    return fig
end


function plot_μ_and_σ_signals_and_roc!(fig, results_df, signals, signal_times, truth_bounds; eeg, rolling_window_s, example_θ, n_signals_used=5, snippets_duration_s, alert_grace_s, title)

    rolling_window = rolling_window_s ÷ snippets_duration_s

    n_signals = size(signals, 1)
    max_significance_μ_weights = motif_weights_based_on_pvalues(results_df, :Δμ, n_signals_used, n_signals)
    detections_mean = apply_rolling_deviation_window(signals, mean, rolling_window) * max_significance_μ_weights

    max_significance_σ_weights = motif_weights_based_on_pvalues(results_df, :Δσ, n_signals_used, n_signals)
    detections_std = apply_rolling_deviation_window(signals, std, rolling_window) * max_significance_σ_weights

    fig[1,1:2] = ax_mean = Axis(fig; ylabel = "mean μ (five most significant channels)")
    fig[2,1:2] = ax_std = Axis(fig; ylabel = "mean σ (five most significant channels)")
    fig[3,1:2] = ax_rev = Axis(fig; ylabel = "# reviewers", xlabel = "time")
    TriCorrApplications.plot_reviewer_consensus!(ax_rev, eeg)
    TriCorrApplications.plot_contribution!(ax_mean, eeg, signal_times, detections_mean)
    TriCorrApplications.plot_contribution!(ax_std, eeg, signal_times, detections_std)
    hlines!(ax_mean, [example_θ], color=:red, linestyle=:dash)
    hlines!(ax_std, [example_θ], color=:red, linestyle=:dash)

    roc_column = fig[:,end+1]
    # save(joinpath(save_dir, "signals_patient$(PAT)_reviewers$(min_reviewers_per_seizure).png"), fig)


    @info "Calculating ROC curve..."
    signal_window_fns = Dict(
        "σ" => (:Δσ, std),
        "μ" => (:Δμ, mean)
    )
    seizures_with_grace_period = add_grace_to_truth_bounds(eeg.seizure_annotations, get_times(eeg, sample_rate=1/snippets_duration_s), alert_grace_s)
    non_seizure_hours = (eeg.duration + mapreduce((x) -> x[1] - x[2], +, seizures_with_grace_period, init=0) + mapreduce((x) -> x[1] - x[2], +, eeg.artifact_annotations, init=0)) / (60 * 60)
    duration_hours = eeg.duration / (60 * 60)
    roc_drws = map(enumerate(["μ", "σ"])) do (i, signal_type)
        signal_sym, window_fn = signal_window_fns[signal_type]
        roc_data = calculate_ROC(results_df, signals, signal_times, truth_bounds; signal_sym=signal_sym, window_fn=window_fn, rolling_window_s=rolling_window_s, alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, n_signals_used=n_signals_used, n_θs=100)
        @info "done. Now plotting."
        
        if !isnothing(roc_data)
            roc_plt = data(roc_data) * mapping((:false_positives,:gt_negative) => ((f, gt) -> f / non_seizure_hours) => "FP/Hour", (:true_positives,:gt_positive)=> ((t, gt) -> t / gt) => "TPR") * visual(Lines, color=:blue, linewidth=5)
            roc_drw = draw!(roc_column[i,1], roc_plt, axis=(title="$(title) ($(signal_type) signal)", limits=((0.,maximum(roc_data.gt_negative)/duration_hours),(0.,1.))))

            roc_drw
        end
    end
    return fig
end

function plot_μ_comparison(args...; resolution=(1300, 1000), kwargs...)
    fig = Figure(resolution=resolution)
    plot_μ_comparison!(fig, args...; kwargs...)
    return fig
end


function plot_μ_comparison!(fig, results_dfs, (aeeg_signal_times, aeeg_signals), (tricorr_signal_times, tricorr_signals), truth_bounds; eeg, rolling_window_s, example_θ, n_signals_used=5, aeeg_snippets_duration_s, tricorr_snippets_duration_s, alert_grace_s, title)

    aeeg_rolling_window = rolling_window_s ÷ aeeg_snippets_duration_s
    tricorr_rolling_window = rolling_window_s ÷ tricorr_snippets_duration_s

    aeeg_detections_mean, tricorr_detections_mean = map(zip(results_dfs, (aeeg_signals, tricorr_signals))) do (results_df, signals)
        n_signals = size(signals, 1)
        max_significance_μ_weights = motif_weights_based_on_pvalues(results_df, :Δμ, n_signals_used, n_signals)
        detections_mean = apply_rolling_deviation_window(signals, mean, rolling_window) * max_significance_μ_weights
    end

    fig[1,1:2] = ax_mean_aeeg = Axis(fig; ylabel = "mean aEEG μ")
    fig[2,1:2] = ax_mean_tricorr = Axis(fig; ylabel = "mean motif-class μ")
    fig[3,1:2] = ax_rev = Axis(fig; ylabel = "# reviewers", xlabel = "time")
    TriCorrApplications.plot_reviewer_consensus!(ax_rev, eeg)
    TriCorrApplications.plot_contribution!(ax_mean_aeeg, eeg, aeeg_signal_times, tricorr_detections_mean)
    TriCorrApplications.plot_contribution!(ax_mean_tricorr, eeg, tricorr_signal_times, aeeg_detections_mean)
    hlines!(ax_mean_aeeg, [example_θ], color=:red, linestyle=:dash)
    hlines!(ax_mean_tricorr, [example_θ], color=:red, linestyle=:dash)

    roc_column = fig[:,end+1]
    # save(joinpath(save_dir, "signals_patient$(PAT)_reviewers$(min_reviewers_per_seizure).png"), fig)


    @info "Calculating ROC curve..."
    signal_window_fns = Dict(
        "σ" => (:Δσ, std),
        "μ" => (:Δμ, mean)
    )
    df_names = ["aEEG", "TriCorr"]
    roc_drws = map(enumerate(results_dfs)) do (i, results_df)
        signal_type = "μ"
        signal_sym, window_fn = signal_window_fns[signal_type]
        roc_data = calculate_ROC(results_df, signals, signal_times, truth_bounds; signal_sym=signal_sym, window_fn=window_fn, rolling_window_s=rolling_window_s, alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, n_signals_used=n_signals_used, n_θs=100)
        @info "done. Now plotting."
        
        if !isnothing(roc_data)
            roc_plt = data(roc_data) * mapping((:false_positives,:gt_negative) => ((f, gt) -> f / gt) => "FPR", (:true_positives,:gt_positive)=> ((t, gt) -> t / gt) => "TPR") * visual(Lines, color=:blue, linewidth=5)
            roc_drw = draw!(roc_column[i,1], roc_plt, axis=(title="$(title) ($(signal_type) $(df_names[i]))", limits=((0.,1.),(0.,1.))))

            roc_drw
        end
    end
    return fig
end