
# function evaluate_detection_posseizure_negepoch(truth_bounds::Vector{<:Tuple}, detections::Vector)
#     # One positive per seizure, one negative per non-seizure epoch
#     true_positives = 0
#     false_positives = 0
#     true_negatives = 0
#     false_negatives = 0
#     valid_truth_bounds = [(on,off) for (on,off) ∈ truth_bounds if !all(ismissing.(detections[on:off]))]
#     for (on, off) ∈ valid_truth_bounds # should combine adjacent seizures
#         detected_seizure = any(skipmissing(detections[on:off]))
#         true_positives += detected_seizure
#         false_negatives += !detected_seizure
#         detections[on:off] .= missing # should exclude adjacent bins
#     end
#     true_negatives = count(skipmissing(detections) .== false)
#     false_positives = count(skipmissing(detections) .== true)
#     return (
#         gt_positive = length(valid_truth_bounds),
#         gt_negative = count(.!ismissing.(detections)), # should combine seconds
#         true_positives = true_positives,
#         false_positives = false_positives,
#         true_negatives = true_negatives,
#         false_negatives = false_negatives
#     )
# end

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
        if ismissing(detect) || (detect == 2) || (detect < 0)
            continue
        elseif detect == 0
            # new grace period; correct negative detection
            true_negatives += 1
            detections[min(i+1,len):min(i+alert_grace,len)] .+= 2
        elseif (detect == 1)
            # new grace period; incorrect positive detection
            false_positives += 1
            detections[min(i+1,len):min(i+alert_grace,len)] .-= 10
        elseif (detect == 3)
            # current grace period previously negative, now positive
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

function apply_rolling_deviation_window(signals::ARR, window_fn, window_len) where ARR
    n_signals = size(signals, 1)
    window_red = mapreduce(hcat, 1:n_signals) do signal_num
        windowed = [window_fn(signals[signal_num,(i-window_len+1):(i)]) for i ∈ window_len:(size(signals,2))]
        (windowed .- mean(skipmissing(windowed))) ./ std(skipmissing(windowed))
    end
    ARR(hcat(zeros(eltype(window_red), n_signals, window_len - 1) .+ missing, window_red'))
end

function detect_deviation_window(signals, window_fn, window_len, motifs_criterion)
    rolling_deviation = apply_rolling_deviation_window(signals, window_fn, window_len)
    detections = map(motifs_criterion, eachslice(rolling_deviation, dims=1))
    return detections
end

function motif_weights_based_on_pvalues(df, effect_sym, n_nonzero, n_signals)
    if isempty(df)
        @warn "Defaulting to all signals b/c no results w/ p-values for this patient"
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

function calculate_ROC(raw_signals::AbstractArray, signal_times, truth_bounds, process_signal_fn::Function; alert_grace_s, snippets_duration_s, n_θs=100, processing_args...)
    rolled_signal = process_signal_fn(raw_signals; snippets_duration_s=snippets_duration_s, processing_args...)
    min_θ, max_θ = calculate_threshold_range(rolled_signal, signal_times,truth_bounds; alert_grace_s=alert_grace_s)
    if !isfinite(min_θ)
        return nothing
    end

    function full_detection_fn(θ)
        evaluate_detection_posseizure_negalerts(truth_bounds, rolled_signal .>= θ, signal_times; snippets_duration_s=snippets_duration_s, alert_grace_s=alert_grace_s)
    end

    calculate_ROC(full_detection_fn, range(min_θ, max_θ, length=n_θs))
end

function reduce_signal_distance(raw_signals::NamedDimsArray; window_fn, rolling_window_s, snippets_duration_s, unused_params...)
    raw_signals = NamedDimsArray{(:_, :time)}(raw_signals)
    window_len = rolling_window_s ÷ snippets_duration_s
    rolled_signal = apply_rolling_deviation_window(raw_signals, window_fn, window_len)
    count_nonmissing_times = count(.!ismissing.(sum(rolled_signal, dims=1)))
    unmissing_rolled_signal = copy(rolled_signal)
    unmissing_rolled_signal[ismissing.(rolled_signal)] .= 0.
    mean_signal = sum(unmissing_rolled_signal, dims=2) ./ count_nonmissing_times
    return vec(sqrt.(sum((rolled_signal .- mean_signal) .^ 2, dims=1)))
end

function reduce_signal_meanall(raw_signals::NamedDimsArray; window_fn, rolling_window_s, snippets_duration_s, unused_params...)
    raw_signals = NamedDimsArray{(:_, :time)}(raw_signals)
    window_len = rolling_window_s ÷ snippets_duration_s
    rolled_signal = apply_rolling_deviation_window(raw_signals, window_fn, window_len)
    return vec(maximum(abs.(rolled_signal), dims=1))
end

function reduce_signal_maxany(raw_signals::NamedDimsArray; window_fn, rolling_window_s, snippets_duration_s, unused_params...)
    raw_signals = NamedDimsArray{(:_, :time)}(raw_signals)
    window_len = rolling_window_s ÷ snippets_duration_s
    rolled_signal = apply_rolling_deviation_window(raw_signals, window_fn, window_len)
    return vec(abs.(mean(rolled_signal, dims=1)))
end

function reduce_signal_meansignificant(raw_signals::NamedDimsArray; 
        signal_sym, window_fn, n_signals_used=5,
        motif_weights_fn=motif_weights_based_on_pvalues,
        rolling_window_s,
        results_df::DataFrame,
        snippets_duration_s, unused_params...
    )
    raw_signals = NamedDimsArray{(:_, :time)}(raw_signals)
    window_len = rolling_window_s ÷ snippets_duration_s
    n_signals = size(raw_signals,1)
    max_significance_weights = motif_weights_fn(results_df, signal_sym, n_signals_used, n_signals)
    rolled_signal = max_significance_weights' * apply_rolling_deviation_window(raw_signals, window_fn, window_len)
    return vec(rolled_signal)
end


function get_reduce_signals_fn(reduction_type)
    Dict(
        "maxany" => reduce_signal_maxany,
        "meansignificant" => reduce_signal_meansignificant,
        "meanall" => reduce_signal_meanall,
        "distance" => reduce_signal_distance
    )[reduction_type]
end


function plot_μ_and_σ_signals_and_roc(args...; resolution=(1300, 1000), kwargs...)
    fig = Figure(resolution=resolution)
    plot_μ_and_σ_signals_and_roc!(fig, args...; kwargs...)
    return fig
end

function plot_μ_and_σ_signals_and_roc!(fig, signals, signal_times, truth_bounds; analysis_eeg, plot_eeg=analysis_eeg, rolling_window_s, snippets_duration_s, alert_grace_s, title, signals_reduction_params, signals_reduction_name, unused_params...)
    reduce_signals_fn = get_reduce_signals_fn(signals_reduction_name)

    if !isempty(unused_params)
        @warn "Unused params: $unused_params"
    end
    
    seizures_with_grace_period = add_grace_to_truth_bounds(analysis_eeg.seizure_annotations, get_times(analysis_eeg, sample_rate=1/snippets_duration_s), alert_grace_s)
    seizure_and_artifact_bounds = merge_bounds(seizures_with_grace_period, analysis_eeg.artifact_annotations)
    non_seizure_hours = (analysis_eeg.duration + mapreduce((x) -> x[1] - x[2], +, seizure_and_artifact_bounds, init=0)) / (60 * 60)

    @info "Calculating ROC curve..."
    signal_window_fns = Dict(
        "σ" => (signal_sym=:Δσ, window_fn=std),
        "μ" => (signal_sym=:Δμ, window_fn=mean)
    )
    map(enumerate(["μ", "σ"])) do (i, rolling_type)
        rolling_reduction_params = signal_window_fns[rolling_type]
        reduced_signal = reduce_signals_fn(signals; rolling_window_s=rolling_window_s, snippets_duration_s=snippets_duration_s, rolling_reduction_params..., signals_reduction_params...)

        fig[i,1:2] = ax = Axis(fig; ylabel = "$(signals_reduction_name) $(rolling_type)")
        TriCorrApplications.plot_contribution!(ax, plot_eeg, signal_times, reduced_signal)
    end
    fig[3,1:2] = ax_rev = Axis(fig; ylabel = "# reviewers", xlabel = "time")
    TriCorrApplications.plot_reviewer_consensus!(ax_rev, plot_eeg)

    roc_column = fig[:,end+1]

    map(enumerate(["μ", "σ"])) do (i, rolling_type)
        rolling_reduction_params = signal_window_fns[rolling_type]
        roc_data = calculate_ROC(signals, signal_times, truth_bounds, reduce_signals_fn; rolling_reduction_params..., rolling_window_s=rolling_window_s, alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, n_θs=100, signals_reduction_params...)
        @info "done. Now plotting."
        
        if !isnothing(roc_data) && any(roc_data.gt_negative .!= 0) && any(roc_data.gt_positive .!= 0)
            roc_plt = data(roc_data) * mapping((:false_positives,:gt_negative) => ((f, gt) -> f / non_seizure_hours) => "FP/Hour", (:true_positives,:gt_positive)=> ((t, gt) -> t / gt) => "TPR") * visual(Lines, color=:blue, linewidth=5)
            roc_drw = draw!(roc_column[i,1], roc_plt, axis=(title="$(title) (Rolling $(rolling_type) of $(signals_reduction_name))", limits=((0.,maximum(roc_data.gt_negative)/non_seizure_hours),(0.,1.))))
            @show "FP/Hr max = $(maximum(roc_data.gt_negative))/$(non_seizure_hours) = $(maximum(roc_data.gt_negative)/non_seizure_hours)"

            roc_drw
        else
            @warn "Cannot plot valid ROC curve: $(title)"
        end
    end
    return fig
end
