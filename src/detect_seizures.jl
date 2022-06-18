
function evaluate_detection_posseizure_negepoch(seizure_bounds::Vector{<:Tuple}, detections::Vector)
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    for (on, off) ∈ seizure_bounds # should combine adjacent seizures
        detected_seizure = any(skipmissing(detections[on:off]))
        true_positives += detected_seizure
        false_negatives += !detected_seizure
        detections[on:off] .= missing # should exclude adjacent bins
    end
    true_negatives = count(skipmissing(detections) .== false)
    false_positives = count(skipmissing(detections) .== true)
    return (
        gt_positive = length(seizure_bounds),
        gt_negative = count(.!ismissing.(detections)), # should combine seconds
        true_positives = true_positives,
        false_positives = false_positives,
        true_negatives = true_negatives,
        false_negatives = false_negatives
    )
end

function add_grace_to_seizure_bounds(seizure_bounds, grace, min_idx, max_idx)
    if isempty(seizure_bounds)
        return eltype(seizure_bounds)[]
    end
    seizure_bounds_with_grace_potentially_overlapping = [(max(min_idx,on-grace),min(max_idx,off+grace)) for (on,off) ∈ seizure_bounds]
    seizure_bounds_with_grace = eltype(seizure_bounds_with_grace_potentially_overlapping)[]
    prev_bounds = first(seizure_bounds_with_grace_potentially_overlapping)
    for bounds ∈ seizure_bounds_with_grace_potentially_overlapping[begin+1:end]
        if prev_bounds[end] >= bounds[begin]
            merged_bounds = (prev_bounds[begin], bounds[end])
            prev_bounds = merged_bounds
        else
            push!(seizure_bounds_with_grace, prev_bounds)
            prev_bounds = bounds
        end
    end
    return seizure_bounds_with_grace
end

function evaluate_detection_posseizure_negalerts(seizure_bounds::Vector{<:Tuple}, detections::Vector; alert_grace=60)
    detections = Vector{Union{Int,Missing}}(detections)
    len = length(detections)
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    seizure_bounds_with_grace = add_grace_to_seizure_bounds(seizure_bounds, alert_grace, firstindex(detections), lastindex(detections))
    gt_positive = length(seizure_bounds_with_grace)
    for (on, off) ∈ seizure_bounds_with_grace # should combine adjacent seizures
        detected_seizure = any(skipmissing(detections[on:off]) .== 1)
        true_positives += detected_seizure
        false_negatives += !detected_seizure
        detections[on:off] .= missing # should exclude adjacent bins
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

function apply_rolling_deviation_window(arr, window_fn, window_len)
    window_red = mapreduce(hcat, 1:14) do motif_num
        stds = [window_fn(arr[motif_num,(i-window_len+1):(i)]) for i ∈ window_len:(size(arr,:time))]
        (stds .- mean(skipmissing(stds))) ./ std(skipmissing(stds))
    end
    vcat(zeros(eltype(window_red), window_len-1, 14) .+ missing, window_red)
end

function detect_deviation_window(arr, window_fn, window_len, motifs_criterion)
    rolling_deviation = apply_rolling_deviation_window(arr, window_fn, window_len)
    detections = map(motifs_criterion, eachslice(rolling_deviation, dims=1))
end

function motif_weights_based_on_pvalues(df, effect_sym, n_nonzero)
    if isempty(df)
        return ones(14) / 14
    end
    descending_significance = sort(zip(df.p, df[:, effect_sym], 1:14) |> collect)
    weights = zeros(14)
    weights[descending_significance[1:n_nonzero] .|> x -> x[3]] .=(descending_significance[1:n_nonzero] .|> x -> sign(x[2]))
    return weights / norm(weights)
end


function calculate_threshold_range(arr, seizure_bounds; alert_grace)
    min_θ = Inf
    max_θ = -Inf
    seizure_bounds_with_grace = add_grace_to_seizure_bounds(seizure_bounds, alert_grace, firstindex(arr), lastindex(arr))
    for bounds ∈ seizure_bounds_with_grace
        non_artifact = skipmissing(arr[bounds[1]:bounds[2]])
        if isempty(non_artifact)
            return (Inf, -Inf)
        end
        bounds_min, bounds_max = extrema(non_artifact)
        min_θ = min(min_θ, bounds_min)
        max_θ = max(max_θ, bounds_max)
    end
    # below min_θ and above max_θ, no seizures are detected
    return (min_θ, max_θ)
end

function calculate_ROC(full_detection_fn::Function, θs)
    rocnums = map(θs) do θ
        rocnum = full_detection_fn(θ)
        merge((θ=θ,), rocnum)
    end
    DataFrame(rocnums)
end

function calculate_ROC(contributions::AbstractArray, test_results_df::DataFrame, seizure_bounds, signal_sym, window_fn, window_len; alert_grace, n_most_significant_motifs=5, n_θs=100)
    max_significance_weights = motif_weights_based_on_pvalues(test_results_df, signal_sym, n_most_significant_motifs)
    signal = apply_rolling_deviation_window(contributions, window_fn, window_len) * max_significance_weights
    min_θ, max_θ = calculate_threshold_range(signal, seizure_bounds; alert_grace=alert_grace)
    if !isfinite(min_θ)
        return nothing
    end

    function full_detection_fn(θ)
        evaluate_detection_posseizure_negalerts(seizure_bounds, signal .>= θ; alert_grace=alert_grace)
    end

    calculate_ROC(full_detection_fn, range(min_θ, max_θ, length=n_θs))
end

function calculate_patient_ROC(patient::Int, all_motif_results_df, signal_sym, window_fn, window_len; snippets_duration=1, target_match_str, 
    min_reviewers_per_seizure, kwargs...)
    motif_results_df = filter(:patient => p -> p == patient, all_motif_results_df)

    jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    contributions = jld_dict["contributions"]

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient; min_reviewers_per_seizure=min_reviewers_per_seizure)

    calculate_ROC(contributions, motif_results_df, seizure_bounds, signal_sym, window_fn, window_len; kwargs...)
end
