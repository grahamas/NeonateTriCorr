
not_nothing(::Nothing, i) = i
not_nothing(x::Any, _) = x


function coarse_grain_binary_signal(signal::AbstractVector{<:Union{Bool,Missing}}, super_bin_len)
    coarse_grained_signal_len = ceil(Int, length(signal) / super_bin_len)
    coarse_grained_signal = zeros(Union{Bool,Missing}, coarse_grained_signal_len)
    for (i_sig, i_sup) ∈ zip(firstindex(signal):super_bin_len:lastindex(signal), 1:coarse_grained_signal_len)
        coarse_grained_signal[i_sup] = any(signal[i_sig:min(i_sig+super_bin_len-1,lastindex(signal))])
    end
    coarse_grained_signal
end
function find_containing_bins(start, stop, bin_start_times)
    # find the bins containing the start:stop interval
    start_idx = not_nothing(findlast(bin_start_times .<= start), firstindex(bin_start_times))
    if start_idx == lastindex(bin_start_times)
        return Int[]
    end
    stop_idx = not_nothing(findfirst(bin_start_times .>= stop)-1, lastindex(bin_start_times))
    @assert stop_idx >= start_idx
    return start_idx:stop_idx
end
function evaluate_detection_posseizure_negepoch(truth_bounds::Vector{<:Tuple}, detections::Vector{T}, detections_times; epoch_bin_len) where T
    # One positive per seizure, one negative per non-seizure epoch
    # 
    detections = coarse_grain_binary_signal(detections, epoch_bin_len)
    detections_times = detections_times[begin:epoch_bin_len:end]
    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    for (on, off) ∈ truth_bounds
        interval_bin_idxs = find_containing_bins(on, off, detections_times)
        interval_bins = @view detections[interval_bin_idxs]
        if all(ismissing.(interval_bins))
            continue
        end
        detected_seizure = any(skipmissing(interval_bins))
        true_positives += detected_seizure
        false_negatives += !detected_seizure
        interval_bins .= missing # should exclude adjacent bins
    end
    true_negatives = count(skipmissing(detections) .== false)
    false_positives = count(skipmissing(detections) .== true)
    return (
        gt_positive = length(truth_bounds),
        gt_negative = count(.!ismissing.(detections)),
        true_positives = true_positives,
        false_positives = false_positives,
        true_negatives = true_negatives,
        false_negatives = false_negatives
    )
end

function add_grace_to_truth_bounds(truth_bounds, times, grace_s)
    if isempty(truth_bounds)
        return eltype(truth_bounds)[]
    end
    min_time = Int(minimum(times))
    max_time = Int(maximum(times))
    truth_bounds_with_grace_potentially_overlapping = [(max(min_time,on-grace_s),min(max_time,off+grace_s)) for (on,off) ∈ truth_bounds]
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

# function evaluate_detection_posseizure_negalerts(truth_bounds::AbstractVector{<:Tuple}, detections::AbstractVector, detections_times; snippets_duration_s, alert_grace_s=60)
#     # One positive per seizure, one negative per alert period
#     alert_grace = alert_grace_s ÷ snippets_duration_s
#     detections = Vector{Union{Int,Missing}}(detections)
#     len = length(detections)
#     true_positives = 0
#     false_positives = 0
#     true_negatives = 0
#     false_negatives = 0
#     valid_truth_bounds = [(on,off) for (on,off) ∈ truth_bounds if !all(ismissing.(detections[on .<= detections_times .< off]))]
#     truth_bounds_with_grace = add_grace_to_truth_bounds(valid_truth_bounds, detections_times, alert_grace_s)
#     gt_positive = length(truth_bounds_with_grace)
#     for (on, off) ∈ truth_bounds_with_grace # should combine adjacent seizures
#         detected_seizure = any(skipmissing(detections[on .<= detections_times .< off]) .== 1)
#         true_positives += detected_seizure
#         false_negatives += !detected_seizure
#         detections[on .<= detections_times .< off] .= missing # should exclude adjacent bins
#     end
#     for i ∈ eachindex(detections)
#         detect = detections[i]
#         if ismissing(detect) || (detect == 2) || (detect < 0)
#             continue
#         elseif detect == 0
#             # new grace period; correct negative detection
#             true_negatives += 1
#             detections[min(i+1,len):min(i+alert_grace-1,len)] .+= 2
#         elseif (detect == 1)
#             # new grace period; incorrect positive detection
#             false_positives += 1
#             detections[min(i+1,len):min(i+alert_grace-1,len)] .-= 10
#         elseif (detect == 3)
#             # current grace period previously negative, now positive
#             false_positives += 1
#             true_negatives -= 1
#             detections[min(i+1,len):min(i+alert_grace-1,len)] .-= 10
#         end
#     end
#     return (
#         gt_positive = gt_positive,
#         gt_negative = true_negatives + false_positives,
#         true_positives = true_positives,
#         false_positives = false_positives,
#         true_negatives = true_negatives,
#         false_negatives = false_negatives
#     )
# end

function apply_rolling_deviation_window(signals::NamedDimsArray{L,T}, window_fn, window_len) where {L, T}
    n_signals = size(signals, 1)
    window_red = mapreduce(hcat, 1:n_signals) do signal_num
        windowed = [window_fn(signals[signal_num,(i-window_len+1):(i)]) for i ∈ window_len:(size(signals,2))]
        (windowed .- mean(skipmissing(windowed))) ./ std(skipmissing(windowed))
    end
    NamedDimsArray{L}(hcat(zeros(Union{eltype(window_red),Missing}, n_signals, window_len - 1) .+ missing, window_red'))
end

function apply_rolling_window(signals::NamedDimsArray{L,T}, window_fn, window_len) where {L, T}
    n_signals = size(signals, 1)
    window_red = mapreduce(hcat, 1:n_signals) do signal_num
        [window_fn(signals[signal_num,(i-window_len+1):(i)]) for i ∈ window_len:(size(signals,2))]
    end
    NamedDimsArray{L}(hcat(zeros(Union{eltype(window_red),Missing}, n_signals, window_len - 1) .+ missing, window_red'))
end

function detect_deviation_window(signals, window_fn, window_len, motifs_criterion)
    rolling_deviation = apply_rolling_deviation_window(signals, window_fn, window_len)
    detections = map(motifs_criterion, eachslice(rolling_deviation, dims=1))
    return detections
end

function motif_weights_based_on_pvalues(df, effect_sym, n_nonzero, n_signals)
    @assert !isempty(df)
    descending_significance = sort(zip(df.p, df[:, effect_sym], 1:n_signals) |> collect)
    weights = zeros(n_signals)
    weights[descending_significance[1:n_nonzero] .|> x -> x[3]] .=(descending_significance[1:n_nonzero] .|> x -> sign(x[2]))
    return weights / norm(weights)
end

function motif0_weight_based_effect(df, effect_sym, n_nonzero, n_signals)
    @assert !isempty(df)
    weights = zeros(n_signals)
    weights[1] = sign(df[1, effect_sym])
    return weights
end

function calculate_threshold_range(arr)
    extrema(skipmissing(arr))
end

function calculate_threshold_range(arrs::AbstractVector{<:AbstractArray})
    (minimum(minimum.(skipmissing.(arrs))) - sqrt(eps()), maximum(maximum.(skipmissing.(arrs))) + sqrt(eps()))
end

function calculate_ROC(full_detection_fn::Function, θs)
    rocnums = map(θs) do θ
        rocnum = full_detection_fn(θ)
        merge((θ=θ,), rocnum)
    end
    DataFrame(rocnums)
end

function calculate_ROC(raw_signals::AbstractArray, signal_times, truth_bounds, reduce_signals_fn::Function; alert_grace_s, snippets_duration_s, n_θs=100, processing_args...)
    reduced_signal = reduce_signals_fn(raw_signals)
    min_θ, max_θ = calculate_threshold_range(reduced_signal)
    if !isfinite(min_θ)
        return nothing
    end

    function full_detection_fn(θ)
        evaluate_detection_posseizure_negalerts(truth_bounds, reduced_signal .>= θ, signal_times; snippets_duration_s=snippets_duration_s, alert_grace_s=alert_grace_s)
    end

    calculate_ROC(full_detection_fn, range(min_θ, max_θ, length=n_θs))
end

function calculate_AUC(roc_df::DataFrame)
    roc_df = sort(roc_df, :θ, rev=true)
    TPRs = roc_df.true_positives ./ roc_df.gt_positive
    FPRs = roc_df.false_positives ./ roc_df.gt_negative
    prev_fpr = 0.
    prev_tpr = 0.
    AUC = 0.
    for (fpr, tpr) ∈ zip(FPRs, TPRs)
        rect_area = (prev_tpr) * (fpr - prev_fpr)
        triangle = 0.5 * (tpr - prev_tpr) * (fpr - prev_fpr) # when TPR and FPR change simulatenously
        AUC += rect_area + triangle
        prev_fpr = fpr
        prev_tpr = tpr
    end
    return AUC
end

function roll_signals(raw_signals::NamedDimsArray; window_fn, rolling_window_s, snippets_duration_s, unused_params...)
    raw_signals = NamedDimsArray{(:_, :time)}(raw_signals)
    window_len = rolling_window_s ÷ snippets_duration_s
    NamedDimsArray{(:_, :time)}(apply_rolling_window(raw_signals, window_fn, window_len))
end

function reduce_signal_distance(raw_signals::NamedDimsArray; window_fn, rolling_window_s, snippets_duration_s, unused_params...)
    raw_signals = NamedDimsArray{(:_, :time)}(raw_signals)
    window_len = rolling_window_s ÷ snippets_duration_s
    rolled_signal = apply_rolling_window(raw_signals, window_fn, window_len)
    count_nonmissing_times = count(.!ismissing.(sum(rolled_signal, dims=1)))
    unmissing_rolled_signal = copy(rolled_signal)
    unmissing_rolled_signal[ismissing.(rolled_signal)] .= 0.
    mean_signal = sum(unmissing_rolled_signal, dims=2) ./ count_nonmissing_times
    return vec(sqrt.(sum((rolled_signal .- mean_signal) .^ 2, dims=1)))
end

function reduce_signal_meanall(signals::NamedDimsArray)
    return vec(mean(rolled_signal, dims=1))
end

function reduce_signal_meanallabs(signals::NamedDimsArray)
    return vec(mean(abs.(signals), dims=1))
end

function reduce_signal_maxany(signals::NamedDimsArray)
    return vec(maximum(signals, dims=1))
end

function reduce_signal_maxanyabs(signals::NamedDimsArray)
    return vec(maximum(abs.(signals), dims=1))
end


function get_reduce_signals_fn(reduction_type)
    Dict(
        "maxany" => reduce_signal_maxany,
        "meanall" => reduce_signal_meanall,
        "maxanyabs" => reduce_signal_maxanyabs,
        "meanallabs" => reduce_signal_meanallabs,
        "distance" => reduce_signal_distance
    )[reduction_type]
end


function plot_μ_and_σ_signals_and_roc(args...; resolution=(1300, 1000), kwargs...)
    fig = Figure(resolution=resolution)
    plot_μ_and_σ_signals_and_roc!(fig, args...; kwargs...)
    return fig
end

function calc_non_seizure_hours(eeg, bounds)
    (eeg.duration + mapreduce((x) -> x[1] - x[2], +, bounds, init=0)) / (60 * 60)
end

function plot_μ_and_σ_signals_and_roc!(fig, signals, signal_times, truth_bounds; analysis_eeg, plot_eeg=analysis_eeg, snippets_duration_s, alert_grace_s, title, signals_reduction_name, unused_params...)
    reduce_signals_fn = get_reduce_signals_fn(signals_reduction_name)

    if !isempty(unused_params)
        @warn "Unused params: $unused_params"
    end
    
    seizures_with_grace_period = add_grace_to_truth_bounds(analysis_eeg.seizure_annotations, get_times(analysis_eeg, sample_rate=1/snippets_duration_s), alert_grace_s)
    seizure_and_artifact_bounds = merge_bounds(seizures_with_grace_period, analysis_eeg.artifact_annotations)
    non_seizure_hours = calc_non_seizure_hours(analysis_eeg, seizure_and_artifact_bounds)

    @info "Calculating ROC curve..."
    signal_window_fns = Dict(
        "σ" => (signal_sym=:Δσ, window_fn=std),
        "μ" => (signal_sym=:Δμ, window_fn=mean)
    )
    map(enumerate(["μ", "σ"])) do (i, rolling_type)
        rolling_reduction_params = signal_window_fns[rolling_type]
        reduced_signal = reduce_signals_fn(signals)

        fig[i,1:2] = ax = Axis(fig; ylabel = "$(signals_reduction_name) $(rolling_type)")
        TriCorrApplications.plot_contribution!(ax, plot_eeg, signal_times, reduced_signal)
    end
    fig[3,1:2] = ax_rev = Axis(fig; ylabel = "# reviewers", xlabel = "time")
    TriCorrApplications.plot_reviewer_consensus!(ax_rev, plot_eeg)

    roc_column = fig[:,end+1]

    map(enumerate(["μ", "σ"])) do (i, rolling_type)
        rolling_reduction_params = signal_window_fns[rolling_type]
        roc_data = calculate_ROC(signals, signal_times, truth_bounds, 
            reduce_signals_fn; 
            alert_grace_s=alert_grace_s, 
            snippets_duration_s=snippets_duration_s, 
            n_θs=100, 
            rolling_reduction_params...
        )
        @info "done. Now plotting."
        
        if !isnothing(roc_data) && any(roc_data.gt_negative .!= 0) && any(roc_data.gt_positive .!= 0)
            plot_seizure_detection_ROC!(roc_column[i,1], roc_data;
                non_seizure_hours=non_seizure_hours, 
                title="$(title) (Rolling $(rolling_type) of $(signals_reduction_name))"
            )
        else
            @warn "Cannot plot valid ROC curve: $(title)"
        end
    end
    return fig
end

function plot_seizure_detection_ROC!(layout, roc_df::DataFrame; non_seizure_hours, title, roc_plt = plot_NODRAW_seizure_detection_ROC!(roc_df, non_seizure_hours), unused_params...)
    roc_drw = draw!(layout, roc_plt, axis=(title=title, 
    #limits=((0.,maximum(roc_df.gt_negative)/non_seizure_hours),(0.,1.))
    ))

    roc_drw
end


function plot_NODRAW_seizure_detection_ROC!(roc_df::DataFrame, non_seizure_hours; color=:blue)
    data(roc_df) * mapping(
        (:false_positives,:gt_negative) => 
            ((f, gt) -> f / non_seizure_hours) => 
            "FP/Hour", 
        (:true_positives,:gt_positive) => 
            ((t, gt) -> t / gt) => 
            "TPR"
    ) * visual(Lines, color=color, linewidth=5)
end

function plot_NODRAW_seizure_detection_ROC_standard!(roc_df::DataFrame, non_seizure_hours; color=:blue)
    data(roc_df) * mapping(
        (:false_positives,:gt_negative) => 
            ((f, gt) -> f / gt) => 
            "FPR", 
        (:true_positives,:gt_positive) => 
            ((t, gt) -> t / gt) => 
            "TPR"
    ) * visual(Lines, color=color, linewidth=5)
end