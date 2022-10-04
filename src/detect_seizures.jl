
not_nothing(::Nothing, i) = i
not_nothing(x::Any, _) = x


function coarse_grain_binary_signal(signal::AbstractVector{<:Union{Bool,Missing}}, super_bin_len)
    coarse_grained_signal_len = ceil(Int, length(signal) / super_bin_len)
    coarse_grained_signal = zeros(Union{Bool,Missing}, coarse_grained_signal_len)
    for (i_sig, i_sup) ∈ zip(firstindex(signal):super_bin_len:lastindex(signal), 1:coarse_grained_signal_len)
        fine_grained = signal[i_sig:min(i_sig+super_bin_len-1,lastindex(signal))]
        if any(ismissing.(fine_grained))
            coarse_grained_signal[i_sup] = missing
        else
            coarse_grained_signal[i_sup] = any(fine_grained)
        end
    end
    coarse_grained_signal
end

function coarse_grain_signal(signal::AbstractVector{<:Union{T,Missing}}, super_bin_len; reduction_fn=maximum) where T
    coarse_grained_signal_len = ceil(Int, length(signal) / super_bin_len)
    coarse_grained_signal = zeros(Union{T,Missing}, coarse_grained_signal_len)
    for (i_sig, i_sup) ∈ zip(firstindex(signal):super_bin_len:lastindex(signal), 1:coarse_grained_signal_len)
        fine_grained = signal[i_sig:min(i_sig+super_bin_len-1,lastindex(signal))]
        if any(ismissing.(fine_grained))
            coarse_grained_signal[i_sup] = missing
        else
            coarse_grained_signal[i_sup] = reduction_fn(fine_grained)
        end
    end
    coarse_grained_signal
end

function find_containing_bins(start, stop, bin_start_times)
    # find the bins containing the start:stop interval
    start_idx = not_nothing(findlast(bin_start_times .<= start), firstindex(bin_start_times))
    if start_idx == lastindex(bin_start_times)
        return Int[]
    end
    stop_idx = not_nothing(findfirst(bin_start_times .>= stop), lastindex(bin_start_times)+1)-1
    @assert stop_idx >= start_idx
    return start_idx:stop_idx
end

function calculate_negepoch_non_seizure_hours(bounds, times, epoch_s, snippet_s)
    epoch = epoch_s ÷ snippet_s
    coarse_times = times[begin:epoch:end]
    seizure_epochs = mapreduce(+, bounds, init=0) do (on, off)
        bin_idxs = find_containing_bins(on, off, coarse_times)
        length(bin_idxs)
    end
    return (length(coarse_times) - seizure_epochs) * epoch_s / (60 * 60)
end

function epochize_bounds(truth_bounds, start, stop; epoch_s, unused...)
    truth_bounds_plus_one_epoch = map(truth_bounds) do (on, off)
        (max(on-epoch_s, start), min(off+epoch_s, stop))
    end
    epoch_truth_bounds = discretize_and_merge_bounds(truth_bounds_plus_one_epoch, epoch_s; min_bound=start, max_bound=stop)
    return epoch_truth_bounds
end

# signal, signal_times, bounds -> (target, non_target)
function pospatient_negepoch(signal::AbstractVector{<:Union{Missing,T}}, signal_times, target_bounds::Vector{<:Tuple}; epoch_s, snippets_duration_s) where T
    # One positive per seizure, one negative per non-seizure epoch
    epoch_bin_len = epoch_s ÷ snippets_duration_s
    @assert epoch_bin_len * snippets_duration_s == epoch_s
    epoch_signal = coarse_grain_signal(signal, epoch_bin_len, reduction_fn=maximum)
    epoch_signal_times = signal_times[begin:epoch_bin_len:end]
    epoch_target_bounds = epochize_bounds(target_bounds, first(epoch_signal_times), last(epoch_signal_times); epoch_s=epoch_s)
    epoch_nontarget_bounds = invert_bounds(epoch_target_bounds, epoch_signal_times[begin], epoch_signal_times[end])

    nonmissing_seizure_targets = skipmissing(map(epoch_target_bounds) do (on, off)
        interval_bin_idxs = find_containing_bins(on, off, epoch_signal_times)
        nonmissings = skipmissing(epoch_signal[interval_bin_idxs])
        if isempty(nonmissings)
            missing
        else
            maximum(nonmissings)
        end
    end)
    targets = if isempty(nonmissing_seizure_targets)
        T[]
    else
        maximum(nonmissing_seizure_targets)
    end
    nonmissing_non_targets = map(epoch_nontarget_bounds) do (on, off)
        interval_bin_idxs = find_containing_bins(on, off, epoch_signal_times)
        collect(skipmissing(epoch_signal[interval_bin_idxs]))
    end
    non_targets = if isempty(nonmissing_non_targets)
        T[]
    else
        reduce(vcat, nonmissing_non_targets)
    end

    return (targets, non_targets)
end

function posseizure_negepoch(signal::AbstractVector{<:Union{Missing,T}}, signal_times, target_bounds::Vector{<:Tuple}; epoch_s, snippets_duration_s) where T
    # One positive per seizure, one negative per non-seizure epoch
    epoch_bin_len = epoch_s ÷ snippets_duration_s
    @assert epoch_bin_len * snippets_duration_s == epoch_s
    epoch_signal = coarse_grain_signal(signal, epoch_bin_len, reduction_fn=maximum)
    epoch_signal_times = signal_times[begin:epoch_bin_len:end]
    epoch_target_bounds = epochize_bounds(target_bounds, first(epoch_signal_times), last(epoch_signal_times); epoch_s=epoch_s)
    epoch_nontarget_bounds = invert_bounds(epoch_target_bounds, epoch_signal_times[begin], epoch_signal_times[end])

    targets = collect(skipmissing(map(epoch_target_bounds) do (on, off)
        interval_bin_idxs = find_containing_bins(on, off, epoch_signal_times)
        nonmissings = skipmissing(epoch_signal[interval_bin_idxs])
        if isempty(nonmissings)
            missing
        else
            maximum(nonmissings)
        end
    end))
    nonmissing_non_targets = map(epoch_nontarget_bounds) do (on, off)
        interval_bin_idxs = find_containing_bins(on, off, epoch_signal_times)
        collect(skipmissing(epoch_signal[interval_bin_idxs]))
    end
    non_targets = if isempty(nonmissing_non_targets)
        T[]
    else
        reduce(vcat, nonmissing_non_targets)
    end

    return (targets, non_targets)
end

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
    extrema(skipmissing(arr)) .+ (-sqrt(eps()), sqrt(eps()))
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

function calculate_ROC(signal::AbstractVector, signal_times, truth_bounds; epoch_s, snippets_duration_s, n_θs=100, unused_args...)
    min_θ, max_θ = calculate_threshold_range(signal)
    if !isfinite(min_θ)
        return nothing
    end

    function full_detection_fn(θ)
        evaluate_detection_posseizure_negepoch(truth_bounds, signal .>= θ, signal_times; snippets_duration_s=snippets_duration_s, epoch_s=epoch_s)
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

function standardize_signals!(signals::AbstractVector{<:AbstractArray}; standardization, unused_params...)
    if standardization == "within"
        for idx ∈ eachindex(signals)
            not_missings = .!ismissing.(signals[idx][1,:])
            signals[idx] .-= mean(signals[idx][:, not_missings], dims=:time)
            signals[idx] ./= std(signals[idx][:, not_missings], dims=:time)
        end
    elseif standardization == "across"
        cat_signals = cat(signals..., dims=:time)
        cat_not_missings = .!ismissing.(cat_signals[1,:])
        signal_means = mean(cat_signals[:,cat_not_missings], dims=:time)
        signal_stds = std(cat_signals[:,cat_not_missings], dims=:time)

        for idx ∈ eachindex(signals)
            signals[idx] .-= signal_means
            signals[idx] ./= signal_stds
        end
    else
        error("What's \"$standardization\" standardization?")
    end
end

function standardize_signals!(signals::AbstractArray{<:Union{Missing,Number}}; standardization, unused_params...)
    if standardization == "within"
        not_missings = .!ismissing.(signals[1,:])
        signals .-= mean(signals[:, not_missings], dims=:time)
        signals ./= std(signals[:, not_missings], dims=:time)
    elseif standardization == "across"
        error("Cannot standardize across single patient.")
    else
        error("What's \"$standardization\" standardization?")
    end
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
    return vec(mean(signals, dims=1))
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

# function calc_non_seizure_hours(eeg, bounds)
#     (eeg.duration + mapreduce((x) -> x[1] - x[2], +, bounds, init=0)) / (60 * 60)
# end

function plot_μ_and_σ_signals_and_roc!(fig, signals, signal_times, truth_bounds; analysis_eeg, plot_eeg=analysis_eeg, snippets_duration_s, epoch_s::Number, title, signals_reduction_name, standard_mean=nothing, standard_std=nothing, unused_params...)
    reduce_signals_fn = get_reduce_signals_fn(signals_reduction_name)

    if !isempty(unused_params)
        @warn "Unused params: $unused_params"
    end
    signal_window_fns = Dict(
        "σ" => (signal_sym=:Δσ, window_fn=std),
        "μ" => (signal_sym=:Δμ, window_fn=mean)
    )
    map(enumerate(["μ", "σ"])) do (i, rolling_type)
        rolling_reduction_params = signal_window_fns[rolling_type]
        rolled_signals = roll_signals(signals; snippets_duration_s=snippets_duration_s, unused_params..., rolling_reduction_params...)
        standardized_signals = if isnothing(standard_mean) && isnothing(standard_std)
            not_missings = .!ismissing.(rolled_signals[1,:])
            (rolled_signals .- mean(rolled_signals[:,not_missings], dims=:time)) ./ std(rolled_signals[:,not_missings], dims=:time)
        elseif !isnothing(standard_mean) && !isnothing(standard_std)
            (rolled_signals .- standard_mean) ./ standard_std
        else
            error("Need both or neither standard_mean and standard_std")
        end
        
        reduced_signal = reduce_signals_fn(standardized_signals)

        fig[i,1:2] = ax = Axis(fig; ylabel = "$(signals_reduction_name) $(rolling_type)")
        TriCorrApplications.plot_contribution!(ax, plot_eeg, signal_times, reduced_signal; epoch_s=epoch_s)
    end
    fig[3,1:2] = ax_rev = Axis(fig; ylabel = "# reviewers", xlabel = "time")
    TriCorrApplications.plot_reviewer_consensus!(ax_rev, plot_eeg)

    roc_column = fig[:,end+1]

    map(enumerate(["μ", "σ"])) do (i, rolling_type)
        rolling_reduction_params = signal_window_fns[rolling_type]
        rolled_signals = roll_signals(signals; snippets_duration_s=snippets_duration_s, unused_params..., rolling_reduction_params...)
        standardized_signals = if isnothing(standard_mean) && isnothing(standard_std)
            not_missings = .!ismissing.(rolled_signals[1,:])
            (rolled_signals .- mean(rolled_signals[:,not_missings], dims=:time)) ./ std(rolled_signals[:,not_missings], dims=:time)
        elseif !isnothing(standard_mean) && !isnothing(standard_std)
            (rolled_signals .- standard_mean) ./ standard_std
        else
            error("Need both or neither standard_mean and standard_std")
        end
        
        reduced_signal = reduce_signals_fn(standardized_signals)

        roc_data = calculate_ROC(reduced_signal, signal_times, truth_bounds; 
            epoch_s=epoch_s, 
            snippets_duration_s=snippets_duration_s, 
            n_θs=100, 
            rolling_reduction_params...
        )
        
        if !isnothing(roc_data) && any(roc_data.gt_negative .!= 0) && any(roc_data.gt_positive .!= 0)
            plot_seizure_detection_ROC!(roc_column[i,1], roc_data;
                epoch_s = epoch_s,
                title="$(title) (Rolling $(rolling_type) of $(signals_reduction_name))"
            )
        else
            @warn "Cannot plot valid ROC curve: $(title)"
        end
    end
    return fig
end

function plot_seizure_detection_ROC!(layout, roc_df::DataFrame; epoch_s=nothing, title, roc_plt = plot_NODRAW_seizure_detection_ROC!(roc_df; epoch_s=epoch_s), unused_params...)
    roc_drw = draw!(layout, roc_plt, axis=(title=title, 
    limits=((0.,60.),(0.,1.))
    ))

    roc_drw
end


function plot_NODRAW_seizure_detection_ROC!(roc_df::DataFrame; epoch_s, color=:blue)
    data(roc_df) * mapping(
        (:false_positives,:gt_negative) => 
            ((f, gt) -> f / (gt * epoch_s / (60 * 60))) => 
            "FP/Hour", 
        (:true_positives,:gt_positive) => 
            ((t, gt) -> t / gt) => 
            "TPR"
    ) * visual(Lines, color=color, linewidth=5)
end

function plot_NODRAW_seizure_detection_ROC_standard!(roc_df::DataFrame; color=:blue)
    data(roc_df) * mapping(
        (:false_positives,:gt_negative) => 
            ((f, gt) -> f / gt) => 
            "FPR", 
        (:true_positives,:gt_positive) => 
            ((t, gt) -> t / gt) => 
            "TPR"
    ) * visual(Lines, color=color, linewidth=5)
end


