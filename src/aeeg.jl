using DSP

function calculate_aeeg(signal::Vector{T}, fs, upper_freq=0.31, window_len_s=15) where T
    f = digitalfilter(Lowpass(upper_freq, fs=fs), Butterworth(5))
    output = copy(signal)
    next_missing = findfirst(ismissing.(signal))
    prev_missing = 0
    if next_missing == 1
        prev_missing = findfirst(.!ismissing.(signal))
        if prev_missing !== nothing
            next_missing = findfirst(ismissing.(signal[prev_missing:end]))
            prev_missing -= 1
            output[1:prev_missing] .= missing
        else
            next_missing = nothing
        end
    end
    while next_missing !== nothing
        next_missing += prev_missing
        idx = (prev_missing+1):(next_missing-1)
        rect::Vector{Float64} = Float64.(abs.(signal[idx]))
        output[idx] .= filtfilt(f, rect)
        prev_missing = findfirst(.!ismissing.(signal[next_missing:end]))
        if prev_missing === nothing
            break
        else
            prev_missing += next_missing-1
            output[next_missing:prev_missing] .= missing
        end
        next_missing = findfirst(ismissing.(signal[prev_missing+1:end]))
    end
    if prev_missing !== nothing
        rect2::Vector{Float64} = abs.(signal[prev_missing+1:end])
        output[prev_missing+1:end] .= filtfilt(f, rect2)
    else
        output[next_missing:end] .= missing
    end

    window_len_idx = floor(Int,window_len_s*fs)
    window_starts = (0:window_len_idx:(length(output)-window_len_idx)) .+ 1
    perc_9 = floor(Int, 0.09 * window_len_idx)
    perc_93 = floor(Int, 0.93 * window_len_idx)
    margins = mapreduce(hcat, window_starts) do start
        window = output[start:(start+window_len_idx-1)]
        if any(ismissing.(window))
            return [missing, missing]
        else
            sort!(window)
            return [window[perc_9], window[perc_93]]
        end
    end
    return margins
end

function calculate_aeeg(eeg::AbstractEEG; upper_freq=0.31, window_len_s=15)
    signals = get_signal(eeg)
    fs = eeg.sample_rate
    map(1:size(signals,1)) do i_channel
        calculate_aeeg(signals[i_channel,:], fs, upper_freq, window_len_s)
    end
end

function aeeg_lower_margin(channels_margins::Vector)
    mapreduce(hcat, channels_margins) do margins
        margins[1,:]
    end
end