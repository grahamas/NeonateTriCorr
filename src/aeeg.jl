using DSP

function calculate_aeeg(eeg::AbstractEEG, signal::AbstractVector{T}, fs, upper_freq=0.31, window_len_s=15) where T
    if any(ismissing.(signal))
        @error "aEEG given incomplete signal; Missings corrupt."
    end
    f = digitalfilter(Lowpass(upper_freq, fs=fs), Butterworth(5))
    output = filt(f, abs.(signal))
    output = set_artifacts_missing(output, eeg)
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
        calculate_aeeg(eeg, signals[i_channel,:], fs, upper_freq, window_len_s)
    end
end

function aeeg_lower_margin(channels_margins::Vector)
    mapreduce(hcat, channels_margins) do margins
        margins[1,:]
    end
end