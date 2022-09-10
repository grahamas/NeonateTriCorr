using DSP

function calculate_aEEG(eeg::AbstractEEG, signal::AbstractVector{T}, fs, lowpass_freq, snippets_duration_s, lower_margin_perc, upper_margin_perc) where T
    if any(ismissing.(signal))
        @error "aEEG given incomplete signal; Missings corrupt."
    end
    f = digitalfilter(Lowpass(lowpass_freq, fs=fs), Butterworth(5))
    output = filt(f, abs.(signal))
    output = set_artifacts_missing(output, eeg)
    window_len_idx = floor(Int,snippets_duration_s*fs)
    window_starts = (0:window_len_idx:(length(output)-window_len_idx)) .+ 1
    perc_9 = floor(Int, lower_margin_perc * window_len_idx)
    perc_93 = floor(Int, upper_margin_perc * window_len_idx)
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

function calculate_aEEG(eeg::AbstractEEG; lowpass_freq, snippets_duration_s, lower_margin_perc, upper_margin_perc, unused_params...)
    signals = get_signal(eeg)
    fs = eeg.sample_rate
    map(1:size(signals,1)) do i_channel
        calculate_aEEG(eeg, signals[i_channel,:], fs, lowpass_freq, snippets_duration_s, lower_margin_perc, upper_margin_perc)
    end
end

function aEEG_lower_margin(channels_margins::Vector)
    NamedDimsArray{(:channel,:time)}(mapreduce(vcat, channels_margins) do margins
        margins[1,:]'
    end)
end