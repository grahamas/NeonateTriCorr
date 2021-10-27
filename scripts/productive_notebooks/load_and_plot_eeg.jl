quickactivate("NeonateTriCorr")

include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics

f_low = 0.1
f_high = 70

function normalize_01!(arr)
    arr .-= minimum(arr)
    arr ./= maximum(arr)
    return arr
end
function process_signal(signal::EDF.Signal, sample_rate)
    int_data = signal.samples
    data = int_data .- mean(int_data)
    notch = digitalfilter(Bandstop(50-5, 50+5; fs=sample_rate), Butterworth(6))
    pass = digitalfilter(Bandpass(f_low, f_high; fs=sample_rate), Butterworth(2))
    data = filtfilt(notch, filtfilt(pass, data))' # want row
end

function ProcessedEEG(edf::EDF.File; exclude=[])
    @assert edf.header.is_contiguous
    n_records = edf.header.record_count
    seconds_per_record = edf.header.seconds_per_record
    samples_per_record = edf.signals[1].header.samples_per_record
    # FIXME should verify sample_rate same for all signals
    sample_rate = samples_per_record / seconds_per_record
    duration = n_records * seconds_per_record
    signals = NamedDimsArray{(:channel,:time)}(
        vcat(
            [process_signal(sig, sample_rate) for sig in edf.signals
             if !any(contains(sig.header.label, ex) for ex in exclude)
            ]...
        )
    )
    labels = [replace(replace(sig.header.label, "-Ref" => ""), "EEG " => "") for sig in edf.signals
            if !any(contains(sig.header.label, ex) for ex in exclude)
        ]
    ProcessedEEGv1(signals, labels, sample_rate, duration)
end

edf = EDF.read(datadir("exp_raw", "helsinki", "eeg50.edf"))
eeg = ProcessedEEG(edf; exclude=["ECG","Resp","Cz"])

draw_eeg_traces(eeg)

trunc_n_seconds = 4
truncated_start = round(Int, 10*eeg.sample_rate)
truncated_len = round(Int, trunc_n_seconds * eeg.sample_rate)
truncated_eeg = ProcessedEEGv1(eeg.signals[:, truncated_start:truncated_start+truncated_len], eeg.labels, eeg.sample_rate, truncated_len / eeg.sample_rate)

draw_eeg_traces(truncated_eeg)