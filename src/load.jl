using MAT, NamedDims, EDF

function load_EEG_snippet(path, key)
    test_data_vars = matread(path)

    expected_vars = Set([key])
    @assert Set(keys(test_data_vars)) == expected_vars "$expected_vars ≠ $(keys(test_data_vars))"

    eeg_snippet_arr = test_data_vars[key]
    return NamedDimsArray{(:channel,:time)}(eeg_snippet_arr)
end

# From TimeAxes
is_time(sym::Symbol) = sym == :time

function download_helsinki_eegs(eeg_nums; target_dir=datadir("exp_raw", "helsinki"))
    mkpath(target_dir)
    for i ∈ eeg_nums
        download("https://zenodo.org/record/2547147/files/eeg$(i).edf?download=1", joinpath(target_dir, "eeg$(i).edf"))
    end
end


struct ProcessedEEGv1{T,SIG<:NamedDimsArray{(:channel,:time),T}}
    signals::SIG
    labels::Vector{String}
    sample_rate::Float64
    duration::Float64
end


function normalize_01!(arr)
    arr .-= minimum(arr)
    arr ./= maximum(arr)
    return arr
end
function process_signal(signal::EDF.Signal, sample_rate; f_low=0.1, f_high=70.)
    int_data = signal.samples
    data = int_data .- mean(int_data)
    notch = digitalfilter(Bandstop(50-5, 50+5; fs=sample_rate), Butterworth(6))
    pass = digitalfilter(Bandpass(f_low, f_high; fs=sample_rate), Butterworth(2))
    data = filtfilt(notch, filtfilt(pass, data))' # want row
    data ./= std(data)
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
    labels = [replace(replace(replace(sig.header.label, "-Ref" => ""), "-REF"=>""), "EEG " => "") for sig in edf.signals
            if !any(contains(sig.header.label, ex) for ex in exclude)
        ]
    ProcessedEEGv1(signals, labels, sample_rate, duration)
end

function load_binary_annotations(eeg_num; filepath=datadir("annotations_2017.mat"))
    sum(matread(filepath)["annotat_new"][eeg_num]; dims=1) .== 3
end

function calc_seizure_bounds(annotations::AbstractVector)
    changes = diff(annotations[:])
    onsets = findall(changes .== 1) .+ 1
    offsets = findall(changes .== -1) .+ 1
    if length(onsets) == length(offsets) && onsets[1] < offsets[1]
        return (onsets, offsets)
    else
        error("Must implement where artifact within seizure.")
    end
end

function load_seizure_index_annotations(eeg_num; sample_rate, kwargs...)
    second_annotations = load_binary_annotations(eeg_num; kwargs...)
    calc_seizure_bounds(second_annotations)
end