using MAT, NamedDims, EDF, Statistics, DSP, Downloads

include(datadir("helsinki_eeg_artifacts.jl"))

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
        eeg_name = "eeg$(i).edf"
        eeg_path = joinpath(target_dir, eeg_name)
        if isfile(eeg_path)
            @warn "File already exists: $(eeg_name)"
            continue
        end
        p = Progress(0; enabled=false, desc="EEG $(i):")
        function update_progress(total, now)
            if p.enabled == false && total > 0
                p.n = total
                p.enabled = true
            end
            ProgressMeter.update!(p, now)
        end
        Downloads.download("https://zenodo.org/record/2547147/files/eeg$(i).edf?download=1", eeg_path; progress=update_progress)
    end
end

abstract type AbstractProcessedEEG end

struct ProcessedEEGv2{T,SIG<:NamedDimsArray{(:channel,:time),T}} <: AbstractProcessedEEG
    signals::SIG
    labels::Vector{String}
    sample_rate::Float64
    duration::Float64
    seizure_annotations::Vector{Tuple{Float64,Float64}}
    artifact_annotations::Vector{Tuple{Float64,Float64}}
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

function ProcessedEEG(edf::EDF.File; exclude=[], seizure_annotations=Tuple{Float64,Float64}[], artifact_annotations=Tuple{Float64,Float64}[])
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
    ProcessedEEGv2(signals, labels, sample_rate, duration, seizure_annotations, artifact_annotations)
end

function load_binary_annotations(eeg_num; filepath=datadir("annotations_2017.mat"))
    (sum(matread(filepath)["annotat_new"][eeg_num]; dims=1) .== 3)[:]
end

function calc_seizure_bounds(annotations::AbstractVector)
    changes = diff(annotations[:])
    onsets = findall(changes .== 1) .+ 1
    offsets = findall(changes .== -1) .+ 1
    if length(onsets) == length(offsets) && onsets[1] < offsets[1]
        return zip(onsets, offsets)
    else
        error("Must implement where artifact within seizure.")
    end
end

function load_seizure_annotations(eeg_num; kwargs...)
    second_annotations = load_binary_annotations(eeg_num; kwargs...)
    calc_seizure_bounds(second_annotations)
end

function load_artifact_annotations(eeg_num)
    helsinki_eeg_artifacts[eeg_num]
end

function load_helsinki_eeg(eeg_num)
    edf = EDF.read(datadir("exp_raw", "helsinki", "eeg$(eeg_num).edf"))
    eeg = ProcessedEEG(edf; exclude=helsinki_eeg_bad_channels[eeg_num], seizure_annotations=Tuple{Float64,Float64}.(load_seizure_annotations(eeg_num)) |> collect, artifact_annotations=load_artifact_annotations(eeg_num))
end