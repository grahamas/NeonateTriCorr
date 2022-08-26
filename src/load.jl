using MAT, EDF, Statistics, DSP, Downloads, CSV, Dates
using DataFrames

include(scriptsdir("helsinki_eeg_bad_channels.jl"))
include(scriptsdir("twente_eeg_bad_channels.jl"))

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

function process_signal(signal::EDF.Signal; seconds_per_record, asserted_samples_per_record, f_low=0.1, f_high=70., mains_hz, mains_bandwidth=10)
    @assert signal.header.samples_per_record == asserted_samples_per_record "Asserted $(asserted_samples_per_record) != $(signal.header.samples_per_record) (channel: $(signal.header.label))"
    sample_rate = asserted_samples_per_record / seconds_per_record
    int_data = signal.samples
    data = int_data .- mean(int_data)
    mains_halfband = mains_bandwidth/2
    notch = digitalfilter(Bandstop(mains_hz-mains_halfband, mains_hz+mains_halfband; fs=sample_rate), Butterworth(6))
    pass = digitalfilter(Bandpass(f_low, f_high; fs=sample_rate), Butterworth(2))
    data = filtfilt(notch, filtfilt(pass, data))' # want row
    data ./= std(data)
end

function ProcessedEEG(edf::EDF.File; 
        exclude_channels=[], 
        seizure_annotations=Tuple{Float64,Float64}[], 
        artifact_annotations=Tuple{Float64,Float64}[], 
        label_replace = identity,
        seizure_reviewers_count = Int[],
        excluded_artifact_grades,
        kwargs...
    )
    @assert edf.header.is_contiguous
    n_records = edf.header.record_count
    seconds_per_record = edf.header.seconds_per_record
    
    edf_signals = [sig for sig ∈ edf.signals
        if (
            sig isa EDF.Signal &&
            !any(contains(lowercase(sig.header.label), lowercase(ex)) for ex in exclude_channels)
        )
    ]
    first_samples_per_record = edf_signals[begin].header.samples_per_record
    # FIXME should verify sample_rate same for all signals
    sample_rate = Int(first_samples_per_record / seconds_per_record)
    duration = n_records * seconds_per_record
    signals = NamedDimsArray{(:channel,:time)}(
        vcat(
            [
                process_signal(sig; 
                    seconds_per_record = seconds_per_record, 
                    asserted_samples_per_record = first_samples_per_record, 
                    kwargs...
                ) for sig in edf_signals
            ]...
        )
    )
    labels = [label_replace(sig.header.label) for sig in edf_signals
        ]
    ProcessedEEGv8(signals, labels, sample_rate, 0, duration, seizure_annotations, artifact_annotations, seizure_reviewers_count, excluded_artifact_grades)
end

function load_count_annotations(eeg_num; filepath=scriptsdir("annotations_2017.mat"))
    map(Int, sum(matread(filepath)["annotat_new"][eeg_num]; dims=1)[:])
end

function load_binary_annotations(eeg_num; filepath=scriptsdir("annotations_2017.mat"))
    (sum(matread(filepath)["annotat_new"][eeg_num]; dims=1) .== 3)[:]
end

function calc_seizure_bounds(annotations::AbstractVector)
    changes = diff(vcat([0], annotations[:], [0]))
    onsets = findall(changes .== 1)
    offsets = findall(changes .== -1)
    if !isempty(offsets) && offsets[end] == length(annotations)+1
        offsets[end] = length(annotations)
    end
    if length(onsets) == length(offsets) && ((length(onsets) == 0) || (onsets[1] < offsets[1]))
        return Vector{Tuple{Int,Int}}(zip(onsets, offsets) |> collect)
    else
        @show annotations
        @show onsets
        @show offsets
        error("Must implement where artifact within seizure.")
    end
end

function load_helsinki_seizure_annotations(eeg_num; min_reviewers_per_seizure=3, kwargs...)
    second_annotations = load_count_annotations(eeg_num; kwargs...)
    consensus_bounds = calc_seizure_bounds(second_annotations .>= min_reviewers_per_seizure)
    return (consensus_bounds, second_annotations)
end

function parse_grade(text)
    text = lowercase(text)
    if occursin("grade i", text) || occursin("grade 1", text)
        return 1
    elseif occursin("grade ii", text) || occursin("grade 2", text)
        return 2
    else
        error("Could not parse grade: $text")
    end
end
function parse_start(text, start_time)
    d1, time = split(text)
    @assert d1 == "d1"
    return Second(Time(time) - start_time).value |> Float64
end
function parse_artifact_duration(text)
    duration, sec = split(text)
    @assert lowercase(sec) == "sec"
    return parse(Int, duration)
end
function artifact_tuple(start, duration, artifact_buffer=1)
    if duration > 0
        return (start-artifact_buffer, start+duration+artifact_buffer)
    else
        return (start-artifact_buffer, start+artifact_buffer)
    end
end
function ordered_overlap(tup1::NTuple{2}, tup2::NTuple{2})
    # Assumed that tup1[1] <= tup2[1]
    if tup2[1] <= tup1[2]
        return true
    else
        return false
    end
end
function collapse_tuples!(tups::Array{T}) where T
    dummy_val = (-1.,-1.)
    sort!(tups)
    prev_extant_i = 1
    for i ∈ eachindex(tups)[begin+1:end]
        if tups[i] != dummy_val && ordered_overlap(tups[prev_extant_i], tups[i])
            tups[prev_extant_i] = (tups[prev_extant_i][1], tups[i][2])
            tups[i] = dummy_val
        else
            prev_extant_i = i
        end
    end
    filter!(!=(dummy_val), tups)
end

function load_helsinki_artifact_annotations(eeg_num, start_time::Time, excluded_grades=(1,))
    df = CSV.read(scriptsdir("helsinki_artifacts.csv"), DataFrame)
    subset!(df, "Patient #" => ByRow(==(eeg_num)))
    output_df = DataFrame(
        grade = parse_grade.(df.Text),
        start = parse_start.(df.Time, Ref(start_time)),
        duration = parse_artifact_duration.(df.Duration)
    )
    possibly_intersecting_tuples = Tuple{Int,Int}[artifact_tuple(t.start, t.duration) for t in eachrow(output_df) if t.grade ∈ excluded_grades]
    collapse_tuples!(possibly_intersecting_tuples)
    return possibly_intersecting_tuples
end

function load_helsinki_eeg(eeg_num::Int; min_reviewers_per_seizure=3, excluded_artifact_grades=[1])
    edf = EDF.read(datadir("exp_raw", "helsinki", "eeg$(eeg_num).edf"))
    seizures_start_stop, seizure_reviewers_count = load_helsinki_seizure_annotations(eeg_num; min_reviewers_per_seizure=min_reviewers_per_seizure)
    ProcessedEEG(edf; 
        exclude_channels=helsinki_eeg_bad_channels[eeg_num], 
        seizure_annotations=seizures_start_stop, 
        artifact_annotations=load_helsinki_artifact_annotations(eeg_num, Time(edf.header.start), excluded_artifact_grades),
        label_replace = (label) -> replace(replace(replace(label, "-Ref" => ""), "-REF"=>""), "EEG " => ""),
        mains_hz=50,
        seizure_reviewers_count = seizure_reviewers_count,
        excluded_artifact_grades = excluded_artifact_grades
    )
end

function get_relevant_annotations(annot::EDF.AnnotationsSignal)
    filter(rec -> rec.annotations != [""], vcat(annot.records...))
end

function get_twente_seizure_annotations(tals::Vector{EDF.TimestampedAnnotationList}; seizure_keywords=["seizure", "epi", "eppi"])
    seizure_bounds = Tuple{Float64,Float64}[]
    map(tals) do tal 
        if any(
                label -> any(
                    occursin.(seizure_keywords, Ref(lowercase(label)))
                ),
                tal.annotations
            )
            onset = tal.onset_in_seconds
            duration = tal.duration_in_seconds
            if tal.duration_in_seconds !== nothing && tal.duration_in_seconds > 0
                push!(seizure_bounds, (onset, onset + duration))
            else
                push!(seizure_bounds, (onset-1, onset+1))
            end
        end
    end
    sort!(seizure_bounds)
    collapse_tuples!(seizure_bounds)
    return seizure_bounds
end

function get_twente_artifact_annotations(tals::Vector{EDF.TimestampedAnnotationList}; artifact_keywords=[])#["ogen", "hyperventilatie"])
    artifact_bounds = Tuple{Float64,Float64}[]
    @show tals
    map(tals) do tal 
        if any(
                label -> any(
                    occursin.(artifact_keywords, Ref(lowercase(label)))
                ), 
                tal.annotations
            )
            onset = tal.onset_in_seconds
            duration = tal.duration_in_seconds
            if tal.duration_in_seconds !== nothing && tal.duration_in_seconds > 0
                push!(artifact_bounds, (onset, onset + duration))
            else
                push!(artifact_bounds, (onset-1, onset+1))
            end
        end
    end
    sort!(artifact_bounds)
    collapse_tuples!(artifact_bounds)
    return artifact_bounds
end

load_twente_eeg(num::Number; kwargs...) = load_twente_eeg("EEG$(lpad(num, 3, "0"))"; kwargs...)
function load_twente_eeg(eeg_name::String; exclude=["ECG", "Stimuli", "Stimulus", "Test", "Unspec", "Buf", "EOG", "EKG"])
    eeg_name = splitext(eeg_name)[1]
    edf = EDF.read(datadir("exp_raw", "twente", "$(eeg_name).edf"))

    annotations_idx = isa.(edf.signals, EDF.AnnotationsSignal)
    annotations = get_relevant_annotations(edf.signals[annotations_idx] |> only)
    seizure_annotations = get_twente_seizure_annotations(annotations)
    artifact_annotations = get_twente_artifact_annotations(annotations)

    @error "Unclear what artifact grades are excluded from Twente dataset"
    ProcessedEEG(edf; 
        exclude_channels=[exclude..., twente_eeg_bad_channels[eeg_name]...], 
        seizure_annotations=seizure_annotations, 
        artifact_annotations=artifact_annotations, 
        mains_hz=50
    )
end

function all_channel_names()
    reduce(map((num) -> load_twente_eeg(num).labels, 1:50); init=[]) do list1, list2
        unique(vcat(list1, list2))
    end
end

function all_annotations()
    reduce(
        map(1:50) do num
            edf = EDF.read(
                datadir("exp_raw", "twente", "EEG$(lpad(num, 3, "0")).edf")
            )
            annotations_idx = isa.(edf.signals, EDF.AnnotationsSignal)
            tals = get_relevant_annotations(edf.signals[annotations_idx] |> only)
            l_annotations = tals .|> (tal) -> tal.annotations
            vcat(l_annotations...)
        end
        ; init=[]) do list1, list2
            unique(lowercase.(vcat(list1, list2)))
    end
end

function load_most_recent_jld2(match_str, dir)
    all_saves = mapreduce(vcat, walkdir(dir)) do (root, dirs, files)
        basenames = filter(Base.Fix1(occursin, "jld2"), files)
        matching_names = filter(Base.Fix1(occursin, match_str), basenames)
        joinpath.(Ref(root), matching_names)
    end
    if isempty(all_saves)
        return nothing
    end
    most_recent_save = sort!(all_saves; by=x -> stat(x).mtime)[end]
    @warn "Loading $most_recent_save..."
    jld_dict = load(most_recent_save)
    @warn "done."
    return jld_dict
end

function save_mat(eeg, matname)
    matdct = Dict(
        "signals" => eeg.signals,
        "channel_names" => eeg.labels,
        "sample_rate" => eeg.sample_rate,
        "seizure_reviewers_count" => eeg.seizure_reviewers_count,
        "seizure_start_sec" => eeg.seizure_annotations .|> first,
        "artifact_start_sec" => eeg.artifact_annotations .|> first,
        "seizure_stop_sec" => eeg.seizure_annotations .|> last,
        "artifact_stop_sec" => eeg.artifact_annotations .|> last
    )
    matwrite(matname, matdct)
end