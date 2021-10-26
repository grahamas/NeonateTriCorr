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
