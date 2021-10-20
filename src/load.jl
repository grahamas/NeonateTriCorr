using MAT, AxisIndices

function load_EEG_snippet(filename)
    test_data_vars = matread(datadir("exp_raw", "patient1_sz1_10secs.mat"))

    expected_vars = Set(["snippetSave"])
    @assert Set(keys(test_data_vars)) == expected_vars

    eeg_snippet_arr = test_data_vars["snippetSave"]
    return NamedAxisArray{(:channel,:time_bin)}(eeg_snippet_arr)
end

