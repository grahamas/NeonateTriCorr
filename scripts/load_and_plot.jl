using MAT, AxisIndices

include(scriptsdir("include_src.jl"))

test_data_vars = matread(datadir("exp_raw", "patient1_sz1_10secs.mat"))

expected_vars = Set(["snippetSave"])
@assert Set(keys(test_data_vars)) == expected_vars

eeg_snippet_arr = test_data_vars["snippetSave"]
eeg_snippet = NamedAxisArray{(:channel,:time_bin)}(eeg_snippet_arr)

plt, drw = plot_eeg_traces(eeg_snippet)
drw