
include(scriptsdir("include_src.jl"))

eeg_snippet = load_EEG_snippet(datadir("exp_raw", "patient1_sz1_10secs.mat"))

plt, drw = plot_eeg_traces(eeg_snippet)
drw