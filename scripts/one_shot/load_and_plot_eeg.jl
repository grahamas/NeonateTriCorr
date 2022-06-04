quickactivate(@__DIR__, "NeonateTriCorr")

using WGLMakie
WGLMakie.activate!()

include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics

f_low = 0.1
f_high = 70

eeg = load_helsinki_eeg(50)

draw_eeg_traces(eeg; downsample_factor=100)

# trunc_n_seconds = 4
# truncated_start = round(Int, 10*eeg.sample_rate)
# truncated_len = round(Int, trunc_n_seconds * eeg.sample_rate)
# #truncated_eeg = ProcessedEEG(eeg.signals[:, truncated_start:truncated_start+truncated_len], eeg.labels, eeg.sample_rate, truncated_len / eeg.sample_rate)

# draw_eeg_traces(truncated_eeg)