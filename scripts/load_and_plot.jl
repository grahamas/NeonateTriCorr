quickactivate("NeonateTriCorr")
include(scriptsdir("include_src.jl"))
using GLMakie
GLMakie.activate!()

function load_and_plot_helsinki_eeg(pat_num; downsample_factor=10)

    eeg = load_helsinki_eeg(pat_num)

    draw_eeg_traces(eeg; downsample_factor=downsample_factor)
end

for pat_num