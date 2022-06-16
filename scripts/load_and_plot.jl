quickactivate("NeonateTriCorr")
include(scriptsdir("include_src.jl"))
using CairoMakie
CairoMakie.activate!()

function load_and_plot_helsinki_eeg(pat_num; downsample_factor=100)
    eeg = load_helsinki_eeg(pat_num)
    draw_eeg_traces(eeg; downsample_factor=downsample_factor, resolution=(600,1400))
end

function load_and_plot_helsinki_eeg(pat_num, snip_start, snip_stop; downsample_factor=100)
    eeg = load_helsinki_eeg(pat_num)
    snipped_eeg = snip(eeg, snip_start, snip_stop)
    draw_eeg_traces(snipped_eeg; downsample_factor=downsample_factor, resolution=(600,1400))
end