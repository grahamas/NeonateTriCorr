quickactivate(@__DIR__, "NeonateTriCorr")

using Dates
using GLMakie
GLMakie.activate!()

include(scriptsdir("include_src.jl"))

output_dir = plotsdir("all_eeg_traces_$(Dates.now())")
mkpath(output_dir)

for pat ∈ keys(helsinki_eeg_bad_channels)
    eeg = load_helsinki_eeg(pat)
    fig = draw_eeg_traces(eeg; downsample_factor=100)
    save(joinpath(output_dir, "pat$(pat).png"), fig)
end