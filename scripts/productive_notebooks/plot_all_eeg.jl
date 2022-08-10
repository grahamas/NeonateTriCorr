quickactivate(@__DIR__, "NeonateTriCorr")

using Dates
using CairoMakie
ext = "png"; CairoMakie.activate!()

include(scriptsdir("include_src.jl"))

output_dir = plotsdir("all_eeg_traces_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))")
mkpath(output_dir)

download_helsinki_eegs(1:75)
for pat ∈ 1:75# keys(helsinki_eeg_bad_channels)
    eeg = load_helsinki_eeg(pat; excluded_artifact_grades=(1,2))
    fig = draw_eeg_traces(eeg; downsample_factor=100, resolution=(800,1800))
    save(joinpath(output_dir, "pat$(pat)_grades12.$(ext)"), fig)
end
