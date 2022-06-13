quickactivate(@__DIR__, "NeonateTriCorr")

using GLMakie, Makie
GLMakie.activate!()
set_window_config!(framerate=3.)

include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics
using TriCorrApplications

f = Figure(resolution=(600,1800))

lbtb = Textbox(f[2,1], placeholder = "Lower bound (Int)", validator=Int, tellwidth=false)
ubtb = Textbox(f[2,2], placeholder = "Upper bound", validator=Float64, tellwidth=false)

lb = Observable(0)
ub = Observable(Inf)

eeg = load_helsinki_eeg(50)

snipped_eeg = Observable(snip(eeg, lb[], ub[]))
times = Observable(get_times(snipped_eeg[]))

on(lbtb.stored_string) do s
    lb[] = parse(Int, s)
    newly_snipped = snip(eeg, lb[], ub[])
    times.val = get_times(newly_snipped)
    snipped_eeg[] = newly_snipped
end
on(ubtb.stored_string) do s
    ub[] = parse(Float64, s)
    newly_snipped = snip(eeg, lb[], ub[])
    times.val = get_times(newly_snipped)
    snipped_eeg[] = newly_snipped
end

draw_eeg_traces!(f[1,:], times, snipped_eeg; downsample_factor=100)
f
# trunc_n_seconds = 4
# truncated_start = round(Int, 10*eeg.sample_rate)
# truncated_len = round(Int, trunc_n_seconds * eeg.sample_rate)
# #truncated_eeg = ProcessedEEG(eeg.signals[:, truncated_start:truncated_start+truncated_len], eeg.labels, eeg.sample_rate, truncated_len / eeg.sample_rate)

# draw_eeg_traces(truncated_eeg)