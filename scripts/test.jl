using Distributed
using DrWatson
@quickactivate "NeonateTriCorr"
# n_cores = length(Sys.cpu_info())
# n_threaded_workers = floor(Int, n_cores / Base.Threads.nthreads())
# #addprocs(n_threaded_workers - nworkers())
# @show nprocs()
# @everywhere using DrWatson
# @everywhere quickactivate("NeonateTriCorr")
# @everywhere using Pkg
# @everywhere Pkg.instantiate()

# @everywhere begin


include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics, StatsBase, CairoMakie
ext = "png"

using Random, JLD2

using HypothesisTests, CSV

eeg = load_helsinki_eeg(50)
eeg_snip = snip(eeg, 3)

using BenchmarkTools

test_contributions = calc_class_contributions(eeg_snip, Periodic(), AN_01norm;
        λ_max = (8,25),
        n_motif_classes = 14,
        snippets_duration=1
    )

save("test_AN_contributions_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).jld2", Dict("contributions" => test_contributions))

fig = plot_contributions(test_contributions; eeg=eeg_snip)
save(plotsdir("test_AN_contributions_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS")).$(ext)", fig)
 