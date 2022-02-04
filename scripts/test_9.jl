using Distributed
using DrWatson
quickactivate("NeonateTriCorr")
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

using EDF, DSP, Statistics, StatsBase

using Random, JLD2

using HypothesisTests, CSV

eeg = load_helsinki_eeg(9)

using BenchmarkTools

first_second_contributions = @btime calc_class_contributions($eeg, $A_01norm;
        λ_max = (8,25),
        n_motif_classes = 14, 
        n_seconds = 2,
        snippets_start_sec=[0],
        snippets_duration=1
    )

 