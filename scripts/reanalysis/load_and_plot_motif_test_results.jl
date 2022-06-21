using DrWatson
@quickactivate "NeonateTriCorr"

# Assumes you have previously run contributions_patPAT.jl
#   and loads most recently calc'd contributions from datadir("exp_pro")

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using HypothesisTests, CSV, Distances, LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

min_reviewers_per_seizure = 3

motif_results_df = load(datadir("motif_results_df_reviewers$(min_reviewers_per_seizure).jld2"))["motif_results_df"] 

fig_significance = draw_significances_plot!(motif_results_df)
fig_Δμ = draw_Δμ_plot(motif_results_df)
fig_Δσ = draw_Δσ_plot(motif_results_df)


now = Dates.now()
save(plotsdir("$(now)_patient_motif_significances_reviewers$(min_reviewers_per_seizure).png"), fig_significance)
save(plotsdir("$(now)_patient_motif_Δμ_reviewers$(min_reviewers_per_seizure).png"), fig_Δμ)
save(plotsdir("$(now)_patient_motif_Δσ_reviewers$(min_reviewers_per_seizure).png"), fig_Δσ)