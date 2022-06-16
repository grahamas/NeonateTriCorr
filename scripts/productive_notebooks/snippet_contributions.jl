
using Distributed
n_cores = length(Sys.cpu_info())
n_threaded_workers = floor(Int, n_cores / Base.Threads.nthreads())
#addprocs(n_threaded_workers - nworkers())
@show nprocs()
@show Base.Threads.nthreads()

#@everywhere begin

#@everywhere using DrWatson
#@everywhere quickactivate(@__DIR__, "NeonateTriCorr")
#@everywhere using Pkg
#@everywhere Pkg.instantiate()

begin

using DrWatson
quickactivate(@__DIR__, "NeonateTriCorr")
using Pkg
#Pkg.instantiate()

include(scriptsdir("include_src.jl"))

using EDF, DSP, Statistics, StatsBase

using TripleCorrelations, ProgressMeter, Random, JLD2

using HypothesisTests, CSV

end # @everywhere begin


using Dates
let n_snippets = 200,
    boundary = Periodic(),
    contribution_desc = "A_znorm_std",#"AN_01norm",
    λ_max=(8,25),
    selected_patients = [PAT...]#,62,75];
sub_dir = if PAT isa Number
    "pat_$(PAT)_snippet_contributions_by_class_$(typeof(boundary))_$(n_snippets)_$(contribution_desc)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
else 
    "snippet_contributions_by_class_$(typeof(boundary))_$(n_snippets)_$(contribution_desc)_$(Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS"))"
end
mkpath(plotsdir(sub_dir))
@everywhere Random.seed!(12345)

eegs = load_helsinki_eeg.(selected_patients) # defined in let

# Plot the traces (optionally histograms)
# for (pat_num, eeg) ∈ zip(selected_patients, eegs)
#     eeg = load_helsinki_eeg(pat_num)

#     # drw = draw_eeg_hists(eeg)
#     # save(plotsdir(sub_dir, "pat$(pat_num)_eeg_hist.png"), drw)

#     drw = draw_eeg_traces(eeg; std_max=7, layout=:layout, downsample_factor=10)
#     save(plotsdir(sub_dir, "pat$(pat_num)_eeg_traces_downsampled_10.png"), drw)
# end


## Plot contributions by case for all classes
@time eeg_contributions_by_case = pmap(zip(selected_patients, eegs)) do (pat_num, eeg)
    # if isempty(eeg.artifact_annotations)
    #     @info "Not replotting $(pat_num): no artifacts noted"
    #     return missing
    # end
    @info "Calculating class contributions for Patient $(pat_num)"
    contributions_by_case = control_vs_seizure_class_contributions(eeg; contributions_desc=contribution_desc, n_snippets=n_snippets, min_dist_to_seizure=60, boundary=boundary, λ_max=λ_max)
    if contributions_by_case |> ismissing
        @warn "Cannot find enough snippets for Patient $(pat_num)"
        return missing
    end
    contributions_by_case
end

@info "Statistical analyses"
eeg_case_statistics_dfs = map(zip(selected_patients, eeg_contributions_by_case)) do (pat_num, contributions_by_case)
    df = control_vs_seizure_all_class_statistics(contributions_by_case.seizure, contributions_by_case.control)
    df.patient_number = ones(Int, nrow(df)) .* pat_num
    df
end

statistics_df = vcat(eeg_case_statistics_dfs...)
CSV.write(plotsdir(sub_dir, "$(contribution_desc)_$(typeof(boundary))_by_class_statistics.csv"), statistics_df)

using CairoMakie
for (pat_num, contributions_by_case) ∈ zip(selected_patients, eeg_contributions_by_case)
    if ismissing(contributions_by_case)
        @info "No analysis for patient $pat_num"
        continue
    end
    @info "Plotting patient $pat_num"
    dct = Dict(map(keys(contributions_by_case)) do case
            contributions = getproperty(contributions_by_case, case)
            signal_rms = putative_signal_rms(contributions)
            (case => signal_rms)
        end...
    )
    fig = boxplot_control_vs_seizure(dct)
    Label(fig[0,:], "Patient $pat_num ($(typeof(boundary)))", tellwidth=false)
    save(plotsdir(sub_dir, "pat$(pat_num)_comparing_$(contribution_desc)_$(typeof(boundary)).png"), fig)
    
    drw = all_class_boxplots_control_vs_seizure(contributions_by_case; show_outliers=true)
   Label(drw.figure[0,:], "Patient $pat_num ($(typeof(boundary)))", tellwidth=false)
   save(plotsdir(sub_dir, "pat$(pat_num)_N$(n_snippets)_by_class_$(contribution_desc)_$(typeof(boundary)).png"), drw)

   drw = all_class_boxplots_control_vs_seizure(contributions_by_case; show_outliers=false)
   Label(drw.figure[0,:], "Patient $pat_num ($(typeof(boundary))", tellwidth=false)
   save(plotsdir(sub_dir, "pat$(pat_num)_N$(n_snippets)_by_class_$(contribution_desc)_no_outliers_$(typeof(boundary)).png"), drw)
end

# res = map([31]) do pat_num #[9,19,21,31,44,47,50,62,75]
#     contributions_by_case = control_vs_seizure_class_contributions(pat_num; n_snippets=n_snippets, min_dist_to_seizure=60, shuffle_norm=true)
#    if contributions_by_case |> ismissing
#        @warn "Cannot find enough snippets for Patient $(pat_num)"
#        return missing
#    end
#    drw = all_class_boxplots_control_vs_seizure(contributions_by_case; show_outliers=true)
#    Label(drw.figure[0,:], "Patient $pat_num", tellwidth=false)
#    save(plotsdir(sub_dir, "pat$(pat_num)_N$(n_snippets)_by_class_$(contribution_desc).png"), drw)

#    drw = all_class_boxplots_control_vs_seizure(contributions_by_case; show_outliers=false)
#    Label(drw.figure[0,:], "Patient $pat_num", tellwidth=false)
#    save(plotsdir(sub_dir, "pat$(pat_num)_N$(n_snippets)_by_class_$(contribution_desc)_no_outliers.png"), drw)

#    contributions_by_case
# end

end
