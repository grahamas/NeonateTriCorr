
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

function calc_seizure_snippet_starts(eeg::AbstractProcessedEEG, snippets_duration)
    return vcat([collect(on:snippets_duration:off) for (on, off) ∈ eeg.seizure_annotations]...)
end

function calc_control_snippet_starts(eeg::AbstractProcessedEEG, snippets_duration, min_dist_to_seizure)
    # should assert that seizure and artifacts do not overlap
    seizure_onsets = first.(eeg.seizure_annotations); artifact_onsets = first.(eeg.artifact_annotations)
    seizure_offsets = last.(eeg.seizure_annotations); artifact_offsets = last.(eeg.artifact_annotations)
    control_onsets = [0, seizure_offsets..., artifact_offsets...]
    sort!(control_onsets)
    control_offsets = [seizure_onsets..., artifact_onsets..., eeg.duration]
    sort!(control_offsets)
    control_bounds = zip(control_onsets, control_offsets)
    return vcat([(on+min_dist_to_seizure):snippets_duration:(off-min_dist_to_seizure) for (on, off) ∈ control_bounds if (off-min_dist_to_seizure)-(on+min_dist_to_seizure)>snippets_duration]...)
end

function calc_class_contributions(eeg::AbstractProcessedEEG, boundary, contributions_desc::String; λ_max,
        n_motif_classes = 14, 
        n_seconds = floor(Int, eeg.duration),
        snippets_start_sec=0:(n_seconds-1),
        snippets_duration=1
    )
    eeg_motif_class_contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Float64, n_motif_classes, length(snippets_start_sec)))
    # sig_lock = ReentrantLock()
    # class_lock = ReentrantLock()
    p = ProgressMeter.Progress(length(snippets_start_sec))
    for (i_sec, snippet_start_sec) ∈ enumerate(snippets_start_sec)
        i_start = round(Int, (snippet_start_sec*eeg.sample_rate)+1)
        i_end = round(Int, (snippet_start_sec+snippets_duration)*eeg.sample_rate)
        snippet = eeg.signals[:,i_start:i_end]
        contributions = snippet_contributions_fns[contributions_desc](snippet, boundary, λ_max)
        eeg_motif_class_contributions[:,i_sec] .= contributions
        
        ProgressMeter.next!(p)
    end
    # jldsave(datadir("eeg_class_actual_$(λ_max)_$(PAT).jld2"); class_contributions=eeg_motif_class_contributions)
    #plot_contributions(eeg_motif_class_contributions; annotations=annotations, title=PAT)

    eeg_motif_class_contributions

end

function control_vs_seizure_class_contributions(eeg::AbstractProcessedEEG; 
    boundary,
    contributions_desc,
    n_snippets, snippets_duration=1, min_dist_to_seizure=60, kwargs...)
    seizure_snippet_starts = calc_seizure_snippet_starts(eeg, snippets_duration)
    control_snippet_starts = calc_control_snippet_starts(eeg, snippets_duration, min_dist_to_seizure)

    if (length(control_snippet_starts) <= n_snippets)
        @warn "Patient lacks sufficiently many separated snippets."
        return missing
    end

    seizure_snippet_starts = sort(sample(seizure_snippet_starts, min(n_snippets, length(seizure_snippet_starts)), replace=false))
    control_snippet_starts = sort(sample(control_snippet_starts, n_snippets, replace=false))

    seizure_snippets_class_contribution = calc_class_contributions(eeg, boundary, contributions_desc; snippets_start_sec = seizure_snippet_starts, snippets_duration = snippets_duration, kwargs...)
    control_snippets_class_contribution = calc_class_contributions(eeg, boundary, contributions_desc; snippets_start_sec = control_snippet_starts, snippets_duration = snippets_duration, kwargs...)

    return (seizure=seizure_snippets_class_contribution, control=control_snippets_class_contribution)
end

function putative_signal_rms(arr; putative_signal_classes=2:14)
    arr = NamedDimsArray{(:motif_class, :time)}(arr)
    signal_arr = arr[putative_signal_classes,:]
    dropdims(sqrt.(mean(signal_arr .^ 2, dims=1)), dims=:motif_class)
end

function all_class_boxplots_control_vs_seizure(dct; show_outliers=true)
    control = dct[:control]; seizure = dct[:seizure]
    header = [offset_motif_numeral.(1:14)..., "case"]
    values = [control'; seizure']
    labels = [[1 for _ ∈ 1:size(control,:time)]...; [2 for _ ∈ 1:size(seizure,:time)]]
    labeled_values = hcat(values, labels)
    df = DataFrame(labeled_values, header)
    stacked_df = stack(df, 1:14)
    plt = data(stacked_df) * mapping(:case, :value => "A/N - 1", layout=:variable) * visual(BoxPlot; show_outliers=show_outliers)
    axis = (xticks = (1:2, ["control", "seizure"]), height=120, width=120)
    facet = (; linkyaxes = :none)
    draw(plt; axis, facet)
end

function boxplot_control_vs_seizure(dct)
    control = NamedDimsArray{(:time,)}(dct[:control]); 
    seizure = NamedDimsArray{(:time,)}(dct[:seizure])
    values = [control...; seizure...]
    labels = [[1 for _ ∈ control]...; [2 for _ ∈ seizure]]
    fig = Figure(); ax = Axis(fig[1,1])
    boxplot!(ax, labels, values; show_outliers=false)
    ax.xticks = (1:2, ["control", "seizure"])
    fig
end

function control_vs_seizure_all_class_statistics(seizure::NamedDimsArray{(:motif_class,:time)}, control::NamedDimsArray{(:motif_class,:time)})
    DataFrame(map(axes(seizure, :motif_class)) do class_i
        seizure_i = seizure[class_i,:]
        control_i = control[class_i,:]
        return (
            class = offset_motif_numeral(class_i),
            pvalue = pvalue(MannWhitneyUTest(seizure_i, control_i)),
            mean_effect = mean(seizure_i) - mean(control_i)
        )
    end)
end 
end # @everywhere begin


using Dates
let n_snippets = 200,
    boundary = Periodic(),
    contribution_desc = "A_znorm",#"AN_01norm",
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
