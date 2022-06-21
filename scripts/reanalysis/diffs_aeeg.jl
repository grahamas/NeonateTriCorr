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

snippets_duration_s = 15
rolling_window_s = 60
rolling_window = ceil(Int, rolling_window_s / snippets_duration_s)
min_reviewers_per_seizure=3
min_snippets_for_comparison = 8

save_dir = plotsdir("aeeg_$(snippets_duration_s)_rolling$(rolling_window_s)_$(Dates.now())")
mkpath(save_dir)
#
channel_results = mapreduce(vcat, [1:15..., 19,31,44,47,50,62]) do PAT

rolling_window = ceil(Int, rolling_window_s / snippets_duration_s)

eeg = load_helsinki_eeg(PAT)

# MUST DEFINE min_reviewers_per_seizure

min_dist_to_seizure = 30 / snippets_duration_s

seizure_bounds, consensus = load_helsinki_seizure_annotations(PAT; min_reviewers_per_seizure=min_reviewers_per_seizure)
raw_signal = get_signal(eeg)
channels = 1:size(raw_signal,1)

@info "calculate aeeg..."
# Get signal in 1 Hz
aeeg = calculate_aeeg(eeg; window_len_s=snippets_duration_s)
@info "...lower margin..."
signal = aeeg_lower_margin(aeeg)'
@info "...times..."
signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)
@info "done."

@show size(signal_times) size(signal)
fig_lower_aeeg = plot_contributions(eeg, signal_times, signal; resolution=(1000,1600), get_label=x -> eeg.labels[x], title="aEEG lower margin ($(snippets_duration_s)s)")
save(joinpath(save_dir, "lower_aeeg_trace_patient$(PAT).png"), fig_lower_aeeg)



@show rolling_window

rolling_std_z = mapreduce(hcat, channels) do channel_num
    stds = [std(signal[channel_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,2))]
    if all(ismissing.(stds))
        missing
    else
        (stds .- mean(skipmissing(stds))) ./ std(skipmissing(stds))
    end
end
rolling_mean_z = mapreduce(hcat, channels) do channel_num
    means = [mean(signal[channel_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,2))]
    if all(ismissing.(means))
        missing
    else
        (means .- mean(skipmissing(means))) ./ std(skipmissing(means))
    end
end

fig_std = plot_contributions(eeg, get_times(eeg,sample_rate=1/snippets_duration_s)[rolling_window:end], rolling_std_z'; title="Rolling σ (zscored per motif; window = $(rolling_window))", resolution=(800,1400), get_label=x -> eeg.labels[x])
save(joinpath(save_dir, "aeeg_lower_margin_standard_deviations_zscored_window$(rolling_window)_pat$(PAT).png"), fig_std)

fig_mean = plot_contributions(eeg, get_times(eeg,sample_rate=1/snippets_duration_s)[rolling_window:end], rolling_mean_z'; title="Rolling μ (zscored per motif; window = $(rolling_window))", resolution=(800,1400), get_label=x -> eeg.labels[x])
save(joinpath(save_dir, "aeeg_lower_margin_means_zscored_window$(rolling_window)_pat$(PAT).png"), fig_mean)

rolling_std = mapreduce(hcat, channels) do motif_num
    stds = [std(signal[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,2))]
end
rolling_mean = mapreduce(hcat, channels) do motif_num
    means = [mean(signal[motif_num,(i-rolling_window+1):(i)]) for i ∈ rolling_window:(size(signal,2))]
end

fig_std = plot_contributions(eeg, get_times(eeg,sample_rate=1/snippets_duration_s)[rolling_window:end], rolling_std'; title="Rolling σ (window = $(rolling_window))", resolution=(800,1400), get_label=x -> eeg.labels[x])
save(joinpath(save_dir, "aeeg_lower_margin_standard_deviations_window$(rolling_window)_pat$(PAT).png"), fig_std)

fig_mean = plot_contributions(eeg, get_times(eeg,sample_rate=1/snippets_duration_s)[rolling_window:end], rolling_mean'; title="Rolling μ (window = $(rolling_window))", resolution=(800,1400), get_label=x -> eeg.labels[x])
save(joinpath(save_dir, "aeeg_lower_margin_means_window$(rolling_window)_pat$(PAT).png"), fig_mean)






if isempty(seizure_bounds)
    @warn "Patient $(PAT): no seizures."
    return []
end

control_snippets = calc_control_snippet_starts_incl_artifacts(seizure_bounds, eeg.duration, signal_times, min_dist_to_seizure)
seizure_snippets = mapreduce((x,y) -> x .|| y, seizure_bounds) do (on, off)
    on .<= (signal_times .+ 1)< off
end


if (count(control_snippets) < min_snippets_for_comparison)
    @warn "Not enough control snippets ($(count(control_snippets))) in Patient $PAT"
    return []
end
if (count(seizure_snippets) < min_snippets_for_comparison) 
    @warn "Not enough seizure snippets ($(count(seizure_snippets))) in Patient $PAT"
    return []
end


control_signal = signal[:, control_snippets]
seizure_signal = signal[:, seizure_snippets]

if sum(.!ismissing.(seizure_signal)) < (14 * min_snippets_for_comparison)
    @warn "Not enough NON-ARTIFACTUAL seizure snippets ($(length(seizure_snippets))) in Patient $PAT"
    return []
end
@info "Calculate effects"
channel_JSDistances_KSpvalues_effect = map(channels) do channel_num
    channel_control_signal = control_signal[channel_num, :] |> skipmissing |> collect
    channel_seizure_signal = seizure_signal[channel_num, :] |> skipmissing |> collect

    (JS=estimate_JSDistance(channel_control_signal, channel_seizure_signal; show_plots=false),
     p=pvalue(ApproximateTwoSampleKSTest(channel_control_signal, channel_seizure_signal)),
     Δμ=mean(channel_seizure_signal) - mean(channel_control_signal),
     Δσ=std(channel_seizure_signal) - std(channel_control_signal),
     patient=PAT,
     motif=eeg.labels[channel_num])
end # map(channels)
@info "done."


@info "plot..."
fig_kdes = plot_estimated_distributions(control_signal, seizure_signal; title="Patient = $(PAT)", get_label = x -> eeg.labels[x])
save(joinpath(save_dir, "estimated_distributions_pat$(PAT).png"), fig_kdes)
@info "done."

channel_JSDistances_KSpvalues_effect

end # mapreduce(PATs)
eeg9 = load_helsinki_eeg(9)
@info "plot summary..."
channel_results_df = DataFrame(channel_results)
fig_significance = draw_significances_plot!(channel_results_df; all_motifs=eeg9.labels)
fig_Δμ = draw_Δμ_plot(channel_results_df; all_motifs=eeg9.labels)
fig_Δσ = draw_Δσ_plot(channel_results_df; all_motifs=eeg9.labels)
@info "done. save..."

now = Dates.now()
save(plotsdir("$(now)_aeeg_lower_margin_significances_reviewers$(min_reviewers_per_seizure)_window$(rolling_window).png"), fig_significance)
save(plotsdir("$(now)_aeeg_lower_margin_Δμ_reviewers$(min_reviewers_per_seizure)_window$(rolling_window).png"), fig_Δμ)
save(plotsdir("$(now)_aeeg_lower_margin_Δσ_reviewers$(min_reviewers_per_seizure)_window$(rolling_window).png"), fig_Δσ)

@save datadir("aeeg_$(snippets_duration_s)_lower_margin_results_df_reviewers$(min_reviewers_per_seizure)_window$(rolling_window).jld2") channel_results_df

@info "done."
