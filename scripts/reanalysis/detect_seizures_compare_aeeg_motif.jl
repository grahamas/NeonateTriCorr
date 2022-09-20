using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

## Saves contributions timeseries to datadir()/exp_pro/motif_class_contribution_timeseries

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

let min_reviewers_per_seizure = 3;

tricorr_snippets_duration_s = 1
aeeg_snippets_duration_s = 15

epoch_s = 60
rolling_window_s = 60

aeeg_rolling_window = Int(rolling_window_s / aeeg_snippets_duration_s)

tricorr_all_patient_results_df = load(datadir("motif_results_df_reviewers$(min_reviewers_per_seizure).jld2"))["motif_results_df"]
aeeg_all_patient_results_df = load(datadir("aeeg_$(aeeg_snippets_duration_s)_lower_margin_results_df_reviewers$(min_reviewers_per_seizure)_window$(aeeg_rolling_window).jld2"))["channel_results_df"]

tricorr_contributions_spec = "tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(tricorr_snippets_duration_s)_lagextents8x25_helsinkiEEG"
save_dir = plotsdir("comparison_$(tricorr_contributions_spec)_aeegsnippets$(aeeg_snippets_duration_s)_epoch$(epoch_s)_rolling$(rolling_window_s)_reviewers$(min_reviewers_per_seizure)_$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, [15]) do patient_num #1:15..., 19,31,44,47,50,62

    aeeg_results_df = filter(:patient => p -> p == patient_num, aeeg_all_patient_results_df)
    tricorr_results_df = filter(:patient => p -> p == patient_num, tricorr_all_patient_results_df)

    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure)
    aeeg = calculate_aeeg(eeg; window_len_s=aeeg_snippets_duration_s)
    aeeg_signals = aeeg_lower_margin(aeeg)'
    aeeg_signal_times = get_times(eeg, sample_rate=1/aeeg_snippets_duration_s)

    target_match_str = "$(tricorr_contributions_spec)$(patient_num)_"

    jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    tricorr_signals = jld_dict["contributions"]
    tricorr_signal_times = get_times(eeg, sample_rate=1/tricorr_snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    fig = plot_μ_comparison((aeeg_results_df, tricorr_results_df), (aeeg_signal_times, aeeg_signals), (tricorr_signal_times, tricorr_signals), seizure_bounds; eeg=eeg, rolling_window_s=rolling_window_s, example_θ=3, n_signals_used=5, epoch_s=epoch_s, aeeg_snippets_duration_s=aeeg_snippets_duration_s, tricorr_snippets_duration_s=tricorr_snippets_duration_s,  title="Patient $(patient_num)")

    save(joinpath(save_dir, "compare_roc_patient$(patient_num)_reviewers$(min_reviewers_per_seizure).png"), fig)
end

end