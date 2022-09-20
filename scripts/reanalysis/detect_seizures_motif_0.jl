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
all_patient_results_df = load(datadir("motif_results_df_reviewers$(min_reviewers_per_seizure).jld2"))["motif_results_df"]
snippets_duration_s = 1
contributions_spec = "tricorr_ts_zscore_zscore_IndStdNormal_None_snippets$(snippets_duration_s)_lagextents8x25_helsinkiEEG"
save_dir = plotsdir("motif_0_$(contributions_spec)_reviewers$(min_reviewers_per_seizure)_$(Dates.now())")
mkpath(save_dir)
epoch_s = 60
rolling_window_s = 60

drws = mapreduce(vcat, [1:15..., 19,31,44,47,50,62]) do patient_num

    results_df = filter(:patient => p -> p == patient_num, all_patient_results_df)

    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure)

    target_match_str = "$(contributions_spec)$(patient_num)_"

    jld_dict = load_most_recent_jld2(target_match_str, datadir("exp_pro"))
    signals = jld_dict["contributions"]
    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    fig = plot_μ_and_σ_signals_and_roc(results_df, signals, signal_times, seizure_bounds; eeg=eeg, rolling_window_s=rolling_window_s, example_θ=3, n_signals_used=5, epoch_s=epoch_s, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)")

    save(joinpath(save_dir, "roc_patient$(patient_num)_motif0_reviewers$(min_reviewers_per_seizure).png"), fig)
end

end