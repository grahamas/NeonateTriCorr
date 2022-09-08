

using DrWatson
@quickactivate "NeonateTriCorr"

# Define PAT before running this script

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))


let snippets_duration_s = 15;
min_reviewers_per_seizure = 3
rolling_window_s = 60
alert_grace_s = 60
rolling_window = Int(rolling_window_s / snippets_duration_s)
all_patient_results_df = load(datadir("aeeg_$(snippets_duration_s)_lower_margin_results_df_reviewers$(min_reviewers_per_seizure)_window$(rolling_window).jld2"))["channel_results_df"]

save_dir = plotsdir("aeeg_artifacts_$(snippets_duration_s)_lower_margin_reviewers$(min_reviewers_per_seizure)_$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, [1:15..., 19,31,44,47,50,62,75]) do patient_num

    results_df = filter(:patient => p -> p == patient_num, all_patient_results_df)

    @info "Loading aEEG $(patient_num)..."
    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=Int[])
    aeeg = calculate_aeeg(eeg; window_len_s=snippets_duration_s)
    signals = aeeg_lower_margin(aeeg)'
    signal_times = get_times(eeg, sample_rate=1/snippets_duration_s)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    # send eeg with artifacts indicated for plotting
    eeg_artifacts_highlighted = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure, excluded_artifact_grades=Int[1])
    
    fig = plot_μ_and_σ_signals_and_roc(signals, signal_times, seizure_bounds; analysis_eeg=eeg, plot_eeg=eeg_artifacts_highlighted, rolling_window_s=rolling_window_s, example_θ=3, n_signals_used=5, alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)")

    save(joinpath(save_dir, "roc_patient$(patient_num)_aeeg_reviewers$(min_reviewers_per_seizure)_artifacts.png"), fig)
end
end