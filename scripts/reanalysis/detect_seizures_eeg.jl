using DrWatson
@quickactivate "NeonateTriCorr"

# Define patient_num before running this script

using EDF, DSP, Statistics, StatsBase, CairoMakie, AlgebraOfGraphics
ext = "png"
using Random, JLD2
using LinearAlgebra
using KernelDensity

include(scriptsdir("include_src.jl"))

let min_reviewers_per_seizure = 3;

all_patient_results_df = load(datadir("channel_results_df_reviewers$(min_reviewers_per_seizure).jld2"))["channel_results_df"]

snippets_duration_s = 1
alert_grace_s=60
rolling_window_s = 60

save_dir = plotsdir("eeg_1Hz_reviewers$(min_reviewers_per_seizure)_$(Dates.now())")
mkpath(save_dir)

drws = mapreduce(vcat, [1:15..., 19,31,44,47,50,62]) do patient_num
    
    results_df = filter(:patient => p -> p == patient_num, all_patient_results_df)

    @info "Loading EEG $(patient_num)..."
    eeg = load_helsinki_eeg(patient_num; min_reviewers_per_seizure = min_reviewers_per_seizure)
    raw_signals = get_signal(eeg)
    channels = 1:size(raw_signals,1)

    seizure_bounds, consensus = load_helsinki_seizure_annotations(patient_num; min_reviewers_per_seizure=min_reviewers_per_seizure)

    # Get signals in 1 Hz
    window = snippets_duration_s * eeg.sample_rate
    signals = mapreduce(hcat, window:window:size(raw_signals,2)) do i
        mean(raw_signals[:,(i-window+1):i]; dims=2)
    end
    @show size(signals)
    signal_times = get_times(eeg, sample_rate=snippets_duration_s)

    fig = plot_μ_and_σ_signals_and_roc(results_df, signals, signal_times, seizure_bounds; eeg=eeg, rolling_window_s=rolling_window_s, example_θ=3, n_signals_used=1, alert_grace_s=alert_grace_s, snippets_duration_s=snippets_duration_s, title="Patient $(patient_num)")

    save(joinpath(save_dir, "roc_patient$(patient_num)_signals_reviewers$(min_reviewers_per_seizure).png"), fig)
end

end