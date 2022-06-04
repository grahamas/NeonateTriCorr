
function TriCorrApplications.plot_reviewer_consensus!(ax, eeg::AbstractProcessedEEG)
    lines!(ax, get_times(eeg, sample_rate=1), eeg.seizure_reviewers_count)
    tightlimits!(ax); hidespines!(ax)
    hidedecorations!(ax, ticklabels=false)
end