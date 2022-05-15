
using TripleCorrelations, Base.Threads

##### Contribution Metrics ######

function A_01norm(snippet, boundary, λ_max)
    normalize_01!(snippet)
    actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
    return actual_contributions
end

function AN_01norm(snippet, boundary, λ_max)
    normalize_01!(snippet)
    bootstrap_normed_sequence_classes(snippet, boundary, λ_max; n_bootstraps=2, bootstraps_step=2) .- 1
end

function AN_01norm_power(snippet, boundary, λ_max)
    snippet .^= 2
    normalize_01!(snippet)
    actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
    noise_contributions = sequence_class_tricorr(shuffle(snippet), boundary, λ_max) # FIXME should cache noise contributions... somehow.
    return ((actual_contributions ./ noise_contributions) .- 1)
end

function AN_znorm(snippet, boundary, λ_max)
   n_bootstrap = 5
   snippet .-= mean(snippet)
   snippet ./= std(snippet)
   actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
   noise_contributions = sequence_class_tricorr(shuffle(snippet), boundary, λ_max)
   for _ in 1:(n_bootstrap-1)
	noise_contributions += sequence_class_tricorr(shuffle(snippet), boundary, λ_max)
   end
   noise_contributions ./= n_bootstrap
   actual_contributions .- noise_contributions
end

function A_znorm(snippet, boundary, λ_max)
    snippet .-= mean(snippet)
    snippet ./= std(snippet)
    actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
    actual_contributions
 end

snippet_contributions_fns = Dict(
    "AN_01norm" => AN_01norm,
    "AN_01norm_power" => AN_01norm_power,
    "AN_znorm" => AN_znorm,
    "A_znorm" => A_znorm
)

##### Calculating #####

# function calc_class_contributions(eeg::AbstractProcessedEEG, boundary, contributions_fn::Function; λ_max,
#         n_motif_classes, 
#         n_seconds = floor(Int, eeg.duration),
#         snippets_start_sec=0:(n_seconds-1),
#         snippets_duration=1
#     )
#     eeg_motif_class_contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Float64, n_motif_classes, length(snippets_start_sec)))
#     # sig_lock = ReentrantLock()
#     # class_lock = ReentrantLock()
#     #p = ProgressMeter.Progress(length(snippets_start_sec))
#     @threads for (i_sec, snippet_start_sec) ∈ enumerate(snippets_start_sec)
#         i_start = round(Int, (snippet_start_sec*eeg.sample_rate)+1)
#         i_end = round(Int, (snippet_start_sec+snippets_duration)*eeg.sample_rate)
#         snippet = eeg.signals[:,i_start:i_end]
#         contributions = contributions_fn(snippet, boundary, λ_max)
#         eeg_motif_class_contributions[:,i_sec] .= contributions
        
#         #ProgressMeter.next!(p)
#     end
#     # jldsave(datadir("eeg_class_actual_$(λ_max)_$(PAT).jld2"); class_contributions=eeg_motif_class_contributions)
#     #plot_contributions(eeg_motif_class_contributions; annotations=annotations, title=PAT)

#     eeg_motif_class_contributions
# end

function calc_class_contributions(eeg::AbstractProcessedEEG, 
        boundary, contributions_fn::Function; 
        λ_max,
        n_motif_classes, 
        snippets_duration=1
    )
    n_seconds = floor(Int, eeg.duration)
    snippets_start_sec=0:snippets_duration:(n_seconds-1)
    eeg_motif_class_contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Union{Float64,Missing}, n_motif_classes, length(snippets_start_sec)))
    #p = ProgressMeter.Progress(length(snippets_start_sec))
    @threads for i_sec ∈ 1:length(snippets_start_sec)
        snippet_start_sec = snippets_start_sec[i_sec]
        if !in_artifact(snippet_start_sec, eeg)
            i_start = floor(Int, (snippet_start_sec*eeg.sample_rate)+1)
            i_end = floor(Int, (snippet_start_sec+snippets_duration)*eeg.sample_rate)
            snippet = eeg.signals[:,i_start:i_end]
            contributions = contributions_fn(snippet, boundary, λ_max)
            eeg_motif_class_contributions[:,i_sec] .= contributions
        else
            eeg_motif_class_contributions[:,i_sec] .= missing
        end
        
        #ProgressMeter.next!(p)
    end
    # jldsave(datadir("eeg_class_actual_$(λ_max)_$(PAT).jld2"); class_contributions=eeg_motif_class_contributions)
    #plot_contributions(eeg_motif_class_contributions; annotations=annotations, title=PAT)

    eeg_motif_class_contributions
end
