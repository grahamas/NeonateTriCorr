
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

function A_01znorm(snippet, boundary, λ_max)
    normalize_01!(snippet)
    snippet .-= mean(snippet)
    snippet ./= std(snippet)
    actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
    actual_contributions
 end

 function A_z01norm(snippet, boundary, λ_max)
    snippet .-= mean(snippet)
    snippet ./= std(snippet)
    normalize_01!(snippet)
    actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
    actual_contributions
 end

function A_znorm(snippet, boundary, λ_max)
    snippet .-= mean(snippet)
    snippet ./= std(snippet)
    actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
    actual_contributions
end

function A_znorm_std(snippet, boundary, λ_max)
    snippet .-= mean(snippet)
    snippet ./= std(snippet)
    actual_contributions = sequence_class_tricorr(snippet, boundary, λ_max)
    # assumes expectation is zero
    actual_contributions ./ sqrt.(variance_of_standard_normals(boundary, λ_max)) 
end


snippet_contributions_fns = Dict(
    "AN_01norm" => AN_01norm,
    "AN_01norm_power" => AN_01norm_power,
    "AN_znorm" => AN_znorm,
    "A_znorm" => A_znorm,
    "A_znorm_std" => A_znorm_std,
    "A_01znorm" => A_01znorm,
    "A_z01norm" => A_z01norm
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
        boundary, preproc!::Function, postproc!::Function, 
        assumption::TripleCorrelations.AbstractDistributionAssumption,
        condition::TripleCorrelations.AbstractConditional; 
        lag_extents,
        n_motif_classes, 
        snippets_duration=1
    )
    n_seconds = floor(Int, eeg.duration)
    snippets_start_sec=0:snippets_duration:(n_seconds-1)
    eeg_motif_class_contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Union{Float64,Missing}, n_motif_classes, length(snippets_start_sec)))

    snippet_generator = (get_snippet(eeg, start, snippets_duration) for start in snippets_start_sec)
    precalced_postproc! = precalculate(postproc!, assumption, condition, snippet_generator, boundary, lag_extents)

    @threads for i_sec ∈ 1:length(snippets_start_sec)
        snippet = get_snippet(eeg, snippets_start_sec[i_sec], snippets_duration)
        if any(ismissing.(snippet))
            eeg_motif_class_contributions[:,i_sec] .= missing
        else
            processed_snippet = copy(snippet)
            preproc!(processed_snippet, snippet)
            contributions = sequence_class_tricorr(processed_snippet, boundary, lag_extents)
            precalced_postproc!(eeg_motif_class_contributions[:,i_sec], contributions, processed_snippet)
        end
        #ProgressMeter.next!(p)
    end
    # jldsave(datadir("eeg_class_actual_$(λ_max)_$(PAT).jld2"); class_contributions=eeg_motif_class_contributions)
    #plot_contributions(eeg_motif_class_contributions; annotations=annotations, title=PAT)

    eeg_motif_class_contributions
end
