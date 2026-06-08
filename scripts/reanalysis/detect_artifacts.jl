using DrWatson
@quickactivate "NeonateTriCorr"

using CairoMakie
using CSV
using DataFrames
using Dates
using DSP
using EDF
using Statistics

include(scriptsdir("include_src.jl"))

patients = isdefined(Main, :ARTIFACT_PATIENTS) ? ARTIFACT_PATIENTS : artifact_labeled_patients()
visibility_graph_sample_stride = isdefined(Main, :ARTIFACT_VISIBILITY_GRAPH_SAMPLE_STRIDE) ? ARTIFACT_VISIBILITY_GRAPH_SAMPLE_STRIDE : 4
methods = isdefined(Main, :ARTIFACT_METHODS) ? ARTIFACT_METHODS : artifact_detection_methods(
    visibility_graph_sample_stride=visibility_graph_sample_stride
)
bin_s = isdefined(Main, :ARTIFACT_BIN_S) ? ARTIFACT_BIN_S : 15
artifact_grades = isdefined(Main, :ARTIFACT_GRADES) ? ARTIFACT_GRADES : [1, 2]
min_reviewers_per_seizure = isdefined(Main, :ARTIFACT_MIN_REVIEWERS_PER_SEIZURE) ? ARTIFACT_MIN_REVIEWERS_PER_SEIZURE : 3
n_thresholds = isdefined(Main, :ARTIFACT_N_THRESHOLDS) ? ARTIFACT_N_THRESHOLDS : 100
fixed_fraction = isdefined(Main, :ARTIFACT_FIXED_FRACTION) ? ARTIFACT_FIXED_FRACTION : 0.05
output_root = isdefined(Main, :ARTIFACT_OUTPUT_ROOT) ? ARTIFACT_OUTPUT_ROOT : datadir("exp_pro", "artifact_detection")
plot_root = isdefined(Main, :ARTIFACT_PLOT_ROOT) ? ARTIFACT_PLOT_ROOT : plotsdir("artifact_detection_$(Dates.now())")
save_outputs = isdefined(Main, :ARTIFACT_SAVE_OUTPUTS) ? ARTIFACT_SAVE_OUTPUTS : true

artifact_detection_results = run_artifact_detection_survey(;
    patients=patients,
    methods=methods,
    bin_s=bin_s,
    artifact_grades=artifact_grades,
    min_reviewers_per_seizure=min_reviewers_per_seizure,
    n_thresholds=n_thresholds,
    fixed_fraction=fixed_fraction,
    output_root=output_root,
    plot_root=plot_root,
    save_outputs=save_outputs
)

artifact_detection_results
