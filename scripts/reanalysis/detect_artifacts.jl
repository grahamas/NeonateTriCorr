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
run_iterations = isdefined(Main, :ARTIFACT_RUN_ITERATIONS) ? ARTIFACT_RUN_ITERATIONS : true
bin_s = isdefined(Main, :ARTIFACT_BIN_S) ? ARTIFACT_BIN_S : 15
artifact_grades = isdefined(Main, :ARTIFACT_GRADES) ? ARTIFACT_GRADES : [1, 2]
min_reviewers_per_seizure = isdefined(Main, :ARTIFACT_MIN_REVIEWERS_PER_SEIZURE) ? ARTIFACT_MIN_REVIEWERS_PER_SEIZURE : 3
n_thresholds = isdefined(Main, :ARTIFACT_N_THRESHOLDS) ? ARTIFACT_N_THRESHOLDS : 100
fixed_fraction = isdefined(Main, :ARTIFACT_FIXED_FRACTION) ? ARTIFACT_FIXED_FRACTION : 0.05
save_outputs = isdefined(Main, :ARTIFACT_SAVE_OUTPUTS) ? ARTIFACT_SAVE_OUTPUTS : true

if run_iterations
    hvg_sample_strides = isdefined(Main, :ARTIFACT_HVG_SAMPLE_STRIDES) ? ARTIFACT_HVG_SAMPLE_STRIDES : (1, 2, 4)
    nvg_sample_strides = isdefined(Main, :ARTIFACT_NVG_SAMPLE_STRIDES) ? ARTIFACT_NVG_SAMPLE_STRIDES : (4, 8, 16)
    context_specs = isdefined(Main, :ARTIFACT_CONTEXT_SPECS) ? ARTIFACT_CONTEXT_SPECS : ((:max, 1), (:max, 2), (:mean, 1))
    methods = isdefined(Main, :ARTIFACT_METHODS) ? ARTIFACT_METHODS : artifact_iteration_methods(
        hvg_sample_strides=hvg_sample_strides,
        nvg_sample_strides=nvg_sample_strides,
        context_specs=context_specs
    )
    classifier_lambda = isdefined(Main, :ARTIFACT_CLASSIFIER_LAMBDA) ? ARTIFACT_CLASSIFIER_LAMBDA : 1.0
    run_classifier = isdefined(Main, :ARTIFACT_RUN_CLASSIFIER) ? ARTIFACT_RUN_CLASSIFIER : true
    output_root = isdefined(Main, :ARTIFACT_OUTPUT_ROOT) ? ARTIFACT_OUTPUT_ROOT : datadir("exp_pro", "artifact_detection_iterations")
    plot_root = isdefined(Main, :ARTIFACT_PLOT_ROOT) ? ARTIFACT_PLOT_ROOT : plotsdir("artifact_detection_iterations_$(Dates.now())")

    artifact_detection_results = run_artifact_iteration_survey(;
        patients=patients,
        methods=methods,
        bin_s=bin_s,
        artifact_grades=artifact_grades,
        min_reviewers_per_seizure=min_reviewers_per_seizure,
        n_thresholds=n_thresholds,
        fixed_fraction=fixed_fraction,
        classifier_lambda=classifier_lambda,
        run_classifier=run_classifier,
        output_root=output_root,
        plot_root=plot_root,
        save_outputs=save_outputs
    )
else
    visibility_graph_sample_stride = isdefined(Main, :ARTIFACT_VISIBILITY_GRAPH_SAMPLE_STRIDE) ? ARTIFACT_VISIBILITY_GRAPH_SAMPLE_STRIDE : 4
    methods = isdefined(Main, :ARTIFACT_METHODS) ? ARTIFACT_METHODS : artifact_detection_methods(
        visibility_graph_sample_stride=visibility_graph_sample_stride
    )
    output_root = isdefined(Main, :ARTIFACT_OUTPUT_ROOT) ? ARTIFACT_OUTPUT_ROOT : datadir("exp_pro", "artifact_detection")
    plot_root = isdefined(Main, :ARTIFACT_PLOT_ROOT) ? ARTIFACT_PLOT_ROOT : plotsdir("artifact_detection_$(Dates.now())")

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
end

artifact_detection_results
