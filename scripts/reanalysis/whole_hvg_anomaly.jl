using DrWatson
@quickactivate "NeonateTriCorr"

using CSV
using DataFrames
using Dates
using JLD2
using Statistics

include(scriptsdir("include_src.jl"))

patients = isdefined(Main, :WHOLE_HVG_PATIENTS) ?
    WHOLE_HVG_PATIENTS :
    artifact_labeled_patients()
bin_s = isdefined(Main, :WHOLE_HVG_BIN_S) ? WHOLE_HVG_BIN_S : 1
artifact_grades = isdefined(Main, :WHOLE_HVG_ARTIFACT_GRADES) ?
    WHOLE_HVG_ARTIFACT_GRADES :
    [1, 2]
min_reviewers_per_seizure = isdefined(Main, :WHOLE_HVG_MIN_REVIEWERS_PER_SEIZURE) ?
    WHOLE_HVG_MIN_REVIEWERS_PER_SEIZURE :
    3
save_node_scores = isdefined(Main, :WHOLE_HVG_SAVE_NODE_SCORES) ?
    WHOLE_HVG_SAVE_NODE_SCORES :
    false
include_edge_span_features = isdefined(Main, :WHOLE_HVG_INCLUDE_EDGE_SPAN_FEATURES) ?
    WHOLE_HVG_INCLUDE_EDGE_SPAN_FEATURES :
    true
include_centrality_features = isdefined(Main, :WHOLE_HVG_INCLUDE_CENTRALITY_FEATURES) ?
    WHOLE_HVG_INCLUDE_CENTRALITY_FEATURES :
    true
include_within_bin_graph_features = isdefined(Main, :WHOLE_HVG_INCLUDE_WITHIN_BIN_GRAPH_FEATURES) ?
    WHOLE_HVG_INCLUDE_WITHIN_BIN_GRAPH_FEATURES :
    true
include_within_bin_topology_features = isdefined(Main, :WHOLE_HVG_INCLUDE_WITHIN_BIN_TOPOLOGY_FEATURES) ?
    WHOLE_HVG_INCLUDE_WITHIN_BIN_TOPOLOGY_FEATURES :
    true
centrality_landmark_count = isdefined(Main, :WHOLE_HVG_CENTRALITY_LANDMARK_COUNT) ?
    WHOLE_HVG_CENTRALITY_LANDMARK_COUNT :
    64
centrality_max_nodes = isdefined(Main, :WHOLE_HVG_CENTRALITY_MAX_NODES) ?
    WHOLE_HVG_CENTRALITY_MAX_NODES :
    WHOLE_HVG_DEFAULT_CENTRALITY_MAX_NODES
output_root = isdefined(Main, :WHOLE_HVG_OUTPUT_ROOT) ?
    WHOLE_HVG_OUTPUT_ROOT :
    datadir("exp_pro", "whole_hvg_anomaly")
save_outputs = isdefined(Main, :WHOLE_HVG_SAVE_OUTPUTS) ?
    WHOLE_HVG_SAVE_OUTPUTS :
    true

whole_hvg_results = run_whole_hvg_anomaly_survey(;
    patients=patients,
    bin_s=bin_s,
    artifact_grades=artifact_grades,
    min_reviewers_per_seizure=min_reviewers_per_seizure,
    save_node_scores=save_node_scores,
    include_edge_span_features=include_edge_span_features,
    include_centrality_features=include_centrality_features,
    include_within_bin_graph_features=include_within_bin_graph_features,
    include_within_bin_topology_features=include_within_bin_topology_features,
    centrality_landmark_count=centrality_landmark_count,
    centrality_max_nodes=centrality_max_nodes,
    output_root=output_root,
    save_outputs=save_outputs
)

whole_hvg_results
