# Reproducing the Neonatal Seizure Detection Results

This note summarizes the concrete results and analysis settings needed to
reproduce the paper "Detecting heterogeneous seizures in newborn infants using
triple correlation" from this repository.

The source manuscript used for this summary is:

- `/home/graham/cloud/Documents/Projects/Manuscripts/NeonatalSeizureDetectionTex/sn-article.tex`

The repo contains the analysis code, but the raw EEG files, generated `.jld2`
results, and generated plots are intentionally not committed. They are expected
under DrWatson-managed `data/` and `plots/` directories.

## Cohorts and Ground Truth

These cohort values are manuscript-reported values.

- Dataset: Helsinki neonatal EEG dataset from Stevenson et al. 2019.
- Sampling rate: 256 Hz.
- Full cohort: 79 patients.
- Reviewer consensus seizure definition: all 3 expert reviewers marked seizure.
- Full-cohort seizure burden:
  - 39 seizing patients.
  - 343 consensus seizures.
  - Seizing-patient seizure count: 8.8 +/- 11.9 seizures.
  - Total recording duration: 402,825 s.
  - Consensus seizure duration: 39,259 s.
  - Per-seizing-patient seizure duration: 1,006 +/- 1,392 s.
- Artifact-annotated subset:
  - 22 patients after discarding patient 21 because 92% of the recording was
    artifactual.
  - 16 seizing patients.
  - 134 consensus seizures.
  - Seizing-patient seizure count: 8.4 +/- 11.9 seizures.
  - Total recording duration: 121,533 s.
  - Consensus seizure duration: 2,898 s.
  - Per-seizing-patient seizure duration: 929 +/- 842 s.
- Excluded channels:
  - Respiratory and ECG channels.
  - Fz and Cz, due to channel-specific whole-recording artifacts in many
    patients.

Audit note: the artifact-annotated subset's reported `2,898 s` consensus
seizure duration does not reconcile with `16` seizing patients at
`929 +/- 842 s` per seizing patient. Treat that value as a manuscript value to
check during a rerun.

The patient sets used by the repo are in `scripts/params/params.jl`:

- `patients_all = 1:79`
- `patients_artifact_annotated = [1:15..., 19,31,44,47,50,62,75]`
- `patients_unannotated = setdiff(patients_all, patients_artifact_annotated)`

## Preprocessing

The EEG loading and preprocessing code is in `src/load.jl`.

For each included signal:

- De-mean per channel.
- Bandstop filter 45-55 Hz with a 6th-order Butterworth filter.
- Bandpass filter 0.1-70 Hz with a 2nd-order Butterworth filter.
- Divide by per-channel standard deviation.

## Core Analysis Parameters

The paper's figure/table parameters are encoded in `scripts/params/params.jl`.

Common parameters:

- `lag_extents = (8,25)`
- `min_reviewers_per_seizure = 3`
- `excluded_artifact_grades = Int[]` for the main full-cohort analysis.
- `epoch_s = 60`
- `rolling_window_s = 60`
- `window_fn = mean`
- `n_θs = 100`
- `discretization_s = 15`
- `standardization = "within"` by default, though Figure 3 evaluates both
  `"within"` and `"across"`.
- `min_snippets_for_comparison = 20`

TriCorr-specific parameters:

- `snippets_duration_s = 1`
- `preproc! = TripleCorrelations.zscore!`
- `postproc! = TripleCorrelations.identity!`
- `assumption = IndStdNormal()`
- `conditioned_on = None()`
- Lags: spatial `-8:8`, temporal `-25:25`, corresponding to about +/-98 ms at
  256 Hz.
- Signal dimensionality: 14 motif-class contribution traces.

aEEG-specific parameters:

- Bandpass: 2-20 Hz.
- Envelope low-pass frequency: 0.32 Hz.
- `snippets_duration_s = 15`
- Lower margin percentile: 0.09.
- Upper margin percentile: 0.93.
- Detector uses only the lower margin.

## Detector Construction

The multi-patient detector path is implemented in
`scripts/meats/detect_all_patients_seizures.jl` and used by
`scripts/figures/fig3_detectors.jl`.

For both TriCorr and aEEG:

- Compute multi-channel signals:
  - TriCorr: 14 motif classes.
  - aEEG: EEG channels' lower margins.
- Z-score either within patient or across all patients.
- Apply a 60 s rolling mean.
- Reduce to one time series with `maxany`:
  - aEEG uses the maximum over channels.
  - TriCorr uses `maxanyabs`, the maximum absolute motif-class deviation.
- Evaluate per 60 s epoch.
- Count a detection when any suprathreshold signal occurs in the epoch.
- Count an epoch as seizing if reviewer consensus occurs within one minute of
  the epoch boundary.

Main detector entry points:

- `scripts/figures/fig3_detectors.jl`: regenerates Figure 3 ROC curves and the
  ROC statistics table.
- `scripts/reanalysis/allpat_allanalyses.jl`: broader sweep across detector
  variants, reviewer thresholds, and patient sets.
- `scripts/reanalysis/allpat_detecttricorr_seizures.jl`: TriCorr-only detector
  run.
- `scripts/reanalysis/allpat_detectaeeg_seizures.jl`: aEEG-only detector run.

Generated outputs are written under `plots/` and cached ROC/results data are
written under `data/`, both of which are ignored by git.

## Published Detector Results

The paper reports ROC statistics for the full 79-patient cohort with
3-reviewer consensus ground truth. `TPR (5.7 FP/Hr)` uses 5.7 false-positive
minutes per hour because that is the single-reviewer false-positive rate when
single-reviewer annotations are compared to reviewer consensus.

| Target | Standardization | Signal | AUC | TPR at 5.7 FP/Hr | FP/Hr at 80% TPR |
| --- | --- | --- | ---: | ---: | ---: |
| patient | across | aEEG | 0.94 | 0.82 | 5.5 |
| patient | across | TriCorr | 0.92 | 0.79 | 6.3 |
| patient | within | aEEG | 0.96 | 0.92 | 2.9 |
| patient | within | TriCorr | 0.98 | 1.00 | 2.1 |
| seizure | across | aEEG | 0.76 | 0.48 | 27.0 |
| seizure | across | TriCorr | 0.84 | 0.58 | 16.0 |
| seizure | within | aEEG | 0.83 | 0.60 | 15.0 |
| seizure | within | TriCorr | 0.81 | 0.52 | 21.0 |

The main comparative claims to reproduce are:

- Per-seizure TriCorr and aEEG performance are broadly similar.
- Across-patient per-seizure standardization is the largest detector separation:
  TriCorr AUC 0.84 vs aEEG AUC 0.76.
- For per-patient detection with within-patient standardization, TriCorr reaches
  perfect sensitivity at a lower false-positive rate than aEEG:
  - TriCorr: 1.00 TPR at 5.7 FP/Hr and 2.1 FP/Hr at 80% TPR.
  - aEEG: 0.92 TPR at 5.7 FP/Hr and 2.9 FP/Hr at 80% TPR.
- Single expert reviewers, treated as detectors against 3-reviewer consensus,
  have about 6 FP/Hr, reported as 5.7 FP/Hr in the table threshold.

## Published Heterogeneity Result

Among the 39 patients with consensus seizures, TriCorr motif-class contribution
distributions differed from interictal distributions as follows by
Kolmogorov-Smirnov testing:

- 13 patients: significant ictal increase in TriCorr.
- 11 patients: significant ictal decrease in TriCorr.
- 3 patients: heterogeneous effects across motif classes.
- 11 patients: no significant motif-class contribution difference.

This is the result behind the paper's claim that TriCorr exposes seizure-signal
heterogeneity that is collapsed by the final simple detector.

## Expected Figure and Table Products

The manuscript's key reproducibility products are:

- Figure 2/example traces:
  - `neonate_chapter/both_contributions_patient7_2022-10-06T23:17:40.090.png`
  - `neonate_chapter/both_contributions_patient62_2022-10-06T23:26:46.703.png`
- Figure 3 ROC curves:
  - `neonate_chapter/compare_aEEG_tricorr_standardizations_patients79_maxany_ROCs_2022-10-05T17:36:19.236.png`
- Figure 4 detector traces:
  - `neonate_chapter/fig4_detection_traces.png`
- ROC statistics table:
  - Generated by `scripts/figures/fig3_detectors.jl` as
    `sensitivity_3rev_all_patients.tex`.

The generated file names include timestamps, so exact names will differ when
rerun. Match the table values above rather than the timestamp.

## Minimal Reproduction Sequence

1. Instantiate the Julia environment as described in `README.md`.
2. Download the Helsinki EDF files with `download_helsinki_eegs(1:79)`.
3. Generate per-patient TriCorr contribution files with
   `scripts/contributions_timeseries/contributions_patPAT.jl` for each patient.
4. Generate per-patient aEEG files with
   `scripts/contributions_timeseries/aEEG.jl` or the corresponding meat function
   `calculate_patient_aEEG`.
5. Run `scripts/figures/fig3_detectors.jl` to regenerate the ROC curves and
   `sensitivity_3rev_all_patients.tex`.
6. Compare the regenerated detector table to the published values in this note.
