# NeonateTriCorr

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> NeonateTriCorr

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Data pipeline

To obtain triple correlations of recordings from the Helsinki dataset use the following steps. Generally, for scripts ending in `.jl`, the instruction "run `X.jl`" means "type `include("path/to/X.jl")`." Best practice is to run `Y/X.jl` as `include(scriptsdir("Y","X.jl"))`.

0. Activate this project (such as by typing `Pkg.activate("path/to/this/project")`)
1. Download patient recordings with `download_helsinki_eegs(patient_numbers::Vector{Int})`
2. Set `PAT` to be the number of a downloaded patient (`PAT=X`), and then run `contributions_timeseries/contributions_patPAT.jl` (alternatively: run `contributions_patPAT_artifacts.jl` to obtain triple correlation for all timepoints, including those annotated as artifacts). I typically ran this using SLURM on a cluster, so that the contributions were computed in parallel jobs.
3. To compare the differences between seizure and non-seizure epochs, run `reanalysis/diffs_tricorr.jl` (alternatively: `reanalysis/diffs_tricorr_artifacts.jl`).
4. To attempt to detect seizures from triple-correlation outputs, run `reanalysis/detecttricorr_seizures.jl` (alternatively: `reanalysis/detecttricorr_seizures_artifacts.jl`).
5. Repeat the previous two steps with `reanalysis/diffs_aeeg.jl` and `reanalysis/detectaeeg_seizures.jl` respectively to run the same analyses on aEEG-transformed recordings (alternatively: `reanalysis/diffs_aEEG_artifacts.jl` and `reanalysis/detectaeeg_seizures_artifacts.jl`).

## Verify checked-in artifacts

The complete data pipeline depends on raw EEG/annotation files and unregistered Julia packages that are not stored in this repository. To verify the repository-local artifacts that are checked in here, install the verifier dependency (`python -m pip install h5py`) if needed and run:

```
python scripts/verification/verify_checked_results.py
```

This script checks that the pipeline scripts referenced above are present, validates the bundled one-second patient snippet (`pat9_snippet_1s.csv`), and verifies that the two checked-in JLD2 contribution result files still expose the expected `10 × 14` `Float64` contribution matrices.

## Reproduce figures

To reproduce the figures, first follow steps 0-2 of the Data Pipeline instructions above. Then the `scripts/figures` scripts will reproduce all figures except the first (which resulted from plotting snippets of preprocessed EEG loaded using the `load_helsinki_eeg` function, exported, and then imported to MATLAB).
