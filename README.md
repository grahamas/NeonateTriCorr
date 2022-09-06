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
1. Set `PAT` to be the number of a downloaded patient (`PAT=X`), and then run `contributions_timeseries/contributions_patPAT.jl` (alternatively: run `contributions_patPAT_artifacts.jl` to obtain triple correlation for all timepoints, including those annotated as artifacts)
2. To compare the differences between seizure and non-seizure epochs, run `reanalysis/diffs_motifs.jl` (alternatively: with `_artifacts` suffix). 
3. To attempt to detect seizures, run `reanalysis/detect_seizures_motif_0.jl`.
4. Repeat the previous two steps with `diffs_aeeg.jl` and `detect_seizures_aeeg.jl` respectively to run the same analyses on aEEG-transformed recordings.