
using ProgressMeter, JLD2

include(srcdir("types.jl"))
include(srcdir("load.jl"))
include(srcdir("motif_class.jl"))
include(srcdir("class_contributions.jl"))
include(srcdir("filename_utils.jl"))
include(srcdir("plots.jl"))
include(srcdir("contribution_comparisons.jl"))
include(srcdir("rolling_estimates.jl"))
include(srcdir("detect_seizures.jl"))
include(srcdir("aeeg.jl"))

include(scriptsdir("meats/detect_patient_seizures.jl"))
include(scriptsdir("meats/detect_all_patients_seizures.jl"))
include(scriptsdir("meats/calculate_patient_tricorr.jl"))
include(scriptsdir("meats/calculate_patient_aEEG.jl"))
include(scriptsdir("meats/epoch_differences.jl"))

include(scriptsdir("params/params.jl"))