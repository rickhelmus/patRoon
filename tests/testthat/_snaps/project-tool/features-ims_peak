
library(patRoon)

# -------------------------
# initialization
# -------------------------

workPath <- "test_temp/test-np/features-ims_peak"
setwd(workPath)

# NOTE: please set to a valid data.frame with analysis information. See ?`analysis-information` for more details.
anaInfo <- data.frame(path_centroid = character(), analysis = character(), replicate = character(), blank = character())

# -------------------------
# features
# -------------------------

# Find all features
# NOTE: see the reference manual for many more options
fList <- findFeatures(anaInfo, "openms", noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, minFWHM = 1, maxFWHM = 30)

# Group and align features between analyses
fGroups <- groupFeatures(fList, "openms", rtalign = TRUE)

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000, relMinReplicateAbundance = 1,
                  maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE, retentionRange = NULL,
                  mzRange = NULL)

# Update group properties
fGroups <- updateGroups(fGroups, what = c("ret", "mz", "mobility"), intWeight = FALSE)

# -------------------------
# mobility assignment
# -------------------------

fGroups <- assignMobilities(fGroups, mobPeakParams = getDefPeakParams(type = "bruker_ims", algorithm = "xcms"),
                            chromPeakParams = getDefPeakParams(type = "chrom", algorithm = "envipick"),
                            fallbackEIC = TRUE, calcArea = "integrate")

# -------------------------
# reporting
# -------------------------

# Advanced report settings can be edited in the report.yml file.
report(fGroups, MSPeakLists = NULL, formulas = NULL, compounds = NULL, components = NULL, settingsFile = "report.yml",
       openReport = TRUE)
