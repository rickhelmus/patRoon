
library(patRoon)

# -------------------------
# initialization
# -------------------------

workPath <- "test_temp/test-np/features-piek"
setwd(workPath)

# NOTE: please set to a valid data.frame with analysis information. See ?`analysis-information` for more details.
anaInfo <- data.frame(path_centroid = character(), analysis = character(), replicate = character(), blank = character())

# -------------------------
# features
# -------------------------

# Find all features
# NOTE: see the reference manual for many more options
genEICParams <- getPiekEICParams(filter = "none", mzRange = c(80, 800), mzStep = 0.02)
fList <- findFeatures(anaInfo, "piek", genEICParams = genEICParams, peakParams = getDefPeakParams("chrom", "piek"))

# Group and align features between analyses
fGroups <- groupFeatures(fList, "openms", rtalign = TRUE)

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000, relMinReplicateAbundance = 1,
                  maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE, retentionRange = NULL,
                  mzRange = NULL)

# Update group properties
fGroups <- updateGroups(fGroups, what = c("ret", "mz", "mobility"), intWeight = FALSE)

# -------------------------
# reporting
# -------------------------

# Advanced report settings can be edited in the report.yml file.
report(fGroups, MSPeakLists = NULL, formulas = NULL, compounds = NULL, components = NULL, settingsFile = "report.yml",
       openReport = TRUE)
