
library(patRoon)

# -------------------------
# initialization
# -------------------------

workPath <- "test_temp/test-np/features-piek_dma_agilent"
setwd(workPath)

# NOTE: please set to a valid data.frame with analysis information. See ?`analysis-information` for more details.
anaInfo <- data.frame(path_centroid = character(), analysis = character(), replicate = character(), blank = character())

CCSParams <- getCCSParams(method = "agilent", calibrant = "")

# -------------------------
# features
# -------------------------

# Find all features
# NOTE: see the reference manual for many more options
genEICParams <- getPiekEICParams(filter = "none", filterIMS = "none", mzRange = c(80, 800), mzStep = 0.08,
                                 mobRange = c(10, 30), mobStep = 1)
fList <- findFeatures(anaInfo, "piek", IMS = TRUE, genEICParams = genEICParams,
                      peakParams = getDefPeakParams("chrom", "piek"))

# Group and align features between analyses
fGroups <- groupFeatures(fList, "greedy", scoreWeights = c(retention = 1, mz = 1, mobility = 1))

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000, relMinReplicateAbundance = 1,
                  maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE, retentionRange = NULL,
                  mzRange = NULL)

# Update group properties
fGroups <- updateGroups(fGroups, what = c("ret", "mz", "mobility"), intWeight = FALSE)

# Assign CCS values
fGroups <- assignMobilities(fGroups, CCSParams = CCSParams)

# -------------------------
# reporting
# -------------------------

# Advanced report settings can be edited in the report.yml file.
report(fGroups, MSPeakLists = NULL, formulas = NULL, compounds = NULL, components = NULL, settingsFile = "report.yml",
       openReport = TRUE)
