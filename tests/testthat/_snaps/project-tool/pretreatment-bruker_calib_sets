
library(patRoon)

# -------------------------
# initialization
# -------------------------

workPath <- "test_temp/test-np/pretreatment-bruker_calib_sets"
setwd(workPath)

# NOTE: please set to a valid data.frame with analysis information. See ?`analysis-information` for more details.
anaInfoPos <- data.frame(path_centroid = character(), analysis = character(), replicate = character(),
                         blank = character())
anaInfoNeg <- data.frame(path_centroid = character(), analysis = character(), replicate = character(),
                         blank = character())

# Set to FALSE to skip data pre-treatment
doDataPretreatment <- TRUE
if (doDataPretreatment)
{
    setDAMethod(anaInfoPos, "method-pos")
    recalibrarateDAFiles(anaInfoPos)

    setDAMethod(anaInfoNeg, "method-neg")
    recalibrarateDAFiles(anaInfoNeg)
}

# -------------------------
# features
# -------------------------

# Find all features
# NOTE: see the reference manual for many more options
fListPos <- findFeatures(anaInfoPos, "openms", noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, minFWHM = 1,
                         maxFWHM = 30)
fListNeg <- findFeatures(anaInfoNeg, "openms", noiseThrInt = 1000, chromSNR = 3, chromFWHM = 5, minFWHM = 1,
                         maxFWHM = 30)
fList <- makeSet(fListPos, fListNeg, adducts = c("[M+H]+", "[M-H]-"))

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
