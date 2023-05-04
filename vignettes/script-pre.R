# Script automatically generated on Mon Jul 26 14:32:58 2021

library(patRoon)

# -------------------------
# initialization
# -------------------------

workPath <- "C:/myproject"
setwd(workPath)

# Example data from patRoonData package (triplicate solvent blank + triplicate standard)
anaInfo <- patRoonData::exampleAnalysisInfo("positive")

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
                  maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE,
                  retentionRange = NULL, mzRange = NULL)

# -------------------------
# annotation
# -------------------------

# Retrieve MS peak lists
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 5, precursorMzWindow = 4,
                               avgFeatParams = avgMSListParams,
                               avgFGroupParams = avgMSListParams)
# Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.
mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL, relMSMSIntThr = 0.05,
                  topMSPeaks = NULL, topMSMSPeaks = 25)

# Calculate formula candidates
formulas <- generateFormulas(fGroups, mslists, "genform", relMzDev = 5, adduct = "[M+H]+", elements = "CHNOP",
                             oc = FALSE, calculateFeatures = TRUE,
                             featThresholdAnn = 0.75)

# Calculate compound structure candidates
compounds <- generateCompounds(fGroups, mslists, "metfrag", dbRelMzDev = 5, fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "[M+H]+", database = "pubchem",
                               maxCandidatesToStop = 2500)
compounds <- addFormulaScoring(compounds, formulas, updateScore = TRUE)

# -------------------------
# reporting
# -------------------------

report(fGroups, MSPeakLists = mslists, formulas = formulas, compounds = compounds,
       components = NULL, settingsFile = "report.yml", openReport = TRUE)
