# Script automatically generated on Tue Feb 17 08:05:32 2026

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
                  maxReplicateIntRSD = 0.75, blankThreshold = 5, removeBlanks = TRUE, retentionRange = NULL,
                  mzRange = NULL)

# Update group properties
fGroups <- updateGroups(fGroups, what = c("ret", "mz", "mobility"), intWeight = FALSE)

# -------------------------
# annotation
# -------------------------

# Retrieve MS peak lists
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroups, avgFeatParams = avgMSListParams, avgFGroupParams = avgMSListParams)
# Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.
mslists <- filter(mslists, MSLevel = 2, absMinIntensity = NULL, relMinIntensity = 0.05, topMostPeaks = 25,
                  maxMZOverPrec = 4)

# Calculate formula candidates
formulas <- generateFormulas(fGroups, mslists, "genform", adduct = "[M+H]+", elements = "CHNOP", oc = FALSE,
                             calculateFeatures = FALSE)
formulas <- estimateIDConfidence(formulas, IDFile = "idlevelrules.yml")

# Calculate compound structure candidates
compounds <- generateCompounds(fGroups, mslists, "metfrag", adduct = "[M+H]+", database = "pubchemlite",
                               maxCandidatesToStop = 2500)
compounds <- addFormulaScoring(compounds, formulas, updateScore = TRUE)

compounds <- estimateIDConfidence(compounds, MSPeakLists = mslists, formulas = formulas, IDFile = "idlevelrules.yml")

# -------------------------
# reporting
# -------------------------

# Advanced report settings can be edited in the report.yml file.
report(fGroups, MSPeakLists = mslists, formulas = formulas, compounds = compounds, components = NULL,
       settingsFile = "report.yml", openReport = TRUE)
