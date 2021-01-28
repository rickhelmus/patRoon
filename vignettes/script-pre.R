## Script automatically generated on Mon Mar 25 08:24:32 2019

library(patRoon)


# -------------------------
# initialization
# -------------------------


workPath <- "C:/myproject"
setwd(workPath)

# Take example data from patRoonData package (triplicate solvent blank + triplicate standard)
anaInfo <- generateAnalysisInfo(paths = patRoonData::exampleDataPath(),
                                groups = c(rep("solvent", 3), rep("standard", 3)),
                                blanks = "solvent")


# -------------------------
# features
# -------------------------


# Find all features.
# NOTE: see manual for many more options
fList <- findFeatures(anaInfo, "openms")

# Group and align features between analysis
fGroups <- groupFeatures(fList, "openms")

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000,
                  relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE,
                  retentionRange = NULL, mzRange = NULL)


# -------------------------
# annotation
# -------------------------


# Retrieve MS peak lists
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 5, precursorMzWindow = 4,
                               avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
# uncomment and configure for extra filtering of MS peak lists
# mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
#                   relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
#                   deIsotopeMS = FALSE, deIsotopeMSMS = FALSE)

# Calculate formula candidates
formulas <- generateFormulas(fGroups, "genform", mslists, relMzDev = 5,
                             adduct = "[M+H]+", elements = "CHNOP",
                             calculateFeatures = TRUE, featThresholdAnn = 0.75)

# Find compound structure candidates
compounds <- generateCompounds(fGroups, mslists, "metfrag", method = "CL", dbRelMzDev = 5,
                               fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "[M+H]+", database = "pubchem", maxCandidatesToStop = 2500)
compounds <- addFormulaScoring(compounds, formulas, TRUE)

# Perform automatic generation of components
components <- generateComponents(fGroups, "ramclustr", ionization = "positive")


# -------------------------
# reporting
# -------------------------


reportCSV(fGroups, path = "report", reportFeatures = FALSE, formulas = formulas,
          compounds = compounds, compoundsNormalizeScores = "max",
          components = components)

reportPDF(fGroups, path = "report", reportFGroups = TRUE, formulas = formulas, reportFormulaSpectra = TRUE,
          compounds = compounds, compoundsNormalizeScores = "max",
          components = components, MSPeakLists = mslists)

reportHTML(fGroups, path = "report", reportPlots = c("chord", "venn", "upset", "eics", "formulas"), formulas = formulas,
           compounds = compounds, compoundsNormalizeScores = "max",
           components = components, MSPeakLists = mslists,
           selfContained = FALSE, openReport = TRUE)
