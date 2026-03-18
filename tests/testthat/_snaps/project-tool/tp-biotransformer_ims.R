
library(patRoon)

# -------------------------
# initialization
# -------------------------

workPath <- "<WORK_PATH>/test-np/tp-biotransformer_ims"
setwd(workPath)

# NOTE: please set to a valid data.frame with analysis information. See ?`analysis-information` for more details.
anaInfo <- data.frame(path_centroid = character(), analysis = character(), replicate = character(), blank = character())

# -------------------------
# features
# -------------------------

# Find all features
# NOTE: see the reference manual for many more options
genEICParams <- getPiekEICParams(filter = "none", filterIMS = "none", mzRange = c(80, 800), mzStep = 0.02,
                                 mobRange = c(0.4, 1.3), mobStep = 0.08)
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

# -------------------------
# transformation products
# -------------------------

# Load parent suspect list
suspListParents <- read.csv("suspects")

# Obtain TPs
TPs <- generateTPs("biotransformer", parents = suspListParents, type = "env", generations = 2,
                   TPStructParams = getDefTPStructParams())

# Screen TPs
suspListTPs <- convertToSuspects(TPs, includeParents = FALSE)
suspListTPs <- assignMobilities(suspListTPs, from = "c3sdb", adducts = c("[M+H]+", "[M-H]-", NA), overwrite = FALSE)
fGroups <- screenSuspects(fGroups, suspListTPs, onlyHits = TRUE)

# -------------------------
# annotation
# -------------------------

# Retrieve MS peak lists
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroups, avgFeatParams = avgMSListParams, avgFGroupParams = avgMSListParams)
# Rule based filtering of MS peak lists. You may want to tweak this. See the manual for more information.
mslists <- filter(mslists, MSLevel = 2, absMinIntensity = NULL, relMinIntensity = 0.05, topMostPeaks = 25,
                  maxMZOverPrec = 4)

convertToMFDB(TPs, "TP-database.csv", includeParents = TRUE)
# Calculate compound structure candidates
compounds <- generateCompounds(fGroups, mslists, "metfrag", adduct = "[M+H]+", database = "csv",
                               extraOpts = list(LocalDatabasePath = "TP-database.csv"), maxCandidatesToStop = 2500)

compounds <- estimateIDConfidence(compounds, MSPeakLists = mslists, formulas = NULL, IDFile = "idlevelrules.yml")

# -------------------------
# Parent and TP linkage
# -------------------------

# You probably want to prioritize the data before componentization. Please see the handbook for more info.
componentsTP <- generateComponents(fGroups, "tp", fGroupsTPs = fGroups, TPs = TPs, MSPeakLists = mslists,
                                   compounds = compounds)

# You may want to configure the filtering step below. See the manuals for more details.
componentsTP <- filter(componentsTP, retDirMatch = FALSE, minSpecSim = NULL, minSpecSimPrec = NULL,
                       minSpecSimBoth = NULL, minFragMatches = NULL, minNLMatches = NULL, formulas = NULL)

# Only keep linked feature groups
fGroups <- fGroups[results = componentsTP]

# -------------------------
# reporting
# -------------------------

# Advanced report settings can be edited in the report.yml file.
report(fGroups, MSPeakLists = mslists, formulas = NULL, compounds = compounds, components = componentsTP, TPs = TPs,
       settingsFile = "report.yml", openReport = TRUE)
