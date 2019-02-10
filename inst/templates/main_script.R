## Script automatically generated on {{ date() }}

library(patRoon)

workPath <- "{{ destination }}"
setwd(workPath)
{{ optionalCodeBlock(generateAnaInfo == "table") }}

# Load analysis table
anaInfo <- read.csv("{{ analysisTableFile }}", stringsAsFactors = FALSE, colClasses = "character")
{{ endCodeBlock() }}
{{ optionalCodeBlock(generateAnaInfo == "script") }}

anaInfo <- generateAnalysisInfo(paths = c({{ paste0("\"", unique(analyses$path), "\"", collapse = ", ") }}),
                                groups = c({{ paste0("\"", analyses$group, "\"", collapse = ", ") }}),
                                refs = c({{ paste0("\"", analyses$ref, "\"", collapse = ", ") }}))
{{ endCodeBlock() }}
{{ optionalCodeBlock(dataPretreatmentOpts$DAMethod != "" || length(dataPretreatmentOpts$steps) > 0) }}

# Set to FALSE to skip data pretreatment (e.g. calibration, export, ...)
doDataPretreatment <- TRUE
if (doDataPretreatment)
{
{{ endCodeBlock() }}
    setDAMethod(anaInfo, "{{ dataPretreatmentOpts$DAMethod }}") {{ optionalLine(dataPretreatmentOpts$DAMethod != "") }}
    recalibrarateDAFiles(anaInfo) {{ optionalLine("recalibrate" %in% dataPretreatmentOpts$steps) }}
    exportDAFiles(anaInfo, format = "mzML") {{ optionalLine("expMzML" %in% dataPretreatmentOpts$steps) }}
    exportDAFiles(anaInfo, format = "mzXML") {{ optionalLine("expMzXML" %in% dataPretreatmentOpts$steps) }}
} {{ optionalLine(dataPretreatmentOpts$DAMethod != "" || length(dataPretreatmentOpts$steps) > 0) }}

# Find all features.
{{ optionalCodeBlock(featFinderOpts$algo == "OpenMS") }}
# NOTE: see manual for many more options
fList <- findFeatures(anaInfo, "openms")
{{ endCodeBlock() }}
{{ optionalCodeBlock(featFinderOpts$algo == "XCMS") }}
# NOTE: see XCMS manual for many more options
fList <- featureFinder(anaInfo, "xcms", method = "centWave")
{{ endCodeBlock() }}
{{ optionalCodeBlock(featFinderOpts$algo == "enviPick") }}
# NOTE: see enviPickWrap manual for many more options
fList <- featureFinder(anaInfo, "envipick")
{{ endCodeBlock() }}
{{ optionalCodeBlock(featFinderOpts$algo == "Bruker") }}
fList <- featureFinder(anaInfo, "bruker", doFMF = "auto")
{{ endCodeBlock() }}

# Group and align features between analysis
{{ optionalCodeBlock(featGrouperOpts$algo == "OpenMS") }}
fGroups <- groupFeatures(fList, "openms")
{{ endCodeBlock() }}
{{ optionalCodeBlock(featGrouperOpts$algo == "XCMS") }}
fGroups <- groupFeatures(fList, "xcms", rtalign = TRUE, retcorArgs = list(method = "obiwarp"))
{{ endCodeBlock() }}

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = {{ filterFGroupsOpts$preIntThr }}, absMinIntensity = {{ filterFGroupsOpts$intThr }},
                  relMinReplicateAbundance = {{ filterFGroupsOpts$repAbundance }}, maxReplicateIntRSD = {{ filterFGroupsOpts$maxRepRSD }},
                  blankThreshold = {{ filterFGroupsOpts$blankThr }}, removeBlanks = {{ filterFGroupsOpts$removeBlanks }},
                  retentionRange = {{ if (is.null(filterFGroupsOpts$retRange)) "NULL" else paste0("c(", paste0(filterFGroupsOpts$retRange, collapse = ", "), ")") }}, mzRange = {{ if (is.null(filterFGroupsOpts$mzRange)) "NULL" else paste0("c(", paste0(filterFGroupsOpts$mzRange, collapse = ", "), ")") }})
{{ optionalCodeBlock(doMSPeakFind) }}

# Retrieve MS peak lists
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind && peakListOpts$algo == "mzR") }}
# NOTE: please check all arguments, especially precursorMzWindow!
plists <- generateMSPeakLists(fGroups, "mzr", maxRtMSWidth = 20, precursorMzWindow = {{ if (precursorMzWindow == 0) "NULL" else precursorMzWindow }},
                              avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind && peakListOpts$algo == "Bruker" && featFinderOpts$algo != "Bruker") }}
plists <- generateMSPeakLists(fGroups, "bruker", bgsubtr = TRUE, MSMSType = "MSMS", avgFGroupParams = avgPListParams)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind && peakListOpts$algo == "Bruker" && featFinderOpts$algo == "Bruker") }}
plists <- generateMSPeakLists(fGroups, "brukerfmf", avgFGroupParams = avgPListParams)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind) }}
# uncomment and configure for extra filtering of MS peak lists
# plists <- filter(plists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
#                  relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
#                  deIsotopeMS = FALSE, deIsotopeMSMS = FALSE)
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo != "") }}

# Calculate all formulas
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo == "GenForm") }}
formulas <- generateFormulas(fGroups, "genform", plists, maxMzDev = 5,
                             adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}", elements = "CHNOP",
                             calculateFeatures = TRUE, featThreshold = 0.75)
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo == "Bruker") }}
formulas <- generateFormulas(fGroups, "bruker", precursorMzSearchWindow = 0.002, featThreshold = 0.75)
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo == "SIRIUS") }}
formulas <- generateFormulas(fGroups, "sirius", plists, maxMzDev = 5,
                             adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}", elements = "CHNOP",
                             profile = "qtof", calculateFeatures = TRUE, featThreshold = 0.75)
{{ endCodeBlock() }}
{{ optionalCodeBlock(identOpts$algo != "") }}

# Perform automatic compound identification
{{ endCodeBlock() }}
{{ optionalCodeBlock(identOpts$algo == "MetFrag") }}
compounds <- generateCompounds(fGroups, plists, "metfrag", method = "CL", dbRelMzDev = 5,
                               fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}", database = "pubchem", maxCandidatesToStop = 2500)
compounds <- addFormulaScoring(compounds, formulas, TRUE) {{ optionalLine(formulaOpts$algo != "") }}
{{ endCodeBlock() }}
{{ optionalCodeBlock(identOpts$algo == "SIRIUS") }}
compounds <- generateCompounds(fGroups, plists, "sirius", maxMzDev = 5,
                               adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}", elements = "CHNOP", profile = "qtof",
                               fingerIDDatabase = "pubchem")
{{ endCodeBlock() }}
{{ optionalCodeBlock(componentOpts$algo != "") }}

# Perform automatic generation of components
{{ endCodeBlock() }}
{{ optionalCodeBlock(componentOpts$algo == "RAMClustR") }}
components <- generateComponents(fGroups, "ramclustr", ionization = "{{ polarity }}")
{{ endCodeBlock() }}
{{ optionalCodeBlock(componentOpts$algo == "CAMERA") }}
components <- generateComponents(fGroups, "camera", ionization = "{{ polarity }}")
{{ endCodeBlock() }}
{{ optionalCodeBlock(componentOpts$algo == "nontarget") }}
components <- generateComponents(fGroups, "nontarget", ionization = "{{ polarity }}", rtRange = c(-120, 120),
                                 mzRange = c(5, 120), elements = c("C", "H", "O"), maxRTDev = 30,
                                 maxMzDev = 0.002)
{{ endCodeBlock() }}
{{ optionalCodeBlock(length(reportFormats) > 0) }}

# Report & export results
{{ endCodeBlock() }}
{{ optionalCodeBlock("CSV" %in% reportFormats) }}
reportCSV(fGroups, path = "report", reportFeatures = FALSE, formulas = {{ if (formulaOpts$algo != "") "formulas" else "NULL" }},
          compounds = {{ if (identOpts$algo != "") "compounds" else "NULL" }}, compoundNormalizeScores = "max",
          components = {{ if (identOpts$algo != "") "components" else "NULL" }})

{{ endCodeBlock() }}
{{ optionalCodeBlock("PDF" %in% reportFormats) }}
reportPDF(fGroups, path = "report", reportFGroups = TRUE, formulas = {{ if (formulaOpts$algo != "") "formulas" else "NULL" }}, reportFormulaSpectra = TRUE,
          compounds = {{ if (identOpts$algo != "") "compounds" else "NULL" }}, compoundNormalizeScores = "max",
          components = {{ if (identOpts$algo != "") "components" else "NULL" }}, MSPeakLists = {{ if (formulaOpts$algo != "" || identOpts$algo != "") "plists" else "NULL" }})

{{ endCodeBlock() }}
{{ optionalCodeBlock("MD" %in% reportFormats) }}
reportMD(fGroups, path = "report", reportPlots = c("chord", "venn", "upset", "eics", "formulas"), formulas = {{ if (formulaOpts$algo != "") "formulas" else "NULL" }},
         compounds = {{ if (identOpts$algo != "") "compounds" else "NULL" }}, compoundNormalizeScores = "max",
         components = {{ if (componentOpts$algo != "") "components" else "NULL" }}, MSPeakLists = {{ if (formulaOpts$algo != "" || identOpts$algo != "") "plists" else "NULL" }},
         selfContained = FALSE, openReport = TRUE)

{{ endCodeBlock() }}
