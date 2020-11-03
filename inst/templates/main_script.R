## Script automatically generated on {{ date() }}

library(patRoon)


{{ header("initialization") }}


workPath <- "{{ destination }}"
setwd(workPath)
{{ optionalCodeBlock(generateAnaInfo == "table") }}

# Load analysis table
anaInfo <- read.csv("{{ analysisTableFile }}")
{{ endCodeBlock() }}
{{ optionalCodeBlock(generateAnaInfo == "script") }}

anaInfo <- generateAnalysisInfo(paths = c({{ paste0("\"", unique(analyses$path), "\"", collapse = ", ") }}),
                                groups = c({{ paste0("\"", analyses$group, "\"", collapse = ", ") }}),
                                blanks = c({{ paste0("\"", analyses$blank, "\"", collapse = ", ") }}))
{{ endCodeBlock() }}
{{ optionalCodeBlock(generateAnaInfo == "example") }}

# Take example data from patRoonData package (triplicate solvent blank + triplicate standard)
anaInfo <- generateAnalysisInfo(paths = patRoonData::exampleDataPath(),
                                groups = c(rep("solvent", 3), rep("standard", 3)),
                                blanks = "solvent")
{{ endCodeBlock() }}
{{ optionalCodeBlock(generateAnaInfo == "none") }}

# NOTE: please set anaInfo to a valid data.frame with analysis information. See ?`analysis-information` for more details.
anaInfo <- data.frame(path = character(), analysis = character(), group = character(), blank = character())
{{ endCodeBlock() }}
{{ optionalCodeBlock(preTreatOpts$do) }}

# Set to FALSE to skip data pre-treatment
doDataPretreatment <- TRUE
if (doDataPretreatment)
{
    {{ endCodeBlock() }}
    setDAMethod(anaInfo, "{{ preTreatOpts$DAMethod }}") {{ optionalLine(nzchar(preTreatOpts$DAMethod)) }}
    recalibrarateDAFiles(anaInfo) {{ optionalLine(preTreatOpts$doDACalib) }}
    {{ optionalCodeBlock(nzchar(preTreatOpts$convAlgo) && "mzML" %in% preTreatOpts$convTo) }}
    convertMSFiles(anaInfo = anaInfo, from = {{ paste0('"', preTreatOpts$convFrom, '"') }},
                   to = "mzML", algorithm = "{{ preTreatOpts$convAlgo }}", centroid = {{ preTreatOpts$centroid }})
    {{ endCodeBlock() }}
    {{ optionalCodeBlock(nzchar(preTreatOpts$convAlgo) && "mzXML" %in% preTreatOpts$convTo) }}
    convertMSFiles(anaInfo = anaInfo, from = {{ paste0('"', preTreatOpts$convFrom, '"') }},
                   to = "mzXML", algorithm = "{{ preTreatOpts$convAlgo }}", centroid = {{ preTreatOpts$centroid }})
    {{ endCodeBlock() }}
} {{ optionalLine(preTreatOpts$do) }}


{{ header("features") }}


# Find all features.
{{ optionalCodeBlock(featFinderOpts$algo == "OpenMS") }}
# NOTE: see manual for many more options
fList <- findFeatures(anaInfo, "openms")
{{ endCodeBlock() }}
{{ optionalCodeBlock(featFinderOpts$algo == "XCMS") }}
# NOTE: see XCMS manual for many more options
fList <- findFeatures(anaInfo, "xcms3", param = xcms::CentWaveParam())
{{ endCodeBlock() }}
{{ optionalCodeBlock(featFinderOpts$algo == "enviPick") }}
# NOTE: see enviPickWrap manual for many more options
fList <- findFeatures(anaInfo, "envipick")
{{ endCodeBlock() }}
{{ optionalCodeBlock(featFinderOpts$algo == "Bruker") }}
fList <- findFeatures(anaInfo, "bruker", doFMF = "auto")
{{ endCodeBlock() }}

# Group and align features between analysis
{{ optionalCodeBlock(featGrouperOpts$algo == "OpenMS") }}
fGroups <- groupFeatures(fList, "openms")
{{ endCodeBlock() }}
{{ optionalCodeBlock(featGrouperOpts$algo == "XCMS") }}
fGroups <- groupFeatures(fList, "xcms3", rtalign = TRUE,
                         groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(fList)$group),
                         retAlignParam = xcms::ObiwarpParam())
{{ endCodeBlock() }}

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = {{ filterFGroupsOpts$preIntThr }}, absMinIntensity = {{ filterFGroupsOpts$intThr }},
                  relMinReplicateAbundance = {{ filterFGroupsOpts$repAbundance }}, maxReplicateIntRSD = {{ filterFGroupsOpts$maxRepRSD }},
                  blankThreshold = {{ filterFGroupsOpts$blankThr }}, removeBlanks = {{ filterFGroupsOpts$removeBlanks }},
                  retentionRange = {{ if (is.null(filterFGroupsOpts$retRange)) "NULL" else paste0("c(", paste0(filterFGroupsOpts$retRange, collapse = ", "), ")") }}, mzRange = {{ if (is.null(filterFGroupsOpts$mzRange)) "NULL" else paste0("c(", paste0(filterFGroupsOpts$mzRange, collapse = ", "), ")") }})
{{ optionalCodeBlock(nzchar(suspectList)) }}

# Filter feature groups by suspects
suspFile <- read.csv("{{ suspectList }}", stringsAsFactors = FALSE)
scr <- screenSuspects(fGroups, suspFile, rtWindow = 12, mzWindow = 0.005,
                      adduct = {{ if (!nzchar(suspectAdduct)) "NULL" else paste0("\"", suspectAdduct, "\"") }})
fGroups <- groupFeaturesScreening(fGroups, scr)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind || formulaOpts$algo != "" || identOpts$algo != "" || componentOpts$algo != "") }}


{{ header("annotation") }}


{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind) }}
# Retrieve MS peak lists
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind && peakListOpts$algo == "mzR") }}
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 5, precursorMzWindow = {{ if (precursorMzWindow == 0) "NULL" else precursorMzWindow }},
                               avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind && peakListOpts$algo == "Bruker" && featFinderOpts$algo != "Bruker") }}
mslists <- generateMSPeakLists(fGroups, "bruker", maxMSRtWindow = 5, bgsubtr = TRUE, MSMSType = "MSMS", avgFGroupParams = avgPListParams)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind && peakListOpts$algo == "Bruker" && featFinderOpts$algo == "Bruker") }}
mslists <- generateMSPeakLists(fGroups, "brukerfmf", avgFGroupParams = avgPListParams)
{{ endCodeBlock() }}
{{ optionalCodeBlock(doMSPeakFind) }}
# uncomment and configure for extra filtering of MS peak lists
# mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
#                  relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
#                  deIsotopeMS = FALSE, deIsotopeMSMS = FALSE)
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo != "") }}

# Calculate formula candidates
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo == "GenForm") }}
formulas <- generateFormulas(fGroups, "genform", mslists, relMzDev = 5,
                             adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}", elements = "CHNOP",
                             calculateFeatures = TRUE, featThreshold = 0.75)
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo == "Bruker") }}
formulas <- generateFormulas(fGroups, "bruker", precursorMzSearchWindow = 0.002, featThreshold = 0.75,
                             adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}")
{{ endCodeBlock() }}
{{ optionalCodeBlock(formulaOpts$algo == "SIRIUS") }}
formulas <- generateFormulas(fGroups, "sirius", mslists, relMzDev = 5,
                             adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}", elements = "CHNOP",
                             profile = "qtof", calculateFeatures = TRUE, featThreshold = 0.75)
{{ endCodeBlock() }}
{{ optionalCodeBlock(identOpts$algo != "") }}

# Find compound structure candidates
{{ endCodeBlock() }}
{{ optionalCodeBlock(identOpts$algo == "MetFrag") }}
compounds <- generateCompounds(fGroups, mslists, "metfrag", method = "CL", dbRelMzDev = 5,
                               fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "{{ if (polarity == 'positive') '[M+H]+' else '[M-H]-' }}", database = "pubchem", maxCandidatesToStop = 2500)
compounds <- addFormulaScoring(compounds, formulas, TRUE) {{ optionalLine(formulaOpts$algo != "") }}
{{ endCodeBlock() }}
{{ optionalCodeBlock(identOpts$algo == "SIRIUS") }}
compounds <- generateCompounds(fGroups, mslists, "sirius", relMzDev = 5,
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
                                 mzRange = c(5, 120), elements = c("C", "H", "O"), rtDev = 30,
                                 absMzDev = 0.002)
{{ endCodeBlock() }}
{{ optionalCodeBlock(nzchar(suspectList) && annotateSus) }}

# Annotate suspects
{{ endCodeBlock() }}
{{ optionalCodeBlock(nzchar(suspectList) && annotateSus && genIDLevelFile) }}
IDLevelRules <- read.csv("idlevelrules.csv", stringsAsFactors = FALSE)
fGroups <- annotateSuspects(fGroups, MSPeakLists = mslists, formulas = formulas,
                            compounds = compounds, IDLevelRules = IDLevelRules)
{{ endCodeBlock() }}
{{ optionalCodeBlock(nzchar(suspectList) && annotateSus && !genIDLevelFile) }}
fGroups <- annotateSuspects(fGroups, MSPeakLists = mslists, formulas = formulas,
                            compounds = compounds)
{{ endCodeBlock() }}
{{ optionalCodeBlock(length(reportFormats) > 0) }}


{{ header("reporting") }}


{{ endCodeBlock() }}
{{ optionalCodeBlock("CSV" %in% reportFormats) }}
reportCSV(fGroups, path = "report", reportFeatures = FALSE, formulas = {{ if (formulaOpts$algo != "") "formulas" else "NULL" }},
          compounds = {{ if (identOpts$algo != "") "compounds" else "NULL" }}, compoundsNormalizeScores = "max",
          components = {{ if (componentOpts$algo != "") "components" else "NULL" }})

{{ endCodeBlock() }}
{{ optionalCodeBlock("PDF" %in% reportFormats) }}
reportPDF(fGroups, path = "report", reportFGroups = TRUE, formulas = {{ if (formulaOpts$algo != "") "formulas" else "NULL" }}, reportFormulaSpectra = TRUE,
          compounds = {{ if (identOpts$algo != "") "compounds" else "NULL" }}, compoundsNormalizeScores = "max",
          components = {{ if (componentOpts$algo != "") "components" else "NULL" }}, MSPeakLists = {{ if (formulaOpts$algo != "" || identOpts$algo != "") "mslists" else "NULL" }})

{{ endCodeBlock() }}
{{ optionalCodeBlock("HTML" %in% reportFormats) }}
reportHTML(fGroups, path = "report", reportPlots = c("chord", "venn", "upset", "eics", "formulas"),
           formulas = {{ if (formulaOpts$algo != "") "formulas" else "NULL" }}, compounds = {{ if (identOpts$algo != "") "compounds" else "NULL" }}, compoundsNormalizeScores = "max",
           components = {{ if (componentOpts$algo != "") "components" else "NULL" }}, MSPeakLists = {{ if (formulaOpts$algo != "" || identOpts$algo != "") "mslists" else "NULL" }},
           selfContained = FALSE, openReport = TRUE)

{{ endCodeBlock() }}
