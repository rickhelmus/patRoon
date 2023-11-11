#' @include main.R
#' @include feature_annotations.R
NULL

#' Compound annotations class
#'
#' Contains data for compound annotations for feature groups.
#'
#' \code{compounds} objects are obtained from \link[=generateCompounds]{compound generators}. This class is derived from
#' the \code{\link{featureAnnotations}} class, please see its documentation for more methods and other details.
#'
#' @param formulas The \code{\link{formulas}} object that should be used for scoring/annotation. For \code{plotSpectrum}
#'   and \code{annotatedPeakList}: set to \code{NULL} to ignore.
#' @param obj,object,compounds,x The \code{compound} object.
#' @param index The numeric index of the candidate structure.
#'
#'   For \code{plotStructure} and \code{getMCS}: multiple indices (\emph{i.e.} vector with length >=2) should be
#'   specified to plot/calculate the most common substructure (MCS). Alternatively, \samp{-1} may be specified to select
#'   all candidates.
#'
#'   For \code{plotSpectrum}: two indices can be specified to compare spectra. In this case \code{groupName} should
#'   specify values for the spectra to compare.
#' @param \dots For \code{plotSpectrum}: Further arguments passed to \code{\link[graphics]{plot}}.
#'
#'   For \code{delete}: passed to the function specified as \code{j}.
#'
#'   for \code{filter}: passed to the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.
#'
#'   For \code{consensus}: any further (and unique) \code{compounds} objects.
#'
#'   \setsPassedArgs1{compounds}
#'
#' @template plotSpec-args
#'
#' @templateVar normParam normalizeScores
#' @templateVar excludeParam excludeNormScores
#' @template norm-args
#'
#' @slot MS2QuantMeta Metadata from \pkg{MS2Quant} filled in by \code{predictRespFactors}.
#' 
#' @note The values ranges in the \code{scoreLimits} slot, which are used for normalization of scores, are based on the
#'   \emph{original} scorings when the compounds were generated (\emph{prior} to employing the \code{topMost} filter to
#'   \code{\link{generateCompounds}}).
#'
#' @seealso The \code{\link{featureAnnotations}} base class for more relevant methods and
#'   \code{\link{generateCompounds}}.
#'
#' @templateVar class compounds
#' @template class-hierarchy
#'
#' @export
compounds <- setClass("compounds",
                      slots = c(MS2QuantMeta = "list"),
                      contains = "featureAnnotations")

setMethod("initialize", "compounds", function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object@groupAnnotations <- lapply(.Object@groupAnnotations, function(ann) ann[, UID := InChIKey1])
    return(.Object)
})


#' @rdname compounds-class
compoundsConsensus <- setClass("compoundsConsensus",
                               slots = c(mergedConsensusNames = "character"),
                               contains = "compounds")
setMethod("mergedConsensusNames", "compoundsConsensus", function(obj) obj@mergedConsensusNames)

#' @describeIn compounds Returns default scorings that are excluded from normalization.
#' @export
setMethod("defaultExclNormScores", "compounds", function(obj) c("score", "individualMoNAScore", "annoTypeCount",
                                                                "annotHitCount", "libMatch"))

setMethod("annScoreNames", "compounds", function(obj, onlyNums)
{
    ret <- callNextMethod()
    if (onlyNums && length(ret) > 0)
    {
        # exclude suspect presence 'scores'
        suspScores <- compoundScorings("metfrag")
        suspScores <- suspScores[suspScores$suspect_list == TRUE, ]
        suspCols <- getAllMergedConsCols(suspScores, ret, mergedConsensusNames(obj))
        ret <- setdiff(ret, suspCols)
    }
    return(ret)
})


#' @describeIn compounds Show summary information for this object.
#' @export
setMethod("show", "compounds", function(object)
{
    callNextMethod()

    mn <- mergedConsensusNames(object, FALSE)
    if (length(mn) > 1)
        printf("Merged: %s\n", paste0(mn, collapse = ", "))

    printf("Number of feature groups with compounds in this object: %d\n", length(annotations(object)))

    cCounts <- if (length(object) == 0) 0 else sapply(annotations(object), nrow)
    printf("Number of compounds: %d (total), %.1f (mean), %d - %d (min - max)\n",
           sum(cCounts), mean(cCounts), min(cCounts), max(cCounts))
})

#' @describeIn compounds Returns a list containing for each feature group a
#'   character vector with database identifiers for all candidate compounds. The
#'   list is named by feature group names, and is typically used with the
#'   \code{identifiers} option of \code{\link{generateCompoundsMetFrag}}.
#' @aliases identifiers
#' @export
setMethod("identifiers", "compounds", function(compounds)
{
    compTable <- annotations(compounds)
    return(sapply(names(compTable),
                  function(grp) unlist(strsplit(as.character(compTable[[grp]]$identifier), ";")), simplify = FALSE))
})

#' @describeIn compounds Provides rule based filtering for generated compounds. Useful to eliminate unlikely candidates
#'   and speed up further processing. Also see the \code{\link[=filter,featureAnnotations-method]{featureAnnotations}}
#'   method.
#'
#' @param minExplainedPeaks,scoreLimits Passed to the
#'   \code{\link[=filter,featureAnnotations-method]{featureAnnotations}} method.
#'
#' @param minScore,minFragScore,minFormulaScore Minimum overall score, in-silico fragmentation score and formula score,
#'   respectively. Set to \code{NULL} to ignore. The \code{scoreLimits} argument allows for more advanced score
#'   filtering.
#'
#' @export
setMethod("filter", "compounds", function(obj, minExplainedPeaks = NULL, minScore = NULL, minFragScore = NULL,
                                          minFormulaScore = NULL, scoreLimits = NULL, ...)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertNumber, . ~ minScore + minFragScore + minFormulaScore, finite = TRUE,
           null.ok = TRUE, fixed = list(add = ac)) # note: negative scores allowed for SIRIUS
    checkmate::assertList(scoreLimits, null.ok = TRUE, types = "numeric", add = ac)
    checkmate::reportAssertions(ac)

    if (!is.null(minScore) || !is.null(minFragScore) || !is.null(minFormulaScore))
    {
        if (is.null(scoreLimits))
            scoreLimits <- list()
        
        minVals <- c(score = minScore, fragScore = minFragScore, formulaScore = minFormulaScore)
        minVals <- minVals[!sapply(minVals, is.null)]
        for (sc in names(minVals))
            scoreLimits[[sc]] <- c(minVals[[sc]], Inf)
    }
    
    return(callNextMethod(obj, minExplainedPeaks, scoreLimits, ...))
})

#' @describeIn compounds Adds formula ranking data from a \code{\link{formulas}}
#'   object as an extra compound candidate scoring (\code{formulaScore} column).
#'   The formula score for each compound candidate is between \samp{0-1}, where
#'   \emph{zero} means no match with any formula candidates, and \emph{one}
#'   means that the compound candidate's formula is the highest ranked.
#'
#' @templateVar scoreName formula score
#' @templateVar scoreWeightName formulaScoreWeight
#' @template update_comp_score-args
#' 
#' @return \code{addFormulaScoring} returns a \code{compounds} object updated
#'   with formula scoring.
#'
#' @aliases addFormulaScoring
#' @export
setMethod("addFormulaScoring", "compounds", function(compounds, formulas, updateScore, formulaScoreWeight)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(formulas, "formulas", add = ac)
    checkmate::assertFlag(updateScore, add = ac)
    checkmate::assertNumber(formulaScoreWeight, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(compounds) == 0)
        return(compounds)

    cTable <- annotations(compounds)
    cGNames <- names(cTable)

    # UNDONE?
    if (length(mergedConsensusNames(compounds, FALSE)) > 0)
        stop("Currently formula scoring cannot be calculated for consensus results. Please add the scorings before calling consensus()")
    
    calculateScores <- function(cr, forms)
    {
        fCount <- nrow(forms)
        # add one to no-match results to make it the worst score
        fRanks <- match(cr$neutral_formula, forms$neutral_formula, nomatch = fCount + 1)
        fRanks <- (fCount - (fRanks - 1)) / fCount # convert to 0-1 where one is best ranked and zero a no match
        
        return(fRanks)
    }

    printf("Adding formula scoring...\n")
    prog <- openProgBar(0, length(cTable))

    cTable <- lapply(seq_along(cTable), function(grpi)
    {
        forms <- formulas[[cGNames[grpi]]]

        ct <- copy(cTable[[grpi]])

        if (nrow(ct) == 0)
            return(ct)

        if (length(forms) == 0 || nrow(forms) == 0)
            ct[, formulaScore := 0]
        else
            ct[, formulaScore := calculateScores(ct, forms)]

        # update overall scoring
        if (updateScore && any(ct$formulaScore > 0))
        {
            normFormScores <- ct$formulaScore / max(ct$formulaScore)
            ct[, score := score + formulaScoreWeight * normFormScores]
        }

        setTxtProgressBar(prog, grpi)
        return(ct)
    })
    names(cTable) <- cGNames

    setTxtProgressBar(prog, length(cTable))
    close(prog)

    compounds@groupAnnotations <- cTable
    
    return(addCompoundScore(compounds, "formulaScore", updateScore, formulaScoreWeight))
})

#' @describeIn compounds Calculates the maximum common substructure (MCS)
#'   for two or more candidate structures for a feature group. This method uses
#'   the \code{\link{get.mcs}} function from \CRANpkg{rcdk}.
#' @return \code{getMCS} returns an \CRANpkg{rcdk} molecule object
#'   (\code{IAtomContainer}).
#' @export
setMethod("getMCS", "compounds", function(obj, index, groupName)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkIntegerish(index, lower = 1, any.missing = FALSE, min.len = 2, unique = TRUE),
        checkmate::checkTRUE(index == -1),
        .var.name = "index"
    )

    assertChoiceSilent(groupName, groupNames(obj), add = ac)
    checkmate::reportAssertions(ac)

    if (length(index) == 1 && index == -1)
        index <- seq_len(nrow(annotations(obj)[[groupName]]))

    mols <- getMoleculesFromSMILES(annotations(obj)[[groupName]][["SMILES"]][index])
    mcons <- mols[[1]]
    if (length(mols) > 1)
    {
        for (i in seq(2, length(mols)))
        {
            if (!isValidMol(mols[[i]]))
                return(emptyMol())

            # might fail if there is no overlap...
            tryCatch(mcons <- rcdk::get.mcs(mcons, mols[[i]]), error = function(e) FALSE)
            if (is.null(mcons) || mcons == FALSE)
                return(emptyMol())
        }
    }

    return(mcons)
})

# NOTE: argument docs 'borrowed' from plotSpec-args.R template

#' @describeIn compounds Plots a structure of a candidate compound using the
#'   \CRANpkg{rcdk} package. If multiple candidates are specified (\emph{i.e.}
#'   by specifying a \code{vector} for \code{index}) then the maximum common
#'   substructure (MCS) of the selected candidates is drawn.
#'
#' @param width,height The dimensions (in pixels) of the raster image that
#'   should be plotted.
#'
#' @references \addCitations{rcdk}{1}
#'
#' @export
setMethod("plotStructure", "compounds", function(obj, index, groupName, width = 500, height = 500)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(index, lower = -1, any.missing = FALSE, min.len = 1, unique = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    aapply(checkmate::assertNumber, . ~ width + height, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    compTable <- annotations(obj)[[groupName]]

    if (is.null(compTable) || nrow(compTable) == 0)
        return(NULL)

    if (length(index) > 1 || index == -1)
        mol <- getMCS(obj, index, groupName)
    else
        mol <- getMoleculesFromSMILES(compTable$SMILES[index], emptyIfFails = TRUE)[[1]]

    img <- getRCDKStructurePlot(mol, width, height)
    plot(img)
})

setMethod("plotStructureHash", "compounds", function(obj, index, groupName, width = 500,
                                                     height = 500)
{
    compTable <- annotations(obj)[[groupName]]
    SMI <- if (is.null(compTable) || nrow(compTable) == 0) NULL else compTable$SMILES[index]
    return(makeHash(SMI, width, height))
})

#' @describeIn compounds Plots a barplot with scoring of a candidate compound.
#'
#' @param onlyUsed If \code{TRUE} then only scorings are plotted that actually
#'   have been used to rank data (see the \code{scoreTypes} argument to
#'   \code{\link{generateCompoundsMetFrag}} for more details).
#'
#' @export
setMethod("plotScores", "compounds", function(obj, index, groupName, normalizeScores = "max",
                                              excludeNormScores = defaultExclNormScores(obj),
                                              onlyUsed = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertChoice(normalizeScores, c("none", "max", "minmax"))
    checkmate::assertCharacter(excludeNormScores, min.chars = 1, null.ok = TRUE, add = ac)
    checkmate::assertFlag(onlyUsed, add = ac)
    checkmate::reportAssertions(ac)

    annTable <- annotations(obj)[[groupName]]
    if (is.null(annTable) || nrow(annTable) == 0 || index > nrow(annTable))
        return(NULL)
    
    mcn <- mergedConsensusNames(obj, FALSE)
    mcnAll <- mergedConsensusNames(obj, TRUE)
    if (length(mcn) == 0)
        mcn <- mcnAll # HACK: only split set scores if not consensus
    
    if (normalizeScores != "none")
        annTable <- normalizeAnnScores(annTable, annScoreNames(obj, TRUE), obj@scoreRanges[[groupName]], mcnAll,
                                       normalizeScores == "minmax", excludeNormScores)
    
    scoreCols <- getAllMergedConsCols(annScoreNames(obj, FALSE), names(annTable), mcnAll)
    if (onlyUsed)
        scoreCols <- intersect(scoreCols, obj@scoreTypes)
    
    makeScoresPlot(annTable[index, scoreCols, with = FALSE], mcn)
})

setMethod("plotScoresHash", "compounds", function(obj, index, groupName, normalizeScores = "max",
                                                  excludeNormScores = defaultExclNormScores(obj),
                                                  onlyUsed = TRUE)
{
    annTable <- annotations(obj)[[groupName]]
    if (is.null(annTable) || nrow(annTable) == 0 || index > nrow(annTable))
        annTable <- NULL
    else if (normalizeScores == "none")
        annTable <- annTable[index]
    
    return(makeHash(index, annTable, normalizeScores, excludeNormScores, onlyUsed))
})

#' @describeIn compounds Returns an MS/MS peak list annotated with data from a
#'   given candidate compound for a feature group.
#'
#' @param onlyAnnotated Set to \code{TRUE} to filter out any peaks that could
#'   not be annotated.
#'
#' @export
setMethod("annotatedPeakList", "compounds", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                     onlyAnnotated = FALSE)
{
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE)
    allFGroups <- groupNames(obj)
    if (!is.null(formulas))
        allFGroups <- union(allFGroups, groupNames(formulas))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    assertChoiceSilent(groupName, allFGroups, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertFlag(onlyAnnotated, add = ac)
    checkmate::reportAssertions(ac)

    if (is.null(formulas) || is.null(formulas[[groupName]]))
        return(doAnnotatePeakList(MSPeakLists[[groupName]][["MSMS"]], annotations(obj)[[groupName]], index,
                                  onlyAnnotated))
    
    # NOTE: onlyAnnotated is set FALSE so we have the complete peaklist for merging
    ret <- doAnnotatePeakList(MSPeakLists[[groupName]][["MSMS"]], annotations(obj)[[groupName]], index,
                              onlyAnnotated = FALSE)
    if (is.null(ret))
        return(ret)
    
    annTable <- annotations(obj)[[groupName]]
    if (!is.null(annTable) && nrow(annTable) >= index)
    {
        formIndex <- which(annTable$neutral_formula[index] == formulas[[groupName]]$neutral_formula)
        if (length(formIndex) != 0)
        {
            # both annotated peak lists should have equal rows etc since onlyAnnotated was set to FALSE, just merge
            # formula columns
            
            annPLForms <- annotatedPeakList(formulas, index = formIndex, groupName = groupName,
                                            MSPeakLists = MSPeakLists, onlyAnnotated = FALSE)

            if (!is.null(annPLForms) && !is.null(annPLForms[["ion_formula"]]))
            {
                if (is.null(ret[["ion_formula"]]))
                {
                    # only formula annotations
                    ret[, c("ion_formula", "annotated") := .(annPLForms$ion_formula, annPLForms$annotated)]
                    ret[annotated == TRUE, mergedBy := algorithm(formulas)]
                }
                else
                {
                    inMe <- ret$annotated
                    inForms <- annPLForms$annotated
                    
                    # UNDONE: handle different formula assignments?
                    set(ret, j = "ion_formula", value = fifelse(!is.na(ret$ion_formula), ret$ion_formula,
                                                                annPLForms$ion_formula))
                    set(ret, j = "mergedBy", value = fcase(inMe & inForms, paste0(algorithm(obj), ",", algorithm(formulas)),
                                                           inMe, algorithm(obj),
                                                           inForms, algorithm(formulas)))
                    set(ret, j = "annotated", value = inMe | inForms)
                }
            }
        }
        else # no formula annotations
            ret[annotated == TRUE, mergedBy := algorithm(obj)]
    }
    
    if (onlyAnnotated)
        ret <- ret[annotated == TRUE]

    return(ret[])
})

#' @describeIn compounds Plots an annotated spectrum for a given candidate compound for a feature group. Two spectra can
#'   be compared by specifying a two-sized vector for the \code{index} and \code{groupName} arguments.
#'
#' @param plotStruct If \code{TRUE} then the candidate structure is drawn in the spectrum. Currently not supported when
#'   comparing spectra.
#' @param title The title of the plot. If \code{NULL} a title will be automatically made.
#' @param maxMolSize Numeric vector of size two with the maximum width/height of the candidate structure (relative to
#'   the plot size).
#' @param molRes Numeric vector of size two with the resolution of the candidate structure (in pixels).
#'
#' @template specSimParams-arg
#' 
#' @template plot-lim
#'
#' @template fsubscript_source
#'
#' @export
setMethod("plotSpectrum", "compounds", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                plotStruct = FALSE, title = NULL, specSimParams = getDefSpecSimParams(),
                                                mincex = 0.9, xlim = NULL, ylim = NULL,
                                                maxMolSize = c(0.2, 0.4), molRes = c(100, 100), ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(index, lower = 1, min.len = 1, max.len = 2, any.missing = FALSE, add = ac)
    checkmate::assertCharacter(groupName, min.len = 1, max.len = 2, min.chars = 1, add = ac)
    if (length(index) != length(groupName))
        stop("Lengths of index and groupName should be equal.")
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE, add = ac)
    checkmate::assertFlag(plotStruct, add = ac)
    checkmate::assertNumber(mincex, lower = 0, finite = TRUE, add = ac)
    assertXYLim(xlim, ylim, add = ac)
    aapply(checkmate::assertNumeric, . ~ maxMolSize + molRes, finite = TRUE, len = 2, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    if (length(groupName) == 1)
    {
        spec <- annotatedPeakList(obj, index, groupName, MSPeakLists, formulas)
        if (is.null(spec))
            return(NULL)

        if (index > nrow(obj[[groupName]]))
            stop(sprintf("Specified candidate index out of range %d/%d", index, nrow(obj[[groupName]])), call. = FALSE)
                
        compr <- obj[[groupName]][index, ]
        mol <- NULL
        if (plotStruct)
        {
            mol <- getMoleculesFromSMILES(compr$SMILES)
            if (!isValidMol(mol))
                mol <- NULL
        }
        
        if (is.null(title))
            title <- getCompoundsSpecPlotTitle(compr$compoundName, compr$neutral_formula)
        
        makeMSPlot(getMSPlotData(spec, 2), mincex, xlim, ylim, main = title, ..., mol = mol,
                   maxMolSize = maxMolSize, molRes = molRes)
    }
    else
    {
        if (plotStruct)
            stop("Cannot plot structure when comparing spectra") # UNDONE?

        for (i in seq_len(2))
        {
            if (is.null(obj[[groupName[i]]]))
                stop(paste("No data for specified feature group:", groupName[i]), call. = FALSE)
            if (index[i] > nrow(obj[[groupName[i]]]))
                stop(sprintf("Specified candidate index out of range %d/%d", index[i], nrow(obj[[groupName[i]]])),
                     call. = FALSE)
        }

        if (is.null(title))
        {
            compr1 <- obj[[groupName[1]]][index[1], ]; compr2 <- obj[[groupName[2]]][index[2], ]
            title <- getCompoundsSpecPlotTitle(compr1$compoundName, compr1$neutral_formula,
                                               compr2$compoundName, compr2$neutral_formula)
        }
        
        binnedPLs <- getBinnedPLPair(MSPeakLists, groupName, NULL, 2, specSimParams, "unique", mustExist = TRUE)
        
        topSpec <- mergeBinnedAndAnnPL(binnedPLs[[1]], annotatedPeakList(obj, index[1], groupName[1], MSPeakLists,
                                                                         formulas), 1)
        
        bottomSpec <- mergeBinnedAndAnnPL(binnedPLs[[2]], annotatedPeakList(obj, index[2], groupName[2], MSPeakLists,
                                                                            formulas), 2)
        plotData <- getMSPlotDataOverlay(list(topSpec, bottomSpec), TRUE, FALSE, 2, "overlap")
        makeMSPlotOverlay(plotData, title, mincex, xlim, ylim, ...)
    }
})

setMethod("plotSpectrumHash", "compounds", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                                    plotStruct = FALSE, title = NULL,
                                                    specSimParams = getDefSpecSimParams(),
                                                    mincex = 0.9, xlim = NULL, ylim = NULL,
                                                    maxMolSize = c(0.2, 0.4), molRes = c(100, 100), ...)
{
    if (length(groupName) > 1)
    {
        # recursive call for both candidates
        args <- list(obj = obj, MSPeakLists = MSPeakLists, formulas = formulas, plotStruct = plotStruct, title = title,
                     specSimParams = specSimParams, mincex = mincex, xlim = xlim, ylim = ylim, maxMolSize = maxMolSize,
                     molRes = molRes, ...)
        return(makeHash(do.call(plotSpectrumHash, c(args, list(index = index[1], groupName = groupName[1]))),
                        do.call(plotSpectrumHash, c(args, list(index = index[2], groupName = groupName[2])))))
    }
    
    compTable <- annotations(obj)[[groupName]]
    cRow <- if (is.null(compTable) || nrow(compTable) == 0) NULL else compTable[index]
    
    fRow <- NULL
    if (!is.null(cRow) && !is.null(formulas) && !is.null(formulas[[groupName]]))
    {
        wh <- match(cRow$neutral_formula, formulas[[groupName]]$neutral_formula)
        fRow <- if (is.na(wh)) NULL else formulas[[groupName]][wh]
    }
    
    return(makeHash(cRow, fRow, getSpec(MSPeakLists, groupName, 2, NULL), plotStruct, title, mincex, xlim, ylim, ...))
})

#' @rdname pred-quant
#' @export
setMethod("predictRespFactors", "compounds", function(obj, fGroups, calibrants, eluent, organicModifier, pHAq,
                                                      concUnit = "ugL", calibConcUnit = concUnit, updateScore = FALSE,
                                                      scoreWeight = 1, parallel = TRUE)
{
    checkPackage("MS2Quant", "kruvelab/MS2Quant")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    assertQuantEluent(eluent, fGroups, add = ac)
    checkmate::assertChoice(organicModifier, c("MeOH", "MeCN"), add = ac)
    checkmate::assertNumber(pHAq, finite = TRUE, add = ac)
    aapply(assertConcUnit, . ~ concUnit + calibConcUnit, fixed = list(add = ac))
    aapply(checkmate::assertFlag, . ~ updateScore + parallel, fixed = list(add = ac))
    checkmate::assertNumber(scoreWeight, finite = TRUE, lower = 1, add = ac)
    checkmate::reportAssertions(ac)
    
    calibrants <- assertAndPrepareQuantCalib(calibrants, calibConcUnit)
    
    if (length(obj) == 0)
        return(obj)
    
    printf("Predicting response factors from SMILES with MS2Quant for %d candidates...\n", length(obj))
    
    MS2QMD <- list()
    
    obj@groupAnnotations <- doApply("Map", parallel, groupNames(obj), annotations(obj), f = function(grp, ann)
    {
        ann <- copy(ann)
        inp <- data.table(group = grp, SMILES = ann$SMILES)
        res <- predictRespFactorsSMILES(inp, groupInfo(fGroups), calibrants, eluent, organicModifier, pHAq, concUnit)
        if (length(MS2QMD) == 0)
            MS2QMD <<- res$MD
        if (!is.null(ann[["RF_SMILES"]]))
            ann[, RF_SMILES := NULL] # clearout for merge below
        ann <- merge(ann, res$RFs[, c("SMILES", "RF_SMILES"), with = FALSE], by = "SMILES", sort = FALSE, all.x = TRUE)
        doProgress()
        return(ann)
    })

    obj@MS2QuantMeta <- MS2QMD
    
    return(addCompoundScore(obj, "RF_SMILES", updateScore, scoreWeight))
})

#' @rdname pred-tox
#' @export
setMethod("predictTox", "compounds", function(obj, LC50Mode = "static", concUnit = "ugL", updateScore = FALSE,
                                              scoreWeight = 1, parallel = TRUE)
{
    checkPackage("MS2Tox", "kruvelab/MS2Tox")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(LC50Mode, c("static", "flow"))
    assertConcUnit(concUnit, add = ac)
    aapply(checkmate::assertFlag, . ~ updateScore + parallel, fixed = list(add = ac))
    checkmate::assertNumber(scoreWeight, finite = TRUE, lower = 1, add = ac)
    checkmate::reportAssertions(ac)
    
    if (length(obj) == 0)
        return(obj)
    
    printf("Predicting LC50 values from SMILES with MS2Tox for %d candidates...\n", length(obj))
    
    obj@groupAnnotations <- doApply("Map", parallel, groupNames(obj), annotations(obj), f = function(grp, ann)
    {
        ann <- copy(ann)
        pr <- predictLC50SMILES(ann$SMILES, LC50Mode, concUnit)
        if (!is.null(ann[["LC50_SMILES"]]))
            ann[, LC50_SMILES := NULL] # clearout for merge below
        ann <- merge(ann, pr, by = "SMILES", sort = FALSE, all.x = TRUE)
        doProgress()
        return(ann)
    })
    
    return(addCompoundScore(obj, "LC50_SMILES", updateScore, scoreWeight))
})


setMethod("prepareConsensusLabels", "compounds", function(obj, ..., labels)
{
    if (is.null(labels))
        labels <- sapply(list(obj, ...), algorithm)
    
    if (anyDuplicated(labels))
    {
        allCompounds <- list(obj, ...)
        # duplicate algorithms used, try to form unique names by adding library
        labels <- sapply(allCompounds, function(cmp)
        {
            db <- if (length(cmp) == 0) "empty" else cmp[[1]]$database[1]
            return(paste0(substr(algorithm(cmp), 1, 3), "_", substr(db, 1, 3)))
        })
    }
    
    return(callNextMethod(obj, ..., labels = labels))
})

#' @templateVar what compounds
#' @template consensus-form_comp
#'
#' @templateVar what compounds
#' @template consensus-common-args
#'
#' @param labels A \code{character} with names to use for labelling. If \code{NULL} labels are automatically generated.
#'
#' @return \code{consensus} returns a \code{compounds} object that is produced by merging multiple specified
#'   \code{compounds} objects.
#'
#' @export
setMethod("consensus", "compounds", function(obj, ..., MSPeakLists,
                                             specSimParams = getDefSpecSimParams(removePrecursor = TRUE),
                                             absMinAbundance = NULL, relMinAbundance = NULL, uniqueFrom = NULL,
                                             uniqueOuter = FALSE, rankWeights = 1, labels = NULL)
{
    # NOTE: keep args in sync with compoundsSet method
    
    allCompounds <- c(list(obj), list(...))

    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allCompounds, types = "compounds", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    assertSpecSimParams(specSimParams, add = ac)
    checkmate::assertNumeric(rankWeights, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCharacter(labels, min.chars = 1, len = length(allCompounds), null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    labels <- prepareConsensusLabels(obj, ..., labels = labels)
    
    assertConsCommonArgs(absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, labels)

    uniqueCols <- c("neutral_formula", "SMILES", "InChI", "InChIKey1", "InChIKey2", "InChIKey", "neutralMass",
                    "RF_SMILES", "RF_SIRFP", "LC50_SMILES", "LC50_SIRFP")
    cons <- doFeatAnnConsensus(obj, ..., rankWeights = rankWeights, annNames = labels,
                               uniqueCols = uniqueCols)
    
    # rename & merge score types and ranges
    scoreTypes <- Reduce(union, mapply(allCompounds, labels, FUN = function(cmp, cn)
    {
        fifelse(cmp@scoreTypes %chin% uniqueCols, cmp@scoreTypes, paste0(cmp@scoreTypes, "-", cn))
    }))

    scRanges <- Reduce(modifyList, Map(allCompounds, labels, f = function(cmp, cn)
    {
        lapply(cmp@scoreRanges, function(scrg)
        {
            n <- names(scrg)
            wh <- !n %chin% uniqueCols
            names(scrg)[wh] <- paste0(n[wh], "-", cn)
            return(scrg)
        })
    }))

    ret <- compoundsConsensus(groupAnnotations = cons, scoreTypes = scoreTypes, scoreRanges = scRanges,
                              algorithm = paste0(unique(sapply(allCompounds, algorithm)), collapse = ","),
                              mergedConsensusNames = labels, MSPeakLists = MSPeakLists, specSimParams = specSimParams)
    
    ret <- filterFeatAnnConsensus(ret, absMinAbundance, relMinAbundance, uniqueFrom, uniqueOuter, FALSE)
    
    return(ret)
})

#' Automatic compound annotation
#'
#' Automatically perform chemical compound annotation for feature groups.
#'
#' Several algorithms are provided to automatically perform compound annotation for feature groups. To this end,
#' measured masses for all feature groups are searched within online database(s) (\emph{e.g.}
#' \href{https://pubchem.ncbi.nlm.nih.gov/}{PubChem}) to retrieve a list of potential candidate chemical compounds.
#' Depending on the algorithm and its parameters, further scoring of candidates is then performed using, for instance,
#' matching of measured and theoretical isotopic patterns, presence within other data sources such as patent databases
#' and similarity of measured and in-silico predicted MS/MS fragments. Note that this process is often quite time
#' consuming, especially for large feature group sets. Therefore, this is often one of the last steps within the
#' workflow and not performed before feature groups have been prioritized.
#'
#' @templateVar func generateCompounds
#' @templateVar what generate compounds
#' @templateVar ex1 generateCompoundsMetFrag
#' @templateVar ex2 generateCompoundsSIRIUS
#' @templateVar algos metfrag,sirius,library
#' @templateVar algosSuffix MetFrag,SIRIUS,Library
#' @templateVar ret compounds
#' @template generic-algo
#'
#' @param fGroups \code{\link{featureGroups}} object which should be annotated. This should be the same or a subset of
#'   the object that was used to create the specified \code{MSPeakLists}. In the case of a subset only the remaining
#'   feature groups in the subset are considered.
#' @param MSPeakLists A \code{\link{MSPeakLists}} object that was generated for the supplied \code{fGroups}.
#' @param \dots Any parameters to be passed to the selected compound generation algorithm.
#'
#' @section Scorings: Each algorithm implements their own scoring system. Their names have been simplified and
#'   harmonized where possible. The \code{\link{compoundScorings}} function can be used to get an overview of both the
#'   algorithm specific and generic scoring names.
#'
#' @templateVar UID first-block \acronym{InChIKey}
#' @template featAnnSets-gen
#' 
#' @return A \code{\link{compounds}} derived object containing all compound annotations.
#'
#' @templateVar what generateCompounds
#' @template main-rd-method
#' @export
setMethod("generateCompounds", "featureGroups", function(fGroups, MSPeakLists, algorithm, specSimParams, ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertChoice(algorithm, c("metfrag", "sirius", "library"), add = ac)
    assertSpecSimParams(specSimParams)
    checkmate::reportAssertions(ac)
    
    f <- switch(algorithm,
                metfrag = generateCompoundsMetFrag,
                sirius = generateCompoundsSIRIUS,
                library = generateCompoundsLibrary)

    f(fGroups, MSPeakLists, specSimParams, ...)
})
