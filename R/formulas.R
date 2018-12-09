#' @include main.R
NULL

#' Formula lists class
#'
#' Contains data of generated chemical formulae for given feature groups.
#'
#' \code{formulas} objects are obtained from \link[=formula-generation]{formula
#' generators}.
#'
#' @slot formulas Lists of all generated formulae. Use the \code{formulaTable}
#'   method for access.
#' @slot algorithm The algorithm that was used for generation of formulae. Use
#'   the \code{algorithm} method for access.
#'
#' @param obj,x,object,formulas The \code{formulas} object.
#' 
#' @templateVar seli analyses
#' @templateVar selOrderi analyses()
#' @templateVar selj feature groups
#' @templateVar selOrderj groupNames()
#' @template sub_op-args
#'
#' @seealso \code{\link{formulaConsensus}}
#' 
#' @export
formulas <- setClass("formulas",
                     slots = c(formulas = "list", groupFormulas = "list",  algorithm = "character"))

#' @describeIn formulas Accessor method to obtain generated formulae.
#' @return \code{formulas} returns a \code{list} containing for each analysis
#'   and each feature group a \code{\link{data.table}} with an overview of all
#'   generated formulae and other data such as candidate scoring and MS/MS
#'   fragments.
#' @export
setMethod("formulaTable", "formulas", function(obj) obj@formulas)

# UNDONE
setMethod("groupFormulas", "formulas", function(obj) obj@groupFormulas)

#' @describeIn formulas Accessor method for the algorithm (a character
#'   string) used to generate formulae.
#' @export
setMethod("algorithm", "formulas", function(obj) obj@algorithm)

#' @templateVar class formulas
#' @templateVar what analyses
#' @template strmethod
#' @export
setMethod("analyses", "formulas", function(obj) names(obj@formulas))

#' @templateVar class formulas
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "formulas", function(obj) unique(unlist(sapply(obj@formulas, names, simplify = FALSE), use.names = FALSE)))

# UNDONE: this seems not so useful now, change to unique formulae?
#' @describeIn formulas Obtain total number of formulae entries.
#' @export
setMethod("length", "formulas", function(x) sum(unlist(recursiveApplyDT(x@formulas, nrow))))

#' @describeIn formulas Show summary information for this object.
#' @export
setMethod("show", "formulas", function(object)
{
    printf("A formulas object (%s)\n", class(object))
    printf("Algorithm: %s\n", algorithm(object))

    printf("Total formula count: %d\n", length(object))

    if (length(object) == 0)
    {
        ma <- mfg <- 0
    }
    else
    {
        ft <- formulaTable(object)
        ft <- ft[lengths(ft) > 0]
        formCounts <- lapply(ft, function(fa) sapply(fa, nrow))
        ma <- mean(sapply(formCounts, function(fc) sum(fc)))
        mfg <- mean(sapply(formCounts, function(fc) mean(fc)))
    }
    
    printf("Average formulas per analysis: %.1f\n", ma)
    printf("Average formulas per feature group: %.1f\n", mfg)

    showObjectSize(object)
})

#' @describeIn formulas Subset on analyses/feature groups.
#' @export
setMethod("[", c("formulas", "ANY", "ANY", "missing"), function(x, i, j, ...)
{
    if (!missing(i))
        assertSubsetArg(i)
    if (!missing(j))
        assertSubsetArg(j)
    
    # non-existing indices result in NULL values --> prune
    
    if (!missing(i))
    {
        if (!is.character(i))
            i <- analyses(x)[i]
        x@formulas <- pruneList(x@formulas[i])
    }
    
    if (!missing(j))
    {
        if (!is.character(j))
            j <- groupNames(x)[j]
        x@formulas <- sapply(x@formulas, function(a)
        {
            ret <- a[j]
            return(pruneList(ret))
        }, simplify = FALSE)
        x@formulas <- pruneList(x@formulas, TRUE)
    }
    
    return(x)
})

#' @describeIn formulas Extract a formula table from a specified analysis/feature group.
#' @export
setMethod("[[", c("formulas", "ANY", "ANY"), function(x, i, j)
{
    assertExtractArg(i)
    assertExtractArg(j)
    
    if (!is.character(i))
        i <- analyses(x)[i]
    if (!is.character(j))
        j <- groupNames(x)[j]
    
    return(x@formulas[[c(i, j)]])
})

#' @describeIn formulas Generates a consensus and other summarizing data for
#'   formulae of feature groups that are present in multiple analyses.
#'
#' @param \dots Any further \code{formulas} objects on top of the object given
#'   for the \code{obj} parameter to which a consensus should be generated.
#'   Using the \dots parameter formulae from various
#'   \link[=formula-generation]{formula generators} can be merged.
#' @param fGroups The \code{\link{featureGroups}} object that was used when
#'   formulae were generated.
#' @param formAnaThreshold,formListThreshold Fractional minimum amount (0-1) to
#'   which a generated formula should be present in all the analyses
#'   (\code{formAnaThreshold}) or in all given formula lists
#'   (\code{formListThreshold}) containing the related feature group. For
#'   instance, a value of \samp{0.5} for \code{formAnaThreshold} means that a
#'   particular formula should be present in at least \samp{50\%} of all
#'   analyses containing the feature group that was used to generate the
#'   formula.
#' @param maxFormulas,maxFragFormulas Maximum amount of candidate formulae (or
#'   fragment formulae) per feature group.
#' @param minIntensity Results from features below this intensity will be
#'   excluded.
#' @param maxIntensity Results from features above this intensity will be
#'   excluded.
#' @param minPreferredFormulas Preferred minimum amount of candidate formulae
#'   per feature group. See \code{minPreferredIntensity}.
#' @param minPreferredIntensity Minimum preferred intensity of formulae
#'   candidates per feature group. Results from features below this intensity
#'   will be removed unless the total count of formulae of the related feature
#'   group is below \code{minPreferredFormulas}. When a higher
#'   \code{minIntensity} is given this parameter has no effect.
#' @param elements An optional character vector containing elements for which
#'   the total count present in each formula should be calculated. For instance,
#'   \code{c("C", "H")} will count all carbon and hydrogen atoms present in each
#'   formula.
#' @param fragElements Same as \code{elements} parameter, but for MS/MS
#'   formulae.
#'
#' @return \code{consensus} returns a \code{\link{formulaConsensus}} object.
#'
#' @export
setMethod("consensus", "formulas", function(obj, ..., formThreshold = 0)
{
    allFormulas <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allFormulas, types = "formulas", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertNumber(formThreshold, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    allFormulas <- allFormulas[lengths(allFormulas) > 0]
    if (length(allFormulas) < 2)
        stop("Need at least two non-empty formulas objects")
    
    allFormNames <- sapply(allFormulas, algorithm)
    # UNDONE?
    # if (length(allFormulas) > length(unique(allFormNames)))
    #     stop("Consensus can only be generated from different algorithms at this moment.")
    allFormNames <- make.unique(allFormNames)

    allFormulasLists <- sapply(seq_along(allFormulas), function(fi)
    {
        return(lapply(groupFormulas(allFormulas[[fi]]), function(ft)
        {
            ret <- copy(ft)
            setnames(ret, paste0(names(ret), "-", allFormNames[fi]))
            return(ret)
        }))
        
    }, simplify = FALSE)
    
    # UNDONE: remove old style columns?
    uniqueCols <- c("neutral_formula", "formula_mz", "error", "dbe", "frag_formula", "frag_mz",
                    "frag_formula_mz", "frag_error", "neutral_loss", "frag_dbe", "min_intensity", "max_intensity",
                    "ana_min_intensity", "ana_max_intensity")
    
    consFormulaList <- allFormulasLists[[1]]
    leftName <- allFormNames[1]
    for (righti in seq(2, length(allFormulasLists)))
    {
        rightName <- allFormNames[righti]
        
        printf("Merging %s with %s... ", paste0(allFormNames[seq_len(righti-1)], collapse = ","), rightName)
        
        rightTable <- allFormulasLists[[righti]]
        
        for (grp in union(names(consFormulaList), names(rightTable)))
        {
            if (is.null(rightTable[[grp]]))
                next # nothing to merge
            else if (is.null(consFormulaList[[grp]])) # not yet present
            {
                mTable <- rightTable[[grp]]
                
                # rename columns that should be unique from right to left
                unColsPresent <- uniqueCols[sapply(uniqueCols, function(uc) !is.null(mTable[[paste0(uc, "-", rightName)]]))]
                setnames(mTable, paste0(unColsPresent, "-", rightName), paste0(unColsPresent, "-", leftName))
            }
            else
            {
                haveMSMS <- paste0("frag_formula-", leftName) %in% names(consFormulaList[[grp]]) &&
                            paste0("frag_formula-", rightName) %in% names(consFormulaList[[grp]])
                
                mergeCols <- "formula"
                if (haveMSMS)
                    mergeCols <- c(mergeCols, "byMSMS", "frag_formula")
                mTable <- merge(consFormulaList[[grp]], rightTable[[grp]], all = TRUE,
                                by.x = paste0(mergeCols, "-", leftName),
                                by.y = paste0(mergeCols, "-", rightName))
                
                # remove duplicate columns that shouldn't
                for (col in uniqueCols)
                {
                    colLeft <- paste0(col, "-", leftName)
                    colRight <- paste0(col, "-", rightName)
                    if (!is.null(mTable[[colRight]]))
                    {
                        mTable[, (colLeft) := ifelse(!is.na(get(colLeft)), get(colLeft), get(colRight))]
                        mTable[, (colRight) := NULL]
                    }
                }
            }
            
            consFormulaList[[grp]] <- mTable
        }
        
        cat("Done!\n")
    }
    
    
    printf("Determining coverage and final scores... ")
    
    # Determine coverage of compounds between objects and the merged score. The formula column can be
    # used for the former as there is guaranteed to be one for each merged object.
    for (grpi in seq_along(consFormulaList))
    {
        # fix up de-duplicated column names
        deDupCols <- c(uniqueCols)
        leftCols <- paste0(deDupCols, "-", leftName)
        deDupCols <- deDupCols[leftCols %in% names(consFormulaList[[grpi]])]
        leftCols <- leftCols[leftCols %in% names(consFormulaList[[grpi]])]
        if (length(leftCols) > 0)
            setnames(consFormulaList[[grpi]], leftCols, deDupCols)
        
        formCols <- grep("formula-", colnames(consFormulaList[[grpi]]), value = TRUE)
        
        if (length(formCols) == 0) # nothing was merged
            consFormulaList[[grpi]][, coverage := 1 / length(allCompounds)]
        else
        {
            for (r in seq_len(nrow(consFormulaList[[grpi]])))
                set(consFormulaList[[grpi]], r, "coverage",
                    sum(sapply(formCols, function(c) !is.na(consFormulaList[[grpi]][[c]][r]))) / length(allCompounds))
        }
        
        if (formThreshold > 0)
            consFormulaList[[grpi]] <- consFormulaList[[grpi]][coverage >= formThreshold]
        
        setcolorder(consFormulaList[[grpi]], formConsensusColOrder(consFormulaList[[grpi]]))
    }
    
    cat("Done!")
        
    return(formulas(formulas = list(), groupFormulas = consFormulaList,
                    algorithm = paste0(unique(sapply(allFormulas, algorithm)), collapse = ",")))
    
    
    
    
    allFGroups <- unique(unlist(sapply(allFormulasLists, names)))
    formCons <- sapply(allFGroups, function(grp)
    {
        ftPresent <- sapply(allFormulasLists, function(fl) !is.null(fl[[grp]]))
        formTables <- sapply(allFormulasLists[ftPresent], "[[", grp, simplify = FALSE)
        ftNames <- allFormNames[ftPresent]
        
        if (length(formTables) < 2)
            return(formTables[[1]])
        
        consTable <- formTables[[1]]
        for (otherFTi in seq(2, length(formTables)))
        {
            algCons <- paste0(ftNames[seq_len(otherFTi-1)], collapse = ",")
            printf("Merging %s and %s... ", algCons, ftNames[otherFTi])
            
            otherFT <- formTables[[otherFTi]]
            
            haveMSMS <- "frag_formula" %in% colnames(consTable) && "frag_formula" %in% colnames(otherFT)
            
            mergeCols <- "formula"
            if (haveMSMS)
                mergeCols <- c(mergeCols, "byMSMS", "frag_formula")
            
            mTable <- merge(consTable, otherFT, all = TRUE, by = mergeCols)
            
            # score and anaCoverage columns may now be double and should be renamed
            # note that double columns get .x/.y suffix on merge
            for (col in c("score", "anaCoverage"))
            {
                cols <- paste0(col, c(".x", ".y"))
                if (!any(sapply(cols, function(cl) is.null(mTable[[cl]]))))
                    setnames(mTable, cols, paste0(col, "-", c(algCons, ftNames[otherFTi])))
            }
            
            # Other (potentially) double columns which should be (more or less) the same, so only keep one of these
            # UNDONE: remove old style columns?
            takeOneCols <- c("ret", "mz", "neutral_formula", "formula_mz", "error", "dbe", "frag_formula", "frag_mz",
                             "frag_formula_mz", "frag_error", "neutral_loss", "frag_dbe", "min_intensity", "max_intensity",
                             "ana_min_intensity", "ana_max_intensity")
            for (col in takeOneCols)
            {
                colLeft <- paste0(col, ".x")
                colRight <- paste0(col, ".y")
                if (colLeft %in% colnames(mTable))
                {
                    mTable[, (col) := ifelse(!is.na(get(colLeft)), get(colLeft), get(colRight))]
                    mTable[, c(colLeft, colRight) := NULL]
                }
            }
            
            cat("Done!\n")
            
            # return(mTable)
            consTable <- mTable
        }
        
        # Calculate coverage of formulae across formula lists. score is a
        # column that is present with all algorithms, so we can use presence of
        # its value as counter
        scoreCols <- grep("score-", colnames(consTable), value = TRUE)
        for (r in seq_len(nrow(consTable)))
            set(consTable, r, "mergeCoverage", sum(sapply(scoreCols, function(c) !is.na(consTable[[c]][r]))) / length(allFormulas))
        
        if (formThreshold > 0)
            consTable <- consTable[mergeCoverage >= formThreshold]
        
        setcolorder(consTable, formConsensusColOrder(consTable))
    })
    
    return(formulas(formulas = list(), groupFormulas = formCons,
                    algorithm = paste0(unique(sapply(allFormulas, algorithm)), collapse = ",")))
})

if (FALSE){
setMethod("consensus", "formulas", function(obj, ..., fGroups, formAnaThreshold = 0.75, formListThreshold = 0,
                                            maxFormulas = 10, maxFragFormulas = 10, minIntensity = NULL,
                                            maxIntensity = NULL, minPreferredFormulas = 1,
                                            minPreferredIntensity = NULL, elements = c(), fragElements = c())
{
    allFLists <- c(list(obj), list(...))
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allFLists, types = "formulas", min.len = 1, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    aapply(checkmate::assertNumber, . ~ formAnaThreshold + formListThreshold, lower = 0, finite = TRUE,
           fixed = list(add = ac))
    aapply(checkmate::assertCount, . ~ maxFormulas + maxFragFormulas + minPreferredFormulas, positive = TRUE,
           null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minIntensity + maxIntensity + minPreferredIntensity,
           lower = 0, finite = TRUE, null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertCharacter, . ~ elements + fragElements, min.chars = 1,
           any.missing = FALSE, null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)
    
    ftind <- groupFeatIndex(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)
    gTable <- groups(fGroups)
    
    
    if (length(allFLists) > length(unique(sapply(allFLists, algorithm))))
        stop("Consensus can only be generated from different algorithms at this moment.")
    
    allFLists <- allFLists[lengths(allFLists) > 0]
    if (length(allFLists) < 1)
        stop("Need at least one non-empty formulas objects")
    
    formConsList <- lapply(allFLists,
                           function (fl) consensusForFormulaList(fl, fGroups, formAnaThreshold, maxFormulas, maxFragFormulas,
                                                                 minIntensity, maxIntensity, minPreferredFormulas,
                                                                 minPreferredIntensity))
    
    if (length(formConsList) == 1)
        retConsTable <- formConsList[[1]]@formulas
    else
    {
        retConsTable <- Reduce(function(fCons1, fCons2)
        {
            printf("Merging %s and %s... ", algorithm(fCons1), algorithm(fCons2))
            
            forms1 <- formulaTable(fCons1)
            forms2 <- formulaTable(fCons2)
            
            haveMSMS <- "frag_formula" %in% colnames(forms1) && "frag_formula" %in% colnames(forms2)
            
            mergeCols <- c("group", "formula")
            if (haveMSMS)
                mergeCols <- c(mergeCols, "byMSMS", "frag_formula")
            
            mTable <- merge(forms1, forms2, all = TRUE, by = mergeCols)
            
            # score and anaCoverage columns may now be double and should be renamed
            # note that double columns get .x/.y suffix on merge
            for (col in c("score", "anaCoverage"))
            {
                cols <- paste0(col, c(".x", ".y"))
                if (!any(sapply(cols, function(cl) is.null(mTable[[cl]]))))
                    setnames(mTable, cols, paste0(col, "-", c(algorithm(fCons1), algorithm(fCons2))))
            }
            
            # Other (potentially) double columns which should be (more or less) the same, so only keep one of these
            takeOneCols <- c("ret", "mz", "neutral_formula", "formula_mz", "error", "dbe", "frag_formula", "frag_mz",
                             "frag_formula_mz", "frag_error", "neutral_loss", "frag_dbe", "min_intensity", "max_intensity",
                             "ana_min_intensity", "ana_max_intensity")
            for (col in takeOneCols)
            {
                colLeft <- paste0(col, ".x")
                colRight <- paste0(col, ".y")
                if (colLeft %in% colnames(mTable))
                {
                    mTable[, (col) := ifelse(!is.na(get(colLeft)), get(colLeft), get(colRight))]
                    mTable[, c(colLeft, colRight) := NULL]
                }
            }
            
            cat("Done!\n")
            
            return(mTable)
        }, formConsList)
        
        # Calculate coverage of formulae across formula lists. anaCoverage is a
        # column that is present with all algorithms, so we can use presence of
        # its value as counter
        anaCovCols <- grep("anaCoverage-", colnames(retConsTable), value = TRUE)
        for (r in seq_len(nrow(retConsTable)))
            set(retConsTable, r, "listCoverage", sum(sapply(anaCovCols, function(c) !is.na(retConsTable[[c]][r]))) / length(formConsList))
        
        if (formListThreshold > 0)
            retConsTable <- retConsTable[listCoverage >= formListThreshold]
    }
    
    setcolorder(retConsTable, formConsensusColOrder(retConsTable))
    
    if (length(elements) > 0)
    {
        # Retrieve element lists from formulas
        el <- getElements(retConsTable$formula, elements)
        retConsTable <- cbind(retConsTable, el)
    }
    if (!is.null(retConsTable[["frag_formula"]]) && length(fragElements) > 0)
    {
        el <- getElements(retConsTable$frag_formula, fragElements)
        colnames(el) <- sapply(colnames(el), function(e) paste0("frag_", e), USE.NAMES = F)
        retConsTable <- cbind(retConsTable, el)
    }
    
    return(formulaConsensus(formulas = retConsTable,
                            algorithm = paste0(unique(sapply(formConsList, algorithm)), collapse = ",")))
})}

#' Formula consensus class
#'
#' Contains all formulae for given feature groups after a consensus was made.
#'
#' Objects of \code{formulaConsensus} are generated by the
#' \code{\link{consensus}} method.
#'
#' @slot formulas A \code{\link{data.table}} containing an overview of all
#'   consensus formulae. Use the \code{formulaTable} method for access.
#' @slot algorithm The algorithm that was used for generation of formulae. Use
#'   the \code{algorithm} method for access.
#'
#' @param obj,x,object The \code{formulaConsensus} object.
#' 
#' @templateVar seli feature groups
#' @templateVar selOrderi groupNames()
#' @templateVar dollarOpName feature group
#' @template sub_op-args
#' 
#' @seealso \code{\link{formulas}}
#' 
#' @export
formulaConsensus <- setClass("formulaConsensus",
                             slots = c(formulas = "data.table", algorithm = "character"),
                             prototype = list(formulas = data.table(), algorithm = "none"))

#' @describeIn formulaConsensus Accessor method to obtain generated formulae.
#' @return \code{formulaTable} returns a \code{\link{data.table}} containing all
#'   generated formulae, scoring information, MS/MS fragments, etc.
#' @export
setMethod("formulaTable", "formulaConsensus", function(obj) obj@formulas)

#' @describeIn formulaConsensus Accessor method for the algorithm (a character
#'   string) used to generate formulae.
#' @export
setMethod("algorithm", "formulaConsensus", function(obj) obj@algorithm)

#' @templateVar class formulaConsensus
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "formulaConsensus", function(obj) unique(obj@formulas$group))

#' @describeIn formulaConsensus Obtain total number of formulae entries.
#' @export
setMethod("length", "formulaConsensus", function(x) nrow(x@formulas))

#' @describeIn formulaConsensus Show summary information for this object.
#' @export
setMethod("show", "formulaConsensus", function(object)
{
    printf("A formulaConsensus object (%s)\n", class(object))
    printf("Algorithm: %s\n", paste0(algorithm(object), collapse = "/"))

    printf("Total formula count: %d (all), %d (MS) and %d (MS/MS)\n", length(object),
           nrow(object@formulas[byMSMS == FALSE]), nrow(object@formulas[byMSMS == TRUE]))

    fCounts <- object@formulas[, .(byMSMS, .N), by = group]
    printf("Average formulas per feature group: %.1f (all), %.1f (MS), %.1f (MS/MS)\n", mean(fCounts$N),
           mean(fCounts[byMSMS == FALSE, N]), mean(fCounts[byMSMS == TRUE, N]))

    showObjectSize(object)
})

#' @describeIn formulaConsensus Subset on feature groups.
#' @export
setMethod("[", c("formulaConsensus", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
        assertSubsetArg(i)
    
    if (!missing(i))
    {
        if (!is.character(i))
            i <- groupNames(x)[i]
        x@formulas <- x@formulas[group %in% i]
    }
    
    return(x)
})

#' @describeIn formulaConsensus Extract formula information a feature group.
#' @export
setMethod("[[", c("formulaConsensus", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    if (!is.character(i))
        i <- groupNames(x)[i]
    return(x@formulas[group == i])
})

#' @describeIn formulaConsensus Extract formula information a feature group.
#' @export
setMethod("$", "formulaConsensus", function(x, name)
{
    eval(substitute(x@formulas[group == NAME_ARG], list(NAME_ARG = name)))
})

#' @describeIn formulaConsensus Plots an annotated spectrum for a given
#'   candidate formula of a feature group.
#'
#' @param precursor The formula of the precursor (in ionic form, \emph{i.e.} as
#'   detected by the MS).
#'
#' @template plotSpec-args
#' 
#' @template useGGplot2
#'
#' @return \code{plotSpec} will return a \code{\link[=ggplot2]{ggplot object}}
#'   if \code{useGGPlot2} is \code{TRUE}.
#'
#' @export
setMethod("plotSpec", "formulaConsensus", function(obj, precursor, groupName, MSPeakLists, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertString(precursor, min.chars = 1, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)
    
    formTable <- formulaTable(obj)[group == groupName & byMSMS == TRUE & formula == precursor]

    if (nrow(formTable) == 0)
        return(NULL)

    spec <- MSPeakLists[[groupName]][["MSMS"]]
    if (is.null(spec))
        return(NULL)

    fi <- getFragmentInfoFromForms(spec, formTable)

    if (useGGPlot2)
        return(makeMSPlotGG(spec, fi) + ggtitle(precursor))

    makeMSPlot(spec, fi, main = precursor)
})

#' @templateVar func generateFormulas
#' @templateVar what generate formulae
#' @templateVar ex1 generateFormulasDA
#' @templateVar ex2 generateFormulasGenForm
#' @templateVar algos bruker,genform,sirius
#' @template generic-algo
#'
#' @param ... Any parameters to be passed to the selected formula generation
#'   algorithm.
#'
#' @rdname formula-generation
#' @aliases generateFormulas
#' @export
setMethod("generateFormulas", "featureGroups", function(fGroups, algorithm, ...)
{
    f <- switch(algorithm,
                bruker = generateFormulasDA,
                genform = generateFormulasGenForm,
                sirius = generateFormulasSirius,
                stop("Invalid algorithm! Should be: bruker, genform or sirius"))

    f(fGroups, ...)
})
