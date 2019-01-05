#' @include main.R
#' @include workflow-step.R
NULL

#' Compound lists class
#'
#' Contains data of generated chemical compounds for given feature groups.
#'
#' \code{compounds} objects are obtained from
#' \link[=compound-generation]{compound generators}.
#'
#' @slot compounds Lists of all generated compounds. Use the \code{compounds}
#'   method for access.
#' @slot algorithm The algorithm that was used for generation of compounds. Use
#'   the \code{algorithm} method for access.
#' @param formulas The \code{\link{formulas}} object that should be
#'   used for scoring/annotation. For \code{plotSpec}: set to \code{NULL} to
#'   ignore.
#'
#' @param obj,object,x,compounds The \code{compound} object.
#' @param index The numeric index of the candidate structure. Multiple indices
#'   (\emph{i.e.} vector with length >=2) may be specified for
#'   \code{plotStructure} and are mandatory for \code{getMCS}. Alternatively,
#'   \samp{-1} may be specified to these methods to select all candidates. When
#'   multiple indices are specified for \code{plotStructure}, their maximum
#'   common substructure will be drawn.
#'
#' @templateVar seli feature groups
#' @templateVar selOrderi groupNames()
#' @templateVar dollarOpName feature group
#' @template sub_op-args
#'
#' @return \code{plotSpec} and \code{plotStructure} will return a
#'   \code{\link[=ggplot2]{ggplot object}} if \code{useGGPlot2} is \code{TRUE}.
#'
#' @template plotSpec-args
#'
#' @template useGGplot2
#'
#' @templateVar class compounds
#' @template class-hierarchy
#'
#' @export
compounds <- setClass("compounds",
                      slots = c(compounds = "list"),
                      contains = "workflowStep")

#' @rdname compounds-class
compoundsConsensus <- setClass("compoundsConsensus",
                               slots = c(mergedCompNames = "character"),
                               contains = "compounds")


#' @describeIn compounds Accessor method to obtain generated compounds.
#' @return \code{compoundTable} returns a \code{list} containing for each feature
#'   group a \code{\link{data.table}} with an overview of all candidate
#'   compounds and other data such as candidate scoring, matched MS/MS
#'   fragments, etc.
#' @aliases compoundTable
#' @export
setMethod("compoundTable", "compounds", function(obj) obj@compounds)

#' @describeIn compounds Accessor method for the algorithm (a character
#'   string) used to generate compounds.
#' @export
setMethod("algorithm", "compounds", function(obj) obj@algorithm)

#' @templateVar class compounds
#' @templateVar what feature groups
#' @template strmethod
#' @export
setMethod("groupNames", "compounds", function(obj) names(obj@compounds))

#' @describeIn compounds Obtain total number of candidate compounds.
#' @export
setMethod("length", "compounds", function(x) if (length(x@compounds) > 0) sum(sapply(x@compounds, nrow)) else 0)

#' @describeIn compounds Show summary information for this object.
#' @export
setMethod("show", "compounds", function(object)
{
    callNextMethod()

    mn <- mergedCompoundNames(object)
    if (length(mn) > 1)
        printf("Merged: %s\n", paste0(mn, collapse = ", "))

    printf("Number of feature groups with compounds in this object: %d\n", length(object@compounds))

    cCounts <- if (length(object) == 0) 0 else sapply(object@compounds, nrow)
    printf("Number of compounds: %d (total), %.1f (mean), %d - %d (min - max)\n",
           sum(cCounts), mean(cCounts), min(cCounts), max(cCounts))
})

#' @describeIn compounds Subset on feature groups.
#' @export
setMethod("[", c("compounds", "ANY", "missing", "missing"), function(x, i, ...)
{
    if (!missing(i))
    {
        assertSubsetArg(i)

        if (!is.character(i))
            i <- groupNames(x)[i]

        i <- i[i %in% groupNames(x)]
        x@compounds <- x@compounds[i]
    }

    return(x)
})

#' @describeIn compounds Extract a compound table for a feature group.
#' @export
setMethod("[[", c("compounds", "ANY", "missing"), function(x, i, j)
{
    assertExtractArg(i)
    return(x@compounds[[i]])
})

#' @describeIn compounds Extract a compound table for a feature group.
#' @export
setMethod("$", "compounds", function(x, name)
{
    eval(substitute(x@compounds$NAME_ARG, list(NAME_ARG = name)))
})

#' @describeIn compounds Returns all MS peak list data in a table.
#'
#' @param fragments If \code{TRUE} then information on annotated fragments will
#'   be included.
#' 
#' @template as_data_table-args
#'
#' @export
setMethod("as.data.table", "compounds", function(x, fGroups = NULL, fragments = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", null.ok = TRUE, add = ac)
    checkmate::assertFlag(fragments, add = ac)
    checkmate::reportAssertions(ac)

    if (fragments)
    {
        ret <- rbindlist(lapply(compoundTable(x), function(ct)
        {
            ct <- copy(ct)
            ct[, row := seq_len(nrow(ct))]

            fragTab <- rbindlist(ct$fragInfo, idcol = "row", fill = TRUE)
            fragTab[, PLIndex := NULL]
            cnames <- setdiff(names(fragTab), "row")
            setnames(fragTab, cnames, paste0("frag_", cnames))

            return(merge(ct, fragTab, by = "row")[, -"row"])
        }), idcol = "group", fill = TRUE)
    }
    else
        ret <- rbindlist(compoundTable(x), idcol = "group", fill = TRUE)

    if (!is.null(fGroups))
    {
        ret[, c("ret", "group_mz") := groupInfo(fGroups)[group, c("rts", "mzs")]]
        setcolorder(ret, c("group", "ret", "group_mz"))
    }

    return(ret[, -"fragInfo"])
})

#' @describeIn compounds Returns a list containing for each feature group a
#'   character vector with database identifiers for all candidate compounds. The
#'   list is named by feature group names, and is typically used with the
#'   \code{identifiers} option of \code{\link{generateCompoundsMetfrag}}.
#' @aliases identifiers
#' @export
setMethod("identifiers", "compounds", function(compounds)
{
    compTable <- compoundTable(compounds)
    return(sapply(names(compTable),
                  function(grp) unlist(strsplit(as.character(compTable[[grp]]$identifier), ";")), simplify = FALSE))
})

#' @describeIn compounds Provides rule based filtering for generated compounds.
#'   Useful to eliminate unlikely candidates and speed up further processing.
#'
#' @param minExplainedPeaks,minScore,minFragScore,minFormulaScore Minimal number
#'   of explained peaks, overall score, in-silico fragmentation score and
#'   formula score, respectively. Set to \code{NULL} to ignore. The
#'   \code{scoreLimits} argument allows for more advanced score filtering.
#' @param scoreLimits Filter results by their scores. Should be a named
#'   \code{list} that contains two-sized numeric vectors with the
#'   minimum/maximum value of a score (use \code{-Inf}/\code{Inf} for no
#'   limits). The names of each element should follow the values returned by
#'   \code{\link{compoundScorings}()$name}. For instance,
#'   \code{scoreLimits=list(numberPatents=c(10, Inf))} specifies that
#'   \code{numberPatents} should be at least \samp{10}. For more details of
#'   scorings see \code{\link{compoundScorings}}. Note that a result without a
#'   specified scoring is never removed. Set to \code{NULL} to skip this filter.
#' @param topMost Only keep a maximum of \code{topMost} candidates with highest
#'   score. Set to \code{NULL} to ignore.
#'
#' @templateVar withLoss FALSE
#' @template element-args
#'
#' @return \code{filter} returns a filtered \code{compounds} object.
#'
#' @export
setMethod("filter", "compounds", function(obj, minExplainedPeaks = NULL, minScore = NULL, minFragScore = NULL,
                                          minFormulaScore = NULL, scoreLimits = NULL, elements = NULL,
                                          fragElements = NULL, topMost = NULL)
{
    ac <- checkmate::makeAssertCollection()
    aapply(checkmate::assertCount, . ~ minExplainedPeaks + topMost, positive = c(FALSE, TRUE),
           null.ok = TRUE, fixed = list(add = ac))
    aapply(checkmate::assertNumber, . ~ minScore + minFragScore + minFormulaScore, finite = TRUE,
           null.ok = TRUE, fixed = list(add = ac)) # note: negative scores allowed for SIRIUS
    checkmate::assertList(scoreLimits, null.ok = TRUE, types = "numeric", add = ac)
    if (!is.null(scoreLimits))
    {
        scCols <- unique(compoundScorings()$name)
        checkmate::assertNames(names(scoreLimits), type = "unique", subset.of = scCols, add = ac)
        checkmate::qassertr(scoreLimits, "N2")
    }
    aapply(checkmate::assertCharacter, . ~ elements + fragElements,
           min.chars = 1, min.len = 1, null.ok = TRUE, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    cat("Filtering compounds... ")

    mCompNames <- mergedCompoundNames(obj)
    filterMinCols <- function(cmpTable, col, minVal)
    {
        cols <- getAllCompCols(col, names(cmpTable), mCompNames)
        for (cl in cols)
            cmpTable <- cmpTable[is.na(get(cl)) | get(cl) >= minVal]
        return(cmpTable)
    }

    fMinVals <- c(explainedPeaks = minExplainedPeaks, score = minScore, fragScore = minFragScore,
                  formulaScore = minFormulaScore)
    fMinVals <- fMinVals[!sapply(fMinVals, is.null)]

    oldn <- length(obj)
    obj@compounds <- sapply(obj@compounds, function(cmpTable)
    {
        for (cl in names(fMinVals))
            cmpTable <- filterMinCols(cmpTable, cl, fMinVals[cl])

        if (!is.null(scoreLimits))
        {
            for (sc in names(scoreLimits))
            {
                cols <- getAllCompCols(sc, names(cmpTable), mCompNames)
                if (length(cols) == 0)
                    next
                cmpTable <- cmpTable[cmpTable[, do.call(pmin, c(.SD, list(na.rm = TRUE))) >= scoreLimits[[sc]][1] &
                                                  do.call(pmax, c(.SD, list(na.rm = TRUE))) <= scoreLimits[[sc]][2],
                                              .SDcols = cols]]
            }
        }
        
        if (!is.null(elements))
            cmpTable <- cmpTable[sapply(formula, checkFormula, elements)]
        if (!is.null(fragElements))
        {
            keep <- sapply(cmpTable$fragInfo, function(fi)
            {
                if (nrow(fi) == 0)
                    return(FALSE)
                if (!is.null(fragElements) && !any(sapply(fi$formula, checkFormula, fragElements)))
                    return(FALSE)
                return(TRUE)
            })
            cmpTable <- cmpTable[keep]
        }
        
        if (!is.null(topMost) && nrow(cmpTable) > topMost)
            cmpTable <- cmpTable[seq_len(topMost)]

        return(cmpTable)
    }, simplify = FALSE)

    if (length(obj) > 0)
        obj@compounds <- obj@compounds[sapply(obj@compounds, function(cm) !is.null(cm) && nrow(cm) > 0)]

    newn <- length(obj)
    printf("Done! Filtered %d (%.2f%%) compounds. Remaining: %d\n", oldn - newn, if (oldn == 0) 0 else (1-(newn/oldn))*100, newn)
    return(obj)
})

#' @describeIn compounds Provides compound scoring data that is based on the
#'   presence of candidate formulae being present in a given
#'   \code{\link{formulas}} object. Matched precursor formulae yield one
#'   point, whereas matched MS/MS fragments yield 0.5 point each.
#'
#' @param updateScore If set to \code{TRUE} then the \code{score} column is
#'   updated by adding the normalized \option{formulaScore} (weighted by
#'   \option{formulaScoreWeight}). Currently, this \strong{only} makes sense for
#'   \command{MetFrag} results!
#' @param formulaScoreWeight Weight used to update scoring (see
#'   \code{updateScore} parameter).
#'
#' @return \code{addFormulaScoring} returns a \code{compounds} object updated
#'   with formula scoring.
#'
#' @aliases addFormulaScoring
#' @export
setMethod("addFormulaScoring", "compounds", function(compounds, formulas, updateScore,
                                                     formulaScoreWeight)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(formulas, "formulas", add = ac)
    checkmate::assertFlag(updateScore, add = ac)
    checkmate::assertNumber(formulaScoreWeight, lower = 0, finite = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(compounds) == 0)
        return(compounds)

    cTable <- compoundTable(compounds)
    cGNames <- names(cTable)

    # one point for formula presence, half a point for each fragment
    calculateScores <- function(cr, forms)
    {
        formScores <- as.integer(cr$formula %in% forms$neutral_formula)

        fragScores <- sapply(seq_len(nrow(cr)), function(i)
        {
            sforms <- forms[neutral_formula == cr[i, formula] & byMSMS == TRUE]
            if (nrow(sforms) == 0)
                return(0)

            return(sum(cr$fragInfo[[i]]$formula %in% sforms$frag_formula) * 0.5)
        })

        return(formScores + fragScores)
    }

    printf("Adding formula scoring...\n")
    prog <- txtProgressBar(0, length(cTable), style = 3)

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
        if (any(ct$formulaScore > 0))
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

    compounds@compounds <- cTable

    return(compounds)
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

    assertChoiceSilent(groupName, names(obj@compounds), add = ac)
    checkmate::reportAssertions(ac)

    if (length(index) == 1 && index == -1)
        index <- seq_len(nrow(compoundTable(obj)[[groupName]]))

    mols <- getMoleculesFromSMILES(compoundTable(obj)[[groupName]][["SMILES"]][index])
    mcons <- mols[[1]]
    if (length(mols) > 1)
    {
        for (i in seq(2, length(mols)))
        {
            if (!isValidMol(mols[[i]]))
                return(emptyMol())

            # might fail if there is no overlap...
            tryCatch(mcons <- rcdk::get.mcs(mcons, mols[[i]]), error = function(e) FALSE)
            if (mcons == FALSE)
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
setMethod("plotStructure", "compounds", function(obj, index, groupName, width = 500, height = 500, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertIntegerish(index, lower = -1, any.missing = FALSE, min.len = 1, unique = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    aapply(checkmate::assertNumber, . ~ width + height, lower = 0, finite = TRUE, fixed = list(add = ac))
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)

    compTable <- compoundTable(obj)[[groupName]]

    if (is.null(compTable) || nrow(compTable) == 0)
        return(NULL)

    if (length(index) > 1 || index == -1)
        mol <- getMCS(obj, index, groupName)
    else
        mol <- getMoleculesFromSMILES(compTable$SMILES[index], emptyIfFails = TRUE)[[1]]

    if (useGGPlot2)
    {
        raster <- rcdk::view.image.2d(mol, rcdk::get.depictor(width, height))
        img <- magick::image_trim(magick::image_read(raster))
        cowplot::ggdraw() + cowplot::draw_image(img)
    }
    else
        rcdkplot(mol, width, height)
})

#' @describeIn compounds Plots a barplot with scoring of a candidate compound.
#'
#' @templateVar normParam normalizeScores
#' @templateVar excludeParam excludeNormScores
#' @template comp_norm
#'
#' @aliases plotScores
#' @export
setMethod("plotScores", "compounds", function(obj, index, groupName, normalizeScores, excludeNormScores, useGGPlot2)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertChoice(normalizeScores, c("none", "max", "minmax"))
    checkmate::assertFlag(useGGPlot2, add = ac)
    checkmate::reportAssertions(ac)

    compTable <- compoundTable(obj)[[groupName]]

    if (is.null(compTable) || nrow(compTable) == 0)
        return(NULL)

    mcn <- mergedCompoundNames(obj)

    if (normalizeScores != "none")
        compTable <- normalizeCompScores(compTable, mcn, normalizeScores == "minmax", excludeNormScores)

    scoreCols <- getAllCompCols(c(getCompScoreColNames(), getCompSuspectListColNames()), names(compTable), mcn)
    scores <- setnames(transpose(compTable[index, scoreCols, with = FALSE]), "score")
    scores[, type := scoreCols]
    scores <- scores[!is.na(score)]

    if (length(mcn) > 1)
    {
        scores[, merged := "both"]
        for (n in mcn)
        {
            withM <- which(grepl(paste0("-", n), scores[["type"]], fixed = TRUE))
            set(scores, withM, "merged", n)
            set(scores, withM, "type", gsub(paste0("-", n), "", scores[["type"]][withM]))
        }
    }

    if (!useGGPlot2)
    {
        oldp <- par(no.readonly = TRUE)
        maxStrW <- max(strwidth(unique(scores$type), units = 'in', cex = 0.9)) + 0.5
        omai <- par("mai")
        par(mai = c(maxStrW, 0.5, omai[3], 0))

        if (length(mcn) > 1)
            cols <- getBrewerPal(length(mcn), "Paired")
        else
            cols <- getBrewerPal(nrow(scores), "Paired")

        bpargs <- list(las = 2, col = cols, border = cols, cex.axis = 0.9, xpd = TRUE)

        if (length(mcn) > 1)
        {
            scSplit <- split(scores, by = "type", keep.by = FALSE)
            scSplit <- sapply(names(scSplit), function(mb) setnames(scSplit[[mb]], "score", mb), simplify = FALSE) # assign column names

            plotTab <- Reduce(function(left, right)
            {
                merge(left, right, by = "merged", all = TRUE)
            }, scSplit)

            plot.new()

            makeLegend <- function(x, y, ...) legend(x, y, plotTab$merged, col = cols, lwd = 1, xpd = NA, ncol = 1,
                                                     cex = 0.75, bty = "n", ...)

            # auto legend positioning: https://stackoverflow.com/a/34624632/9264518
            leg <- makeLegend(0, 0, plot = FALSE)
            # lh <- (grconvertY(leg$rect$h, to = "ndc") - grconvertY(0, to = "ndc")) * 1.75
            # par(omd = c(0, 1, 0, 1 - lh), new = TRUE)
            lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
            par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
            bpvals <- as.matrix(plotTab[, -"merged"])
            bp <- do.call(barplot, c(list(bpvals, beside = TRUE), bpargs))
            bpsc <- as.vector(bpvals)
            # legend("top", plotTab$merged, col = cols, lwd = 1, inset = c(0, -0.2), xpd = NA, horiz = TRUE, cex = 0.75)
            # makeLegend(0, par("usr")[4] + lh)
            makeLegend(par("usr")[2], par("usr")[4])
        }
        else
        {
            bp <- do.call(barplot, c(list(scores$score, names.arg = scores$type), bpargs))
            bpsc <- scores$score
        }

        text(bp, bpsc, labels = round(bpsc, 2), pos = 3, cex = 0.8, xpd = TRUE)

        par(oldp)
    }
    else
    {
        scorePlot <- ggplot(scores, aes_string(x = "type", y = "score")) +
            cowplot::theme_cowplot(font_size = 12) +
            theme(axis.title.y = element_blank(), axis.title.x = element_blank(), # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  legend.position = "top", legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
            guides(colour = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))


        if (length(mcn) > 1)
            scorePlot <- scorePlot + geom_bar(stat = "identity", position = "dodge",
                                              aes_string(colour = "merged", fill = "merged"))
        else
        {
            scorePlot <- scorePlot +
                geom_bar(stat = "identity", aes_string(colour = "type", fill = "type")) +
                theme(legend.position = "none")
        }

        # scorePlot <- scorePlot + coord_flip()

        return(scorePlot)
        # scorePlot <- ggplot(scores, aes_string(x = "type", y = "score")) +
        #     geom_bar(stat = "identity", position = "dodge", aes_string(colour = "merged", fill = "merged")) +
        #     theme_cowplot(font_size = 12) +
        #     theme(axis.title.y = element_blank(), axis.title.x = element_blank(), # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        #           legend.position = "top", legend.title = element_blank()) +
        #     # geom_text(aes(label = round(score, 2)), position = position_dodge(width = 0.5)) +
        #     #guides(colour = guide_legend(nrow = 3, ncol = 3, byrow = TRUE)) +
        #     coord_flip()
    }
})

#' @describeIn compounds Plots an annotated spectrum for a given candidate
#'   compound of a feature group.
#'
#' @param plotStruct If \code{TRUE} then the candidate structure is drawn in the
#'   spectrum.
#'
#' @export
setMethod("plotSpec", "compounds", function(obj, index, groupName, MSPeakLists, formulas = NULL,
                                            plotStruct = TRUE, useGGPlot2 = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertCount(index, positive = TRUE, add = ac)
    checkmate::assertString(groupName, min.chars = 1, add = ac)
    checkmate::assertClass(MSPeakLists, "MSPeakLists", add = ac)
    checkmate::assertClass(formulas, "formulas", null.ok = TRUE, add = ac)
    aapply(checkmate::assertFlag, . ~plotStruct + useGGPlot2, fixed = list(add = ac))
    checkmate::reportAssertions(ac)

    compTable <- compoundTable(obj)[[groupName]]

    if (is.null(compTable) || nrow(compTable) == 0)
        return(NULL)

    if (!is.null(formulas) && groupName %in% groupNames(formulas))
        fTable <- formulas[[groupName]][byMSMS == TRUE]
    else
        fTable <- NULL

    compr <- compTable[index, ]
    spec <- MSPeakLists[[groupName]][["MSMS"]]

    # merge formulas
    fi <- compr$fragInfo[[1]]
    if (!is.null(fTable) && nrow(fTable) > 0)
    {
        ft <- fTable[neutral_formula == compr$formula]
        if (is.null(fi))
        {
            fi <- ft
            fi[, mergedBy := list(list(algorithm(formulas)))]
        }
        else
            fi <- mergeFragInfo(fi, getFragmentInfoFromForms(spec, ft),
                                algorithm(obj), algorithm(formulas))
    }

    if (plotStruct)
        mol <- getMoleculesFromSMILES(compr$SMILES)

    if (!useGGPlot2)
    {
        oldp <- par(mar = par("mar") * c(1, 1, 0, 0))

        if (plotStruct && isValidMol(mol))
        {
            molHInch <- 1.5
            makeMSPlot(spec, fi, extraHeightInch = molHInch)
        }
        else
            makeMSPlot(spec, fi)

        # draw structure
        if (plotStruct && isValidMol(mol))
        {
            raster <- rcdk::view.image.2d(mol[[1]], rcdk::get.depictor(100, 100))
            if (!isEmptyMol(mol[[1]]))
            {
                img <- magick::image_trim(magick::image_read(raster))
                img <- magick::image_transparent(img, "white")
            }

            dpi <- (par("cra")/par("cin"))[1]

            startx <- par("usr")[1]
            xlim <- par("usr")[2]
            ylim <- par("usr")[4]
            imgInfo <- magick::image_info(img)

            imgPlotW <- xinch(imgInfo$width / dpi)
            imgPlotH <- yinch(imgInfo$height / dpi)

            maxW <- 0.2 * xlim
            if (imgPlotW > maxW)
            {
                hresize <- imgPlotW / maxW
                imgPlotH <- imgPlotH / hresize
                imgPlotW <- maxW
            }

            maxH <- yinch(molHInch)
            if (imgPlotH > maxH)
            {
                wresize <- imgPlotH / maxH
                imgPlotW <- imgPlotW / wresize
                imgPlotH <- maxH
            }

            # offset a little
            startx <- startx + 0.01 * xlim
            ylim <- ylim * 0.99

            rasterImage(img, startx, ylim - imgPlotH, startx + imgPlotW, ylim)
        }

        par(oldp)
    }
    else
    {
        MSPlot <- makeMSPlotGG(spec, fi)

        if (plotStruct && isValidMol(mol))
        {
            raster <- rcdk::view.image.2d(mol[[1]], rcdk::get.depictor(100, 100))
            img <- magick::image_trim(magick::image_read(raster))

            # positioning if legend is on top... doesn't work too well :(
            # mar <- 7
            # pos <- grid::convertUnit(grid::unit(0.8, "npc") - grid::unit(mar, "pt"), "npc", valueOnly = TRUE)
            # MSPlot <- MSPlot + theme(plot.margin = margin(mar, mar, mar, mar)) +
            #     ylim(0, max(fi$intensity) * (1.2 + (0.8 - pos)))

            pos <- 0.8; size <- 1 - pos
            MSPlot <- MSPlot +
                ylim(0, max(spec$intensity) * (1 + size + 0.1)) # add a bit of space between most intense point+label
            MSPlot <- cowplot::ggdraw(MSPlot) +
                cowplot::draw_image(img, pos, pos, size, size)
        }

        return(MSPlot)
        # mcn <- mergedCompoundNames(obj)
        # scoreCols <- getAllCompCols(getCompScoreColNames(), names(compTable), mcn)
        # scores <- setnames(transpose(compr[, scoreCols, with = FALSE]), "score")
        # scores[, type := scoreCols]
        # scores <- scores[!is.na(score)]
        #
        # if (length(mcn) > 1)
        # {
        #     scores[, merged := "both"]
        #     for (n in mcn)
        #     {
        #         withM <- which(grepl(paste0("-", n), scores[["type"]], fixed = TRUE))
        #         set(scores, withM, "merged", n)
        #         set(scores, withM, "type", gsub(paste0("-", n), "", scores[["type"]][withM]))
        #     }
        # }

        # scorePlot <- ggplot(scores, aes_string(x = "type", y = "score")) +
        #     theme_cowplot(font_size = 12) +
        #     theme(axis.title.y = element_blank(), axis.title.x = element_blank(), # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        #           legend.position = "top", legend.title = element_blank()) +
        #     guides(colour = guide_legend(nrow = 3, ncol = 2, byrow = TRUE))
        #
        #
        # if (length(mcn) > 1)
        #     scorePlot <- scorePlot + geom_bar(stat = "identity", position = "dodge",
        #                                       aes_string(colour = "merged", fill = "merged"))
        # else
        #     scorePlot <- scorePlot + geom_bar(stat = "identity", aes_string(colour = "type", fill = "type"))
        #
        # scorePlot <- scorePlot + coord_flip()

        # ap <- align_plots(MSPlot, scorePlot, align = "h")
        # plot_grid(ap[[1]], plot_grid(structPlot, ap[[2]], nrow = 2, rel_heights = c(1, 2)), rel_widths = c(2, 1))
        # gridExtra::grid.arrange(MSPlot, structPlot, scorePlot, layout_matrix = matrix(c(1, 2, 1, 3, 1, 3), ncol = 2, byrow = TRUE), widths = c(2, 1))

    }
})


setMethod("mergedCompoundNames", "compounds", function(compounds) character(0))
setMethod("mergedCompoundNames", "compoundsConsensus", function(compounds) compounds@mergedCompNames)

#' @describeIn compounds Generates a consensus from multiple \code{compounds}
#'   objects. The compound results will be merged and optionally compounds will
#'   be filtered by their occurrence throughout the specified \code{compounds}
#'   objects.
#'
#' @param \dots \code{compounds} objects that should be used to generate the
#'   consensus.
#' @param compThreshold Minimal fraction (0-1) that a candidate must be present
#'   within all given \code{compounds} objects.
#' @param mergeScoresFunc Function used to calculate the total score for all
#'   (merged) score columns.
#' @param minMaxNormalization Set to \code{TRUE} to apply min-max normalization
#'   of (merged) scoring columns. \code{FALSE} will apply normalization to the
#'   maximum value. Scorings with negative values will always be min-max
#'   normalized.
#'
#' @return \code{consensus} returns a \code{compounds} object that is produced
#'   by merging multiple specified \code{compounds} objects.
#'
#' @export
setMethod("consensus", "compounds", function(obj, ..., compThreshold = 0.0, minMaxNormalization = TRUE,
                                             mergeScoresFunc = sum)
{
    allCompounds <- c(list(obj), list(...))

    ac <- checkmate::makeAssertCollection()
    checkmate::assertList(allCompounds, types = "compounds", min.len = 2, any.missing = FALSE,
                          unique = TRUE, .var.name = "...", add = ac)
    checkmate::assertNumber(compThreshold, lower = 0, finite = TRUE, add = ac)
    checkmate::assertFunction(mergeScoresFunc, add = ac)
    checkmate::reportAssertions(ac)

    compNames <- sapply(allCompounds, algorithm)
    if (anyDuplicated(compNames))
    {
        # duplicate algorithms used, try to form unique names by adding library
        dbs <- lapply(allCompounds, function(cmp) # UNDONE: make this a method
        {
            for (res in compoundTable(cmp))
            {
                if (nrow(res) > 0)
                    return(res$database[[1]])
            }
        })

        compNames <- sapply(seq_along(allCompounds), function(cmpi) paste0(substr(algorithm(allCompounds[[cmpi]]), 1, 3), "-",
                                                                           substr(dbs[[cmpi]], 1, 3)))

        # in case names are still duplicated
        compNames <- make.unique(compNames)
    }


    # initialize all compound objects for merge: copy them, rename columns to
    # avoid duplicates and set merged by field of fragInfo.
    allCompTables <- lapply(seq_along(allCompounds), function(cmpi)
    {
        mergedBy <- list(list(compNames[[cmpi]])) # wrap in a double list as we want it to be a list of characters

        return(lapply(compoundTable(allCompounds[[cmpi]]), function(ct)
        {
            ret <- copy(ct)

            for (r in seq_len(nrow(ret)))
            {
                fi <- ret[["fragInfo"]][[r]]
                if (!is.null(fi) && nrow(fi) > 0 && is.null(fi[["mergedBy"]]))
                {
                    fi <- copy(fi)
                    set(fi, j = "mergedBy", value = mergedBy)
                    set(ret, r, "fragInfo", list(list(fi)))
                }
            }

            setnames(ret, paste0(names(ret), "-", compNames[[cmpi]]))

            return(ret)
        }))
    })

    # columns that should be unique (fragInfo and InChIKey1 are dealt separately)
    uniqueCols <- c("SMILES", "formula", "InChi", "InChIKey2", "InChIKey", "neutralMass")

    leftName <- compNames[[1]]
    mCompList <- allCompTables[[1]]
    for (compIndex in seq(2, length(allCompTables)))
    {
        rightName <- compNames[[compIndex]]

        printf("Merging %s with %s... ", paste0(compNames[seq_len(compIndex-1)], collapse = ","), rightName)

        rightTable <- allCompTables[[compIndex]]

        for (grp in union(names(mCompList), names(rightTable)))
        {
            if (is.null(rightTable[[grp]]))
                next # nothing to merge
            else if (is.null(mCompList[[grp]])) # not yet present
            {
                mCompounds <- rightTable[[grp]]

                # rename columns that should be unique from right to left
                unCols <- c(uniqueCols, "fragInfo", "InChIKey1")
                unCols <- unCols[sapply(unCols, function(uc) !is.null(mCompounds[[paste0(uc, "-", rightName)]]))]
                setnames(mCompounds, paste0(unCols, "-", rightName), paste0(unCols, "-", leftName))
            }
            else
            {
                # UNDONE: (merely) InChIKey1 is a sensible choice?
                mCompounds <- merge(mCompList[[grp]], rightTable[[grp]],
                                    by.x = paste0("InChIKey1-", leftName),
                                    by.y = paste0("InChIKey1-", rightName),
                                    all = TRUE)

                # merge fragment info
                fiColLeft <- paste0("fragInfo-", leftName)
                fiColRight <- paste0("fragInfo-", rightName)

                if (!is.null(mCompounds[[fiColLeft]]))
                {
                    for (r in seq_len(nrow(mCompounds)))
                    {
                        # use copy as a workaround for buggy nested data.tables
                        fiRLeft <- copy(mCompounds[[fiColLeft]][[r]])
                        fiRRight <- copy(mCompounds[[fiColRight]][[r]])
                        hasLeft <- length(fiRLeft) > 0 && nrow(fiRLeft) > 0
                        hasRight <- length(fiRRight) > 0 && nrow(fiRRight) > 0

                        if (hasLeft && hasRight)
                        {
                            # both have fraginfo
                            fiMerged <- mergeFragInfo(fiRLeft, fiRRight, leftName, rightName)
                            set(mCompounds, r, fiColLeft, list(list(fiMerged)))
                        }
                        else if (hasRight) # only right
                            set(mCompounds, r, fiColLeft, list(list(fiRRight)))
                    }

                    mCompounds[, (fiColRight) := NULL]
                }

                # remove duplicate columns that shouldn't
                for (col in uniqueCols)
                {
                    colLeft <- paste0(col, "-", leftName)
                    colRight <- paste0(col, "-", rightName)
                    if (!is.null(mCompounds[[colRight]]))
                    {
                        if (is.null(mCompounds[[colLeft]]))
                            setnames(mCompounds, colRight, colLeft)
                        else
                        {
                            mCompounds[, (colLeft) := ifelse(!is.na(get(colLeft)), get(colLeft), get(colRight))]
                            mCompounds[, (colRight) := NULL]
                        }
                    }
                }
            }

            mCompList[[grp]] <- mCompounds
        }

        cat("Done!\n")
    }

    printf("Determining coverage and final scores... ")

    # Determine coverage of compounds between objects and the merged score. The score column can be
    # used for the former as there is guaranteed to be one for each merged object.
    for (grpi in seq_along(mCompList))
    {
        # fix up de-duplicated column names
        deDupCols <- c(uniqueCols, c("fragInfo", "InChIKey1"))
        leftCols <- paste0(deDupCols, "-", leftName)
        deDupCols <- deDupCols[leftCols %in% names(mCompList[[grpi]])]
        leftCols <- leftCols[leftCols %in% names(mCompList[[grpi]])]
        if (length(leftCols) > 0)
            setnames(mCompList[[grpi]], leftCols, deDupCols)

        scoreCols <- grep("score-", colnames(mCompList[[grpi]]), value = TRUE)

        if (length(scoreCols) == 0) # nothing was merged
            mCompList[[grpi]][, coverage := 1 / length(allCompounds)]
        else
        {
            for (r in seq_len(nrow(mCompList[[grpi]])))
                set(mCompList[[grpi]], r, "coverage",
                    sum(sapply(scoreCols, function(c) !is.na(mCompList[[grpi]][[c]][r]))) / length(allCompounds))
        }

        if (compThreshold > 0)
            mCompList[[grpi]] <- mCompList[[grpi]][coverage >= compThreshold]

        # set scores after filtering: otherwise normalized values are out of date
        if (length(scoreCols) > 0 && nrow(mCompList[[grpi]]) > 0)
        {
            scores <- mCompList[[grpi]][, lapply(.SD, normalize, minMax = minMaxNormalization), .SDcols = scoreCols]

            for (r in seq_len(nrow(mCompList[[grpi]])))
            {
                scoreRow <- scores[r]
                scoreRow[is.na(scoreRow)] <- 0
                set(mCompList[[grpi]], r, "score", mergeScoresFunc(scoreRow))
            }
            setorder(mCompList[[grpi]], -score)
        }
    }

    cat("Done!\n")

    # prune empty/NULL results
    if (length(mCompList) > 0)
        mCompList <- mCompList[sapply(mCompList, function(r) !is.null(r) && nrow(r) > 0, USE.NAMES = FALSE)]

    return(compoundsConsensus(compounds = mCompList, algorithm = paste0(unique(sapply(allCompounds, algorithm)), collapse = ", "),
                              mergedCompNames = compNames))
})

#' @templateVar func generateCompounds
#' @templateVar what generate compounds
#' @templateVar ex1 generateCompoundsMetfrag
#' @templateVar ex2 generateCompoundsSirius
#' @templateVar algos metfrag,sirius
#' @template generic-algo
#'
#' @param ... Any parameters to be passed to the selected compound generation
#'   algorithm.
#'
#' @rdname compound-generation
#' @aliases generateCompounds
#' @export
setMethod("generateCompounds", "featureGroups", function(fGroups, MSPeakLists, algorithm, ...)
{
    f <- switch(algorithm,
                metfrag = generateCompoundsMetfrag,
                sirius = generateCompoundsSirius,
                stop("Invalid algorithm! Should be: metfrag or sirius"))

    f(fGroups, MSPeakLists, ...)
})
