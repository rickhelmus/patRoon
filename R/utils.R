#' @include main.R
NULL

printf <- function(...) cat(sprintf(...), sep = "")
fprintf <- function(file, ..., append = FALSE) cat(sprintf(...), sep = "", file = file, append = append)

if (!exists("hasName")) # should be defined in latest R versions
    hasName <- function(x, name) return(name %in% names(x))

# backslashes are converted to forwardslashes for *nix compat
baseName <- function(...) basename(gsub("\\\\", "/", ...))

simplifyAnalysisNames <- function(slist)
{
    # simplify sample file names: remove extension and path
    return(as.character(sapply(slist, function(s) baseName(tools::file_path_sans_ext(s)), USE.NAMES = F)))
}

getMzMLAnalysisPath <- function(file, path) file.path(path, paste0(file, ".mzML"))
getMzXMLAnalysisPath <- function(file, path) file.path(path, paste0(file, ".mzXML"))
getMzDataAnalysisPath <- function(file, path) file.path(path, paste0(file, ".mzData"))
getAnalysisPath <- function(file, path, ext) file.path(path, paste0(file, ".", ext))

getMzMLOrMzXMLAnalysisPath <- function(file, path)
{
    ret <- getMzMLAnalysisPath(file, path)
    if (!file.exists(ret))
        ret <- getMzXMLAnalysisPath(file, path)
    return(ret)
}

makeFGroupName <- function(id, ret, mz) sprintf("M%.0f_R%.0f_%d", mz, ret, id)

mkdirp <- function(path)
{
    for (p in path)
    {
        if (!dir.exists(p))
            dir.create(p, recursive = TRUE)
    }
}

insertDFColumn <- function(df, col, d, before)
{
    if (before < ncol(df))
    {
        ind <- c()

        if (before > 1)
            ind <- append(ind, 1:(before-1))

        ind <- append(ind, c(ncol(df), before:(ncol(df)-1)))
        df <- df[, ind]
    }

    return(df)
}

insertDTColumn <- function(dt, col, d, before)
{
    dt <- copy(dt)

    dt[, (col) := d]

    if (before < ncol(dt))
    {
        ind <- c()

        if (before > 1)
            ind <- append(ind, 1:(before-1))

        ind <- append(ind, c(ncol(dt), before:(ncol(dt)-1)))
        setcolorder(dt, ind)
    }

    return(dt)
}

# Fix from R DescTools
#' Internal fix for \pkg{RDCOMClient}, ignore.
#' @param ref,className ignore
#' @keywords internal
#' @export
createCOMReference <- function(ref, className)
{
    RDCOMClient::createCOMReference(ref, className)
}

#' Internal fix for \pkg{RDCOMClient}, ignore.
#' @param ... ignore
#' @keywords internal
#' @export
COMStop <- function(...)
{
    RDCOMClient::COMStop(...)
}

checkPackage <- function(pkg, gh = NULL)
{
    # from http://stackoverflow.com/a/20333756
    if (!requireNamespace(pkg, quietly = TRUE))
    {
        if (!is.null(gh))
            stop(sprintf("Please install %s from github: remotes::install_github('%s')", pkg, gh))
        else
            stop(sprintf("Please install %s: install.packages('%s')", pkg, pkg))
    }
}

executeCommand <- function(cmd, args = character(), ...)
{
    return(system2(cmd, sapply(args, shQuote), ...))
}

# NOTE: keep in sync with install-patRoon version
getCommandWithOptPath <- function(cmd, opt, verify = TRUE)
{
    if (Sys.info()[["sysname"]] == "Windows")
        cmd <- paste0(cmd, ".exe") # add file extension for Windows

    opt <- paste0("patRoon.path.", opt)
    path <- getOption(opt)
    if (!is.null(path) && nzchar(path))
    {
        ret <- file.path(path.expand(path), cmd)
        if (!file.exists(ret))
        {
            if (verify)
                stop(sprintf("Cannot find '%s'. Is the option '%s' set correctly?", ret, opt))
            return(NULL)
        }

        return(ret)
    }

    # assume command is in PATH --> no need to add path
    if (!nzchar(Sys.which(cmd)))
    {
        if (verify)
            stop(sprintf("Cannot find '%s'. Either add the correct file location to the PATH environment variable or set '%s' with options().", cmd, opt))
        return(NULL)
    }

    return(cmd)
}

# NOTE: keep in sync with install-patRoon version
findPWizPath <- function()
{
    # try to find ProteoWizard
    # order: options --> win registry --> PATH
    # the PATH is searched last because OpenMS might have added its own old version.

    path <- getOption("patRoon.path.pwiz")
    if (!is.null(path) && nzchar(path))
        return(path)

    if (Sys.info()[["sysname"]] == "Windows")
    {
        # Inspired by scan_registry_for_rtools() from pkgload
        key <- "Directory\\shell\\Open with SeeMS\\command"
        reg <- tryCatch(utils::readRegistry(key, "HCR"), error = function(e) NULL)

        # not sure if this might occur
        if (is.null(reg))
            reg <- tryCatch(utils::readRegistry(key, "HLM"), error = function(e) NULL)

        if (!is.null(reg))
        {
            path <- tryCatch(dirname(sub("\"([^\"]*)\".*", "\\1", reg[[1]])), error = function(e) NULL)
            if (!is.null(path) && file.exists(file.path(path, "msconvert.exe"))) # extra check: see if msconvert is there
                return(path)
        }
    }

    # check PATH
    msc <- if (Sys.info()[["sysname"]] == "Windows") "msconvert.exe" else "msconvert"
    path <- dirname(Sys.which(msc))
    if (nzchar(path))
        return(path)

    return(NULL)
}

# From http://stackoverflow.com/a/30835971
line2user <- function(line, side)
{
    lh <- par('cin')[2] * par('cex') * par('lheight')
    x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
           `1` = grconvertY(-line * y_off, 'npc', 'user'),
           `2` = grconvertX(-line * x_off, 'npc', 'user'),
           `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
           `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
           stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

#' @details \code{generateAnalysisInfo} is an utility function that
#'   automatically generates an analysis information object. It will collect all
#'   datafiles from given file paths and convert the filenames into valid
#'   analysis names (\emph{i.e.} without extensions such as \file{.d} and
#'   \file{.mzML}). Duplicate analyses, which may appear when datafiles with
#'   different file extension (\file{.d}, \file{.mzXML} and/or \file{.mzML}) are
#'   present, will be automatically removed.
#'
#' @param paths A character vector containing one or more file paths that should
#'   be used for finding the analyses.
#' @param groups,blanks An (optional) character vector containing replicate
#'   groups and references, respectively (will be recycled). If \code{groups} is
#'   an empty character string (\code{""}) the analysis name will be set as
#'   replicate group.
#' @param concs An optional numeric vector containing concentration values for
#'   each analysis. Can be \code{NA} if unknown. If the length of \code{concs}
#'   is less than the number of analyses the remainders will be set to
#'   \code{NA}. Set to \code{NULL} to not include concentration data.
#' @param formats A character vector of analyses file types. Valid values are:
#'   \code{Bruker}, \code{mzXML} and \code{mzML}.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfo <- function(paths, groups = "", blanks = "", concs = NULL,
                                 formats = MSFileFormats())
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDirectoryExists(paths, access = "r", add = ac)
    checkmate::assertCharacter(groups, min.len = 1, add = ac)
    checkmate::assertCharacter(blanks, min.len = 1, add = ac)
    checkmate::assertNumeric(concs, finite = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertSubset(formats, MSFileFormats(), empty.ok = FALSE, add = ac)
    checkmate::reportAssertions(ac)

    files <- listMSFiles(paths, formats)

    if (length(files) == 0)
    {
        warning(sprintf("No analyses found in %s!", paste(paths, collapse = ", ")))
        return(NULL)
    }

    ret <- data.frame(path = dirname(files), analysis = simplifyAnalysisNames(files), stringsAsFactors = FALSE)
    ret <- ret[!duplicated(ret[, c("path", "analysis")]), ]

    # set after duplicate removal
    groups <- rep(groups, length.out = nrow(ret))
    ret$group <- ifelse(!nzchar(groups), ret$analysis, groups)
    ret$blank <- blanks

    if (!is.null(concs))
    {
        if (length(concs) >= nrow(ret))
            ret$conc <- concs[seq_len(nrow(ret))]
        else
        {
            ret$conc <- NA
            ret$conc[seq_along(concs)] <- concs
        }
    }

    return(ret)
}

#' @details \code{generateAnalysisInfoFromEnviMass} loads analysis information
#'   from an \pkg{enviMass} project. Note: this funtionality has only been
#'   tested with older versions of \pkg{enviMass}.
#'
#' @param path The path of the enviMass project.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfoFromEnviMass <- function(path)
{
    checkmate::assertDirectoryExists(path, access = "r")

    enviSInfo <- fread(file.path(path, "dataframes", "measurements"))[Type %in% c("sample", "blank")]

    enviSInfo[, DT := as.POSIXct(paste(Date, Time, sep = " "))]
    enviSInfo[, group := tag3]

    blanks <- enviSInfo[enviSInfo$Type == "blank", ]
    if (nrow(blanks) > 0)
    {
        blanks[, blank := paste0("blank", match(DT, unique(DT)))] # add unique date+time identifier
        enviSInfo[, blank := sapply(DT, function(dt)
        {
            bls <- blanks[DT <= dt]
            if (nrow(bls) > 0)
                return(bls[which.max(bls$DT), blank])
            else
                return("")
        })]

        enviSInfo[Type == "blank", group := blank]
    }
    else
        enviSInfo[, blank := ""]

    enviSInfo[group == "FALSE", group := ""]

    ret <- data.frame(path = file.path(path, "files"), analysis = enviSInfo$ID, group = enviSInfo$group,
                      blank = enviSInfo$blank, stringsAsFactors = FALSE)

    return(ret)
}

# based on http://www.cureffi.org/2013/09/23/a-quick-intro-to-chemical-informatics-in-r/
getRCDKStructurePlot <- function(molecule, width = 500, height = 500, trim = TRUE, transparent = TRUE)
{
    # drawing unconnected is not supported. Assume these are salts and we can just draw largest instead
    # UNDONE: is this still relevant?
    # if (!is.connected(molecule))
    #     molecule <- get.largest.component(molecule)

    img <- rcdk::view.image.2d(molecule, rcdk::get.depictor(width, height)) # get Java representation into an image matrix.
    img <- magick::image_read(img)

    if (!isEmptyMol(molecule))
    {
        if (trim)
            img <- magick::image_trim(img)
        if (transparent)
            img <- magick::image_transparent(img, "white")
    }

    return(img)
}


unFactorDF <- function(df)
{
    # from http://stackoverflow.com/a/2853231
    i <- sapply(df, is.factor)
    df[i] <- lapply(df[i], as.character)
    return(df)
}

# from http://stackoverflow.com/a/25455968
loadRData <- function(fileName)
{
    load(fileName)
    get(ls()[ls() != "fileName"])
}

# from https://stackoverflow.com/a/34788691
jgc <- function()
{
    gc()
    rJava::.jcall("java/lang/System", method = "gc")
}

# based on https://stackoverflow.com/a/17972280
recursiveApplyDT <- function(l, f, appl = lapply, ...)
{
    rec <- function(x)
    {
        if (inherits(x, "data.table"))
            f(x)
        else if (is.recursive(x))
            appl(x, rec, ...)
        else
            x
    }

    return(rec(l))
}

prepareDTForComparison <- function(dt)
{
    setattr(dt, ".internal.selfref", NULL)
    setindex(dt, NULL)
}

getCSVStr <- function(df, ...)
{
    tcon <- textConnection("ret", "w", TRUE)
    write.csv(df, tcon, ...)
    close(tcon)
    return(paste0(ret, collapse = "\n"))
}

readAllFile <- function(f) readChar(f, file.size(f))

# NOTE: args converted to list, so sep argument to paste doesn't make sense
pasteNonEmpty <- function(..., collapse = NULL)
{
    strl <- list(...)
    strl <- as.character(strl[!sapply(strl, is.null)])
    paste0(strl[nzchar(strl) > 0], collapse = collapse)
}

getBrewerPal <- function(n, name)
{
    maxn <- RColorBrewer::brewer.pal.info[name, "maxcolors"]
    if (n > maxn)
        return(colorRampPalette(RColorBrewer::brewer.pal(maxn, name))(n))

    if (n < 3)
        return(RColorBrewer::brewer.pal(3, name)[seq_len(n)])

    return(RColorBrewer::brewer.pal(n, name))
}

getArgNames <- function(..., def = NULL)
{
    args <- sapply(substitute(list(...))[-1], deparse)
    n <- names(args)
    if (is.null(def))
        def <- args
    if (!is.null(n))
        args <- ifelse(nzchar(n), n, def)
    else if (!is.null(def))
        args <- def
    return(args)
}

getStrListWithMax <- function(l, m, collapse)
{
    if (length(l) > m)
    {
        l <- l[seq_len(m + 1)]
        l[m + 1] <- "..."
    }
    paste0(l, collapse = collapse)
}

showAnaInfo <- function(anaInfo)
{
    rGroups <- unique(anaInfo$group)
    blGroups <- unique(anaInfo$blank)
    printf("Analyses: %s (%d total)\n", getStrListWithMax(anaInfo$analysis, 6, ", "), nrow(anaInfo))
    printf("Replicate groups: %s (%d total)\n", getStrListWithMax(rGroups, 8, ", "), length(rGroups))
    printf("Replicate groups used as blank: %s (%d total)\n", getStrListWithMax(blGroups, 8, ", "), length(blGroups))
}

showObjectSize <- function(object) printf("Object size (indication): %s\n", format(utils::object.size(object), "auto", "SI"))

allSame <- function(l)
{
    if (length(l) > 1)
    {
        if (all(is.na(l)))
            return(TRUE)
        if (any(is.na(l)))
            return(FALSE)

        return(all(l[[1]] == l))
    }

    return(TRUE)
}

normalize <- function(x, minMax, xrange = range(x, na.rm = TRUE))
{
    xn <- x[!is.na(x)]
    if (length(xn) == 0 || all(xn == 0))
        return(x) # all NA or all values zero

    if (allSame(xn))
        xn <- rep(as(1, typeof(xn[[1]])), length(xn))
    else
    {
        minv <- xrange[1]
        if (!minMax)
            minv <- min(minv, 0) # force minMax if min <0

        xn <- (xn - minv) / (xrange[2] - minv)
    }

    x[!is.na(x)] <- xn

    return(x)
}

# from https://stackoverflow.com/a/38228840
countCharInStr <- function(str, ch) sum(charToRaw(str) == charToRaw(ch))

makeVennPlot <- function(plotObjects, categories, areas, intersectFunc,
                         intersectLenFunc, ...)
{
    nobj <- length(plotObjects)
    areas <- unname(areas) # do.call below won't work with names

    if (all(areas == 0))
        stop("Cannot plot Venn when all objects are empty")

    fill <- getBrewerPal(nobj, "Paired")
    vennArgs <- list(category = categories, lty = rep("blank", nobj), alpha = rep(0.5, nobj),
                     cex = 1.5, cat.cex = 1.5, fill = fill)
    vennArgs <- modifyList(vennArgs, list(...))

    getIntersectCounts <- function(inters) sapply(inters,
                                                  function(i) intersectLenFunc(Reduce(intersectFunc, plotObjects[i])))

    grid::grid.newpage() # need to clear plot region manually
    # plot.new()

    if (nobj == 1)
        gRet <- do.call(draw.single.venn, c(list(area = areas), vennArgs))
    else if (nobj == 2)
    {
        icounts <- getIntersectCounts(list(c(1, 2)))
        gRet <- do.call(VennDiagram::draw.pairwise.venn,
                        c(areas, icounts, list(rotation.degree = if (areas[1] < areas[2]) 180 else 0), vennArgs))
    }
    else if (nobj == 3)
    {
        icounts <- getIntersectCounts(list(c(1, 2), c(2, 3), c(1, 3), c(1, 2, 3)))
        gRet <- do.call(VennDiagram::draw.triple.venn, c(areas, icounts, vennArgs))
    }
    else if (nobj == 4)
    {
        icounts <- getIntersectCounts(list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4),
                                           c(1, 2, 3), c(1, 2, 4), c(1, 3, 4), c(2, 3, 4), c(1, 2, 3, 4)))
        gRet <- do.call(VennDiagram::draw.quad.venn, c(areas, icounts, vennArgs))
    }
    else if (nobj == 5)
    {
        icounts <- getIntersectCounts(list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(2, 3), c(2, 4), c(2, 5), c(3, 4),
                                           c(3, 5), c(4, 5), c(1, 2, 3), c(1, 2, 4), c(1, 2, 5), c(1, 3, 4), c(1, 3, 5),
                                           c(1, 4, 5), c(2, 3, 4), c(2, 3, 5), c(2, 4, 5), c(3, 4, 5),
                                           c(1, 2, 3, 4), c(1, 2, 3, 5), c(1, 2, 4, 5), c(1, 3, 4, 5), c(2, 3, 4, 5),
                                           c(1, 2, 3, 4, 5)))
        gRet <- do.call(VennDiagram::draw.quintuple.venn, c(areas, icounts, vennArgs))
    }

    invisible(list(gList = gRet, areas = areas, intersectionCounts = icounts))
}

getMSPlotData <- function(spec)
{
    hasFragInfo <- !is.null(spec[["formula"]])
    plotData <- copy(spec)

    # default colour/line width
    plotData[, c("colour", "lwd", "legend") := .("grey", 1, "unassigned")]

    if (hasFragInfo)
    {
        isAnnotated <- plotData[!is.na(formula), which = TRUE]
        if (!is.null(spec[["mergedBy"]]))
        {
            plotData[isAnnotated, legend := sapply(mergedBy, wrapStr, width = 10)]

            mbsUnique <- unique(plotData$legend)
            # order from small to big based on number of commas
            mbsUnique <- mbsUnique[order(sapply(mbsUnique, countCharInStr, ch = ",", USE.NAMES = FALSE))]
            mbCombCols <- setNames(getBrewerPal(length(mbsUnique), "Paired"), mbsUnique)

            plotData[isAnnotated, c("colour", "lwd") := .(mbCombCols[match(legend, mbsUnique)], 2)]
        }
        else
            plotData[isAnnotated, c("colour", "lwd", "legend") := .("blue", 2, "assigned")] # nothing merged, just mark all annotated blue
    }

    # mark precursor
    plotData[precursor == TRUE, c("colour", "lwd", "legend") := .("red", 2, "precursor")]

    return(plotData)
}

makeScoresPlot <- function(scoreTable, mcn, useGGPlot2)
{
    scores <- setnames(transpose(scoreTable), "score")
    scores[, type := names(scoreTable)]
    scores <- scores[!is.na(score)]

    if (length(mcn) > 1)
    {
        scores[, merged := "consensus"]
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
            cols <- getBrewerPal(length(unique(scores$merged)), "Paired")
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

            makeLegend <- function(x, y, ...) legend(x, y, unique(scores$merged), col = cols, lwd = 1, xpd = NA, ncol = 1,
                                                     cex = 0.75, bty = "n", ...)

            # auto legend positioning: https://stackoverflow.com/a/34624632/9264518
            leg <- makeLegend(0, 0, plot = FALSE)
            lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
            par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
            bpvals <- as.matrix(plotTab[, -"merged"])
            bp <- do.call(barplot, c(list(bpvals, beside = TRUE), bpargs))
            bpsc <- as.vector(bpvals)
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

        return(scorePlot)
    }
}

# spec may be annotated
makeMSPlot <- function(spec, xlim, ylim, ..., extraHeightInch = 0)
{
    plotData <- getMSPlotData(spec)
    if (!is.null(plotData[["formula"]]))
        plotData[!is.na(formula), formWidth := strwidth(formula, units = "inches")]

    if (is.null(xlim))
        xlim <- range(plotData$mz) * c(0.9, 1.1)
    else
        plotData <- plotData[numGTE(mz, xlim[1]) & numLTE(mz, xlim[2])] # remove any peaks outisde plotting range

    doLegend <- !is.null(plotData[["legend"]]) && any(!is.na(plotData[["legend"]]) & nzchar(plotData[["legend"]]))
    if (doLegend)
    {
        makeLegend <- function(x, y, ...)
        {
            legTab <- unique(plotData, by = "legend")[!is.na(legend)]
            # order from small to big based on number of commas
            legTab <- legTab[order(sapply(legTab$legend, countCharInStr, ch = ",", USE.NAMES = FALSE))]
            legend(x, y, legTab$legend, col = legTab$colour, xpd = NA, bty = "n",
                   text.col = legTab$colour, lty = 1, cex = 0.75, ...)
        }

        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc")
        oldp <- par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }

    if (is.null(ylim))
    {
        formWidths <- if (!is.null(plotData[["formWidth"]])) plotData$formWidth else rep(0, nrow(plotData))
        formWidths[is.na(formWidths)] <- 0

        # see how much extra vertical space is needed by formula labels
        # get character widths (assuming that height of vertically plotted text is the same)
        pheight <- par("din")[2] # plot dev height in inches
        relFormHeights <- formWidths / pheight

        ym <- max(plotData$intensity) # 'regular' plot height in user coordinates, not considering any labels etc
        relIntHeights <- plotData$intensity / ym
        maxRelH <- max(relFormHeights + relIntHeights)

        ym <- ym * maxRelH * 1.05 # enlarge y limit and add some extra spacing

        if (extraHeightInch > 0)
            ym <- ym * (1 + (extraHeightInch / pheight))

        ylim <- c(0, ym)
    }

    plot(0, xlab = "m/z", ylab = "Intensity", xlim = xlim, ylim = ylim,
         type = "n", bty = "l", ...)

    segments(plotData[["mz"]], 0, plotData[["mz"]], plotData[["intensity"]],
             col = plotData[["colour"]], lwd = plotData[["lwd"]])
        
    if (!is.null(plotData[["formula"]]))
    {
        pdf <- plotData[!is.na(formula)]
        if (nrow(pdf) > 0)
            text(pdf[["mz"]], pdf[["intensity"]] + (ylim[2] * 0.02),
                 subscriptFormula(pdf[["formula"]]), srt = 90, adj = 0)
    }

    if (doLegend)
    {
        makeLegend(par("usr")[2], par("usr")[4])
        par(oldp)
    }
}

makeMSPlotGG <- function(spec, ...)
{
    plotData <- getMSPlotData(spec)

    # BUG: throws errors when parse=TRUE and all labels are empty
    if (!is.null(plotData$formula) && any(!is.na(plotData$formula)))
    {
        # convert NAs to empty strings to avoid warnings with ggrepel
        plotData[!is.na(formula), formula := subscriptFormula(formula, parse = FALSE)]
        # BUG: ggrepel warns about NAs, but when replacing this with empty
        # characters or space the text entries are simply removed. for now just
        # use a dummy expression...
        plotData[is.na(formula), formula := "plain()"]
        ret <- ggplot(plotData, aes_string(x = "mz", y = 0, label = "formula")) +
            ggrepel::geom_text_repel(aes_string(y = "intensity", angle = 0), min.segment.length = 0.1, parse = TRUE,
                                     nudge_y = grid::convertUnit(grid::unit(5, "mm"), "npc", valueOnly = TRUE), size = 3.2)
    }
    else
        ret <- ggplot(plotData, aes_string(x = "mz", y = 0))

    ret <- ret + xlim(range(spec$mz) * c(0.9, 1.1)) +
        geom_segment(aes_string(xend = "mz", yend = "intensity", colour = "legend", size = "lwd")) +
        scale_size(range = c(0.5, 2), guide = FALSE) + xlab("m/z") + ylab("Intensity") +
        cowplot::theme_cowplot(font_size = 12) + theme(legend.position = "bottom", legend.title = element_blank())

    return(ret)
}

curTimeMS <- function() as.numeric(Sys.time()) * 1000

doDynamicTreeCut <- function(dendro, maxTreeHeight, deepSplit,
                             minModuleSize)
{
    if (minModuleSize == 1)
    {
        # workaround adapted from RAMClustR (ramclustR.R)
        ret <- dynamicTreeCut::cutreeDynamicTree(dendro = dendro, maxTreeHeight = maxTreeHeight,
                                                 deepSplit = deepSplit, minModuleSize = 2)
        single <- which(ret == 0) # all unassigned have length = 1
        ret[single] <- max(ret) + seq_len(length(single))
    }
    else
        ret <- dynamicTreeCut::cutreeDynamicTree(dendro = dendro, maxTreeHeight = maxTreeHeight,
                                                 deepSplit = deepSplit, minModuleSize = minModuleSize)

    return(ret)
}

plotDendroWithClusters <- function(dendro, ct, pal, colourBranches, showLegend, ...)
{
    if (colourBranches)
    {
        ct <- ct[order.dendrogram(dendro)] # re-order for dendrogram
        nclust <- length(unique(ct[ct != 0])) # without unassigned (in case of dynamicTreeCut)
        cols <- getBrewerPal(nclust, pal)
        dendro <- dendextend::color_branches(dendro, clusters = ct, col = cols[unique(ct)]) # unique: fixup colour order
        lcols <- dendextend::get_leaves_branches_col(dendro)
        dendextend::labels_colors(dendro) <- lcols
    }

    if (colourBranches && showLegend)
    {
        mr <- par("mar")
        mr[4] <- max(5.5, mr[4])
        withr::with_par(list(mar = mr),
        {
            plot(dendro, ...)
            legend("topright", legend = seq_len(nclust),
                   bty = "n", cex = 1, fill = cols, inset = c(-0.22, 0), xpd = NA,
                   ncol = 2, title = "cluster")
        })
    }
    else
        plot(dendro, ...)
}

doPlotSilhouettes <- function(clust, distm, kSeq, pch = 16, type = "b", ...)
{
    if (length(distm) == 1)
        stop("Need >2 clustered compounds")

    meanws <- sapply(kSeq, function(k)
    {
        sil <- silhouette(cutree(clust, k = k), distm)
        summary(sil)$avg.width
    })

    plot(x = kSeq, y = meanws, pch = pch, type = type, xaxt = "none",
         xlab = "cluster k", ylab = "mean silhouette width", ...)
    axis(1, kSeq)
    abline(v = kSeq[which.max(meanws)], lty = 2)
}


isValidMol <- function(mol) !is.null(mol) # && !is.na(mol)
emptyMol <- function() rcdk::parse.smiles("")[[1]]
isEmptyMol <- function(mol) rcdk::get.atom.count(mol) == 0

getMoleculesFromSMILES <- function(SMILES, doTyping = FALSE, doIsotopes = FALSE, emptyIfFails = FALSE)
{
    # vectorization doesn't work if any of the SMILES are not OK
    # mols <- rcdk::parse.smiles(SMILES)
    mols <- lapply(SMILES, function(sm)
    {
        ret <- rcdk::parse.smiles(sm)[[1]]
        if (!isValidMol(ret))
            ret <- rcdk::parse.smiles(sm, kekulise = FALSE)[[1]] # might work w/out kekulization
        if (!isValidMol(ret))
        {
            warning(paste("Failed to parse SMILES:", sm))
            if (emptyIfFails)
                ret <- emptyMol()
        }
        else
        {
            if (doTyping)
            {
                rcdk::do.typing(ret)
                rcdk::do.aromaticity(ret)
            }
            if (doIsotopes)
                rcdk::do.isotopes(ret)
        }
        return(ret)
    })
    
    return(mols)
}

getNeutralMassFromSMILES <- function(SMILES, mustWork = TRUE)
{
    sapply(getMoleculesFromSMILES(SMILES, doTyping = TRUE, doIsotopes = TRUE), function(m)
    {
        if (!mustWork)
        {
            if(!isValidMol(m))
                return(NA)
            
            # this could fail for some edge cases
            return(tryCatch(rcdk::get.exact.mass(m), error = function(e) NA))
        }
        return(rcdk::get.exact.mass(m))
    })
}

getMostIntenseAnaWithMSMS <- function(fGroups, MSPeakLists, groupName)
{
    gTable <- groups(fGroups)
    anaInfo <- analysisInfo(fGroups)
    pLists <- peakLists(MSPeakLists)

    ogind <- order(-gTable[[groupName]])
    oanalyses <- anaInfo$analysis[ogind] # analyses ordered from highest to lowest intensity

    # filter out analyses without MS/MS
    oanalyses <- sapply(oanalyses, function(a) if (!is.null(pLists[[a]][[groupName]][["MSMS"]])) a else "", USE.NAMES = FALSE)
    oanalyses <- oanalyses[oanalyses != ""]

    if (length(oanalyses) < 1)
        return(NULL)

    return(oanalyses[[1]])
}

# don't use all.equal here so functions are vectorized
# NOTE: actually all.equal seems to also take relative
# tolerances into account and thus may give different results.
numEQ <- function(x, y, tol = sqrt(.Machine$double.eps)) abs(x - y) <= tol
numGTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x > y
numLTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x < y

wrapStr <- function(s, width, sep = "\n") paste0(strwrap(s, width), collapse = sep)

pruneList <- function(l, checkEmptyElements = FALSE, checkZeroRows = FALSE)
{
    ret <- l[!sapply(l, is.null)]
    if (checkEmptyElements)
        ret <- ret[lengths(ret) > 0]
    if (checkZeroRows)
        ret <- ret[sapply(ret, nrow) > 0]
    return(ret)
}

# based on tabular() from formatting vignette of roxygen
tabularRD <- function(df, ...)
{
    align <- function(x) if (is.numeric(x)) "r" else "l"
    col_align <- vapply(df, align, character(1))

    # add headers
    df <- rbind(sprintf("\\strong{%s}", names(df)), df)

    cols <- lapply(df, format, ...)

    contents <- do.call("paste",
                        c(cols, list(sep = " \\tab ", collapse = "\\cr\n  ")))

    paste("\\tabular{", paste(col_align, collapse = ""), "}{\n  ",
          contents, "\n}\n", sep = "")
}

# make a S4 class inheritance tree in a format compatible with data.tree::FromListSimple()
makeClassHierarchy <- function(class, showParents)
{
    cldef <- getClassDef(class)
    subcl <- list()
    if (length(cldef@subclasses) > 0)
        subcl <- cldef@subclasses[sapply(cldef@subclasses, function(sc) sc@distance == 1)]

    ret <- c(list(name = class),
             lapply(names(subcl), makeClassHierarchy, showParents = FALSE))

    if (showParents)
    {
        pars <- selectSuperClasses(class)
        if (length(pars) > 0)
            ret <- c(list(name = paste0(pars, collapse = ", ")), list(ret))
    }

    return(ret)
}

printClassHierarchy <- function(class, showParents = TRUE, RD = FALSE)
{
    doPrintTxt <- function(cl, level)
    {
        indent <- strrep(" ", (level + 1) * 2)

        if (level > 0)
            cat(paste0(indent, "|-- ", cl$name))
        else
            cat(cl$name)
        cat("\n")

        for (clsub in cl)
        {
            if (is.list(clsub))
                doPrintTxt(clsub, level + 1)
        }
    }

    doPrintRD <- function(cl, level)
    {
        indent <- strrep(" ", (level + 1) * 2)

        printf("%s%s\\item{%s}\n", if (level == 0) "\\itemize{\n" else "", indent, cl$name)

        more <- any(sapply(cl, is.list))
        if (more)
            cat(paste0(indent, "\\itemize{\n"))

        for (clsub in cl)
        {
            if (is.list(clsub))
                doPrintRD(clsub, level + 1)
        }

        if (more)
            cat(paste0(indent, "}\n"))
        if (level == 0)
            cat("}\n")
    }

    hier <- makeClassHierarchy(class, showParents)
    hasParents <- hier$name != class

    if (RD)
    {
        hier <- rapply(hier, function(h) sprintf("\\code{\\link{%s}}", h),
                       classes = "character", how = "replace")

        # if parents are shown make the second line (ie this class) bold
        if (hasParents)
            hier[[2]]$name <- sprintf("\\strong{%s}", hier[[2]]$name)
        else
            hier$name <- sprintf("\\strong{%s}", hier$name)

        doPrintRD(hier, 0)
    }
    else
        doPrintTxt(hier, 0)

    invisible(NULL)
}

getAllMethods <- function(gen)
{
    # automatically retrieve defined methods for a generic and create document
    # links. This only works if the arguments of the method are named obj, objX or x.

    cl <- showMethods(gen, where = "package:patRoon", inherited = FALSE, printTo = FALSE,
                      classes = getClasses(asNamespace("patRoon")))
    cl <- cl[grepl("obj.*|x=", cl)]
    cl <- gsub("[^\"]*\"([^\"]*)\"[^\"]*", "\\1,", cl)
    # cl <- cl[!grepl("ANY", cl)]
    cl <- gsub(" ", "", cl)
    cl <- gsub(",$", "", cl)

    return(cl[order(tolower(cl))])
}

NULLToZero <- function(x) if (is.null(x)) 0 else x
zeroToNULL <- function(x) if (is.numeric(x) && x == 0) NULL else x

#' Verifies if all dependencies are installed properly and instructs the user if
#' this is not the case.
#' @export
verifyDependencies <- function()
{
    # UNDONE: for now just check one command-line tool of a software package
    # UNDONE: skip GenForm for now? Should be present as embedded binary.

    OK <- TRUE
    check <- function(name, path, opt, isDir = FALSE)
    {
        pleaseSet <- sprintf("Please set the '%s' option.", opt)
        printf("Checking %s... ", name)
        if (is.null(path) || !nzchar(path))
        {
            cat("not found or configured.", pleaseSet, "\n")
            OK <<- FALSE
        }
        else if (isDir)
        {
            if (!dir.exists(path))
            {
                cat("configured directory path does not exist!", pleaseSet, "\n")
                OK <<- FALSE
            }
            else
                printf("found! (directory '%s')\n", path)
        }
        else
        {
            # NOTE: dirname point to current path if getCommandWithOptPath() found it in PATH
            dn <- dirname(path)
            # if ((nzchar(dn) && dn != "." && !file.exists(path)) || !nzchar(Sys.which(path)))
            if (nzchar(dn) && dn != "." && !file.exists(path))
            {
                cat("configured path does not exist!", pleaseSet, "\n")
                OK <<- FALSE
            }
            else if (nzchar(dn) && dn != ".")
                printf("found! (in '%s')\n", dn)
            else
                cat("found!\n")
        }
    }

    check("ProteoWizard", findPWizPath(), "patRoon.path.pwiz", isDir = TRUE)
    check("OpenMS", getCommandWithOptPath("FeatureFinderMetabo", "OpenMS", verify = FALSE), "patRoon.path.OpenMS")
    check("pngquant", getCommandWithOptPath("pngquant", "pngquant", verify = FALSE), "patRoon.path.pngquant")
    check("SIRIUS", getCommandWithOptPath(getSiriusBin(), "SIRIUS", verify = FALSE), "patRoon.path.SIRIUS")
    check("MetFrag CL", getOption("patRoon.path.MetFragCL"), "patRoon.path.MetFragCL")
    check("MetFrag CompTox Database", getOption("patRoon.path.MetFragCompTox"), "patRoon.path.MetFragCompTox")
    check("MetFrag PubChemLite Database", getOption("patRoon.path.MetFragPubChemLite"), "patRoon.path.MetFragPubChemLite")
    check("OpenBabel", getCommandWithOptPath("obabel", "obabel", verify = FALSE), "patRoon.path.obabel")

    if (!OK)
        cat("\nSome dependencies were not found. Please make sure that their file locations are configured properly.",
            "For instance, run the following to set the location of MetFragCL:",
            sprintf("options(patRoon.path.MetFragCL = \"C:/MetFrag2.4.5-CL.jar\")"),
            "\nPlease see ?patRoon for more information on how to configure patRoon options.",
            sep = "\n")

    invisible(NULL)
}

# From https://stackoverflow.com/a/47955845
allArgs <- function(origValues = FALSE)
{
    # get formals for parent function
    parent_formals <- formals(sys.function(sys.parent(n = 1)))

    # Get names of implied arguments
    fnames <- names(parent_formals)

    # Remove '...' from list of parameter names if it exists
    fnames <- fnames[-which(fnames == '...')]

    # Get currently set values for named variables in the parent frame
    args <- evalq(as.list(environment()), envir = parent.frame())

    # Get the list of variables defined in '...'
    args <- c(args[fnames], evalq(list(...), envir = parent.frame()))

    if(origValues)
    {
        # get default values
        defargs <- as.list(parent_formals)
        defargs <- defargs[unlist(lapply(defargs, FUN = function(x) class(x) != "name"))]
        args[names(defargs)] <- defargs
        setargs <- evalq(as.list(match.call())[-1], envir = parent.frame())
        args[names(setargs)] <- setargs
    }

    return(args)
}

openProgBar <- function(min = 0, max, style = 3, ...)
{
    progOpts <- list(min = min, max = max, style = style, ...)
    progOpts <- modifyList(progOpts, getOption("patRoon.progress.opts", list()))
    return(do.call(txtProgressBar, progOpts))
}

verboseCall <- function(f, a, v) if (v) do.call(f, a) else suppressMessages(invisible(do.call(f, a)))

setListNamesIfPresent <- function(l, o, n)
{
    pr <- o %in% names(l)
    o <- o[pr]; n <- n[pr]
    names(l)[match(o, names(l))] <- n
    return(l)
}

babelConvert <- function(input, inFormat, outFormat, mustWork = TRUE)
{
    # Use batch conversion with a single input/output file. Note that obabel
    # will stop after an error. This can be overidden, however, then it is
    # unclear which entries failed. Hence, this option is not used, and when an
    # error occurs the batch conversion is simply restarted with the subsequent
    # entry.

    inputFile <- tempfile("obabel_inp", fileext = ".txt")
    outputFile <- tempfile("obabel_out", fileext = ".txt")
    doConversion <- function(inp)
    {
        cat(inp, file = inputFile, sep = "\n")
        executeCommand(getCommandWithOptPath("obabel", "obabel"),
                       c(paste0("-i", inFormat), inputFile,
                         paste0("-o", outFormat), "-O", outputFile, "-xw"),
                       stderr = FALSE)
        # each conversion is followed by a tab (why??) and newline. Read line
        # by line and remove tab afterwards.
        ret <- readLines(outputFile)
        return(trimws(ret, which = "right", whitespace = "\t"))
    }

    inputLen <- length(input)
    ret <- character(inputLen)
    curIndex <- 1
    while(TRUE)
    {
        curRange <- seq(curIndex, inputLen)
        out <- doConversion(input[curRange])
        outl <- length(out)

        if (outl > 0)
            ret[seq(curIndex, curIndex + outl - 1)] <- out

        curIndex <- curIndex + outl + 1

        if (curIndex <= inputLen)
        {
            msg <- sprintf("Failed to convert %d ('%s')", curIndex - 1, input[curIndex - 1])
            if (mustWork)
                stop(msg)
            else
                warning(msg)
        }
        else
            break
    }

    return(ret)
}

prepareSuspectList <- function(suspects, adduct, skipInvalid)
{
    # UNDONE: check if/make name column is file safe/unique

    if (is.data.table(suspects))
        suspects <- copy(suspects)
    else
        suspects <- as.data.table(suspects)
    
    suspects[, name := as.character(name)] # in case factors are given
    
    checkValidColumn <- function(column, check, printLine = FALSE)
    {
        if (!skipInvalid)
            return(suspects) # checks should've been done by prior assert
        
        valid <- check(suspects[[column]])
        if (any(!valid))
        {
            warning(paste0("Ignored following suspects with invalid data for column '", column, "': ",
                           if (printLine) sprintf("%s (line %d)", suspects$name[!valid], which(!valid)) else suspects$name[!valid]),
                    call. = FALSE)
            suspects <- suspects[valid]
        }
        return(suspects)
    }
    checkValidNumColumn <- function(column, ...) checkValidColumn(column, function(x) !is.na(x), ...)
    checkValidStrColumn <- function(column, ...) checkValidColumn(column, function(x) !is.na(x) & nzchar(x), ...)
    
    # only print line here, next ones might be invalid as the table may have
    # been subset by subsequent calls to this function
    suspects <- checkValidStrColumn("name", printLine = TRUE)
    
    if (!is.null(suspects[["mz"]]))
    {
        suspects <- checkValidNumColumn("mz")
        return(suspects) # no further need for calculation of ion masses
    }

    # neutral mass given?
    if (!is.null(suspects[["neutralMass"]]))
    {
        suspects <- checkValidNumColumn("neutralMass")
        neutralMasses <- suspects[["neutralMass"]]
    }
    else
    {
        # otherwise calculate
        
        printf("Calculating ion masses for each suspect...\n")
        prog <- openProgBar(0, nrow(suspects))
        
        if (!is.null(suspects[["formula"]]))
        {
            suspects <- checkValidStrColumn("formula")
            neutralMasses <- sapply(seq_len(nrow(suspects)), function(i)
            {
                ret <- rcdk::get.formula(suspects$formula[i])@mass
                setTxtProgressBar(prog, i)
                return(ret)
            })
        }
        else
        {
            if (!is.null(suspects[["SMILES"]]))
            {
                suspects <- checkValidStrColumn("SMILES")
                SMI <- suspects$SMILES
            }
            else
            {
                suspects <- checkValidStrColumn("InChI")
                SMI <- babelConvert(suspects$InChI, "inchi", "smi")
            }
            
            neutralMasses <- sapply(seq_along(SMI), function(i)
            {
                ret <- getNeutralMassFromSMILES(SMI[i], mustWork = !skipInvalid)[[1]]
                setTxtProgressBar(prog, i)
                return(ret)
            })
        }
        
        close(prog)
    }
    
    if (!is.null(adduct))
        addMZs <- adductMZDelta(adduct)
    else
        addMZs <- sapply(suspects[["adduct"]], function(a) adductMZDelta(as.adduct(a)))
    
    if (skipInvalid)
    {
        # skip any suspects without proper mass info
        isNA <- is.na(neutralMasses)
        if (any(isNA))
        {
            warning(paste("Ignored following suspects for which no mass could be calculated:",
                          paste0(suspects[["name"]][isNA], collapse = ", ")))
            suspects <- suspects[!isNA]; neutralMasses <- neutralMasses[!isNA]
        }
    }
    
    suspects[, mz := neutralMasses + addMZs][]

    return(suspects)
}
