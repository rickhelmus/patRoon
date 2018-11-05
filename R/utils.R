#' @include main.R
NULL

printf <- function(...) cat(sprintf(...), sep = "")

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
    if (!dir.exists(path))
        dir.create(path, recursive = TRUE)
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

checkPackage <- function(pkg, gh = NULL)
{
    # from http://stackoverflow.com/a/20333756
    if (!requireNamespace(pkg, quietly = TRUE))
    {
        if (!is.null(gh))
            stop(sprintf("Please install %s from github: devtools::install_github('%s')", pkg, gh))
        else
            stop(sprintf("Please install %s: install.packages('%s')", pkg, pkg))
    }
}

executeCommand <- function(cmd, args = character(), ...)
{
    return(system2(cmd, sapply(args, shQuote), ...))
}

getCommandWithOptPath <- function(cmd, opt)
{
    if (Sys.info()[["sysname"]] == "Windows")
        cmd <- paste0(cmd, ".exe") # add file extension for Windows

    opt <- paste0("patRoon.path.", opt)
    path <- getOption(opt)
    if (!is.null(path) && nzchar(path))
    {
        ret <- file.path(path.expand(path), cmd)
        if (!file.exists(ret))
            stop(sprintf("Cannot find '%s'. Is the option '%s' set correctly?", ret, opt))
        return(ret)
    }

    # assume command is in PATH --> no need to add path
    if (!nzchar(Sys.which(cmd)))
        stop(sprintf("Cannot find '%s'. Either add the correct file location to the PATH environment variable or set '%s' with option().", cmd, opt))

    return(cmd)
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

#' @details \code{generateAnalysisInfo} is an utility function that automatically
#'   generates an analysis information object. It will collect all datafiles
#'   from given file paths and convert the filenames into valid analysis names
#'   (\emph{i.e.} without extensions such as \file{.d} and \file{.mzML}).
#'   Duplicate analyses, which may appear when datafiles with different file
#'   extension (\file{.d}, \file{.mzXML} and/or \file{.mzML}) are present, will
#'   be automatically removed.
#'
#' @param paths A character vector containing one or more file paths that should
#'   be used for finding the analyses.
#' @param groups,refs An (optional) character vector containing replicate groups
#'   and references, respectively (will be recycled). If \code{groups} is an
#'   empty character string (\code{""}) the analysis name will be set as
#'   replicate group.
#' @param fileTypes A character vector of analyses file types. Valid values are:
#'   \code{Bruker}, \code{mzXML} and \code{mzML}.
#'
#' @rdname analysis-information
#' @export
generateAnalysisInfo <- function(paths, groups = "", refs = "", fileTypes = c("Bruker", "mzXML", "mzML"))
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertDirectoryExists(paths, access = "r", add = ac)
    checkmate::assertCharacter(groups, min.len = 1, add = ac)
    checkmate::assertCharacter(refs, min.len = 1, add = ac)
    checkmate::assertSubset(fileTypes, c("Bruker", "mzXML", "mzML"), empty.ok = FALSE, add = ac)
    checkmate::reportAssertions(ac)
    
    fileExts <- c(Bruker = ".d", mzXML = ".mzXML", mzML = ".mzML")
    fileTypes <- fileExts[fileTypes] # rename to file extensions

    # escape dots for regex
    exts <- gsub(".", "\\.", fileTypes, fixed = TRUE)
    files <- as.vector(unlist(Vectorize(list.files)(paths, exts, full.names = TRUE)))

    if (length(files) == 0)
    {
        warning(sprintf("No analyses found in %s!", paste(paths, collapse = ", ")))
        return(NULL)
    }

    if (".d" %in% fileTypes)
    {
        # filter out directories unless they end with .d
        files <- files[grepl("(\\.d)$", files) | !file.info(files, extra_cols = FALSE)$isdir]

        # filter out any non directory files that end with .d
        files <- files[!grepl("(\\.d)$", files) | file.info(files, extra_cols = FALSE)$isdir]
    }
    else
        files <- files[!file.info(files, extra_cols = FALSE)$isdir] # filter out directories

    ret <- data.frame(path = dirname(files), analysis = simplifyAnalysisNames(files), group = groups, ref = refs, stringsAsFactors = FALSE)
    ret <- ret[!duplicated(ret[, c("path", "analysis")]), ]
    ret$group <- ifelse(!nzchar(ret$group), ret$analysis, ret$group)

    return(ret)
}

#' @details \code{generateAnalysisInfoFromEnviMass} loads analysis information
#'   from an \pkg{enviMass} project.
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
        blanks[, ref := paste0("blank", match(DT, unique(DT)))] # add unique date+time identifier
        enviSInfo[, ref := sapply(DT, function(dt)
        {
            bls <- blanks[DT <= dt]
            if (nrow(bls) > 0)
                return(bls[which.max(bls$DT), ref])
            else
                return("")
        })]

        enviSInfo[Type == "blank", group := ref]
    }
    else
        enviSInfo[, ref := ""]

    enviSInfo[group == "FALSE", group := ""]

    ret <- data.frame(path = file.path(path, "files"), analysis = enviSInfo$ID, group = enviSInfo$group,
                      ref = enviSInfo$ref, stringsAsFactors = FALSE)

    return(ret)
}

deIsotopeMSPeakList <- function(MSPeakList)
{
    if (nrow(MSPeakList) == 0)
        return(MSPeakList)

    # make sure most intense ions top the table
    MSPeakList <- MSPeakList[order(mz, -intensity)]

    unique_iso <- sapply(seq_along(MSPeakList$cmp), function(i)
    {
        # first and unassigned compounds are always unique
        if (i == 1 || MSPeakList$cmp[i] == "")
            return(TRUE)

        # peak may belong to multiple isotope compounds (separated by whitespace)
        cmp <- strsplit(MSPeakList$cmp[i], "\\s+")

        # isotope compound present before this entry?
        othercmp <- MSPeakList[seq_len(i - 1)][cmp != ""]$cmp
        for (ocmp in othercmp)
        {
            if (any(cmp %in% strsplit(ocmp, "\\s+")))
                return(FALSE)
        }

        return(TRUE)
    }, USE.NAMES = FALSE)

    return(MSPeakList[unique_iso])
}

# from: http://www.cureffi.org/2013/09/23/a-quick-intro-to-chemical-informatics-in-r/
rcdkplot <- function(molecule, width = 500, height = 500)
{
    # drawing unconnected is not supported. Assume these are salts and we can just draw largest instead
    # UNDONE: is this still relevant?
    # if (!is.connected(molecule))
    #     molecule <- get.largest.component(molecule)

    temp <- rcdk::view.image.2d(molecule, rcdk::get.depictor(width, height)) # get Java representation into an image matrix.
    # plot(NA, NA, xlim = c(1, 10), ylim = c(1, 10), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') # create an empty plot
    # rasterImage(temp, 1, 1, 10, 10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
    # grid::grid.newpage()
    # grid::grid.raster(temp)

    img <- magick::image_read(temp)
    # trim image before plotting it
    if (!isEmptyMol(molecule))
        img <- magick::image_transparent(magick::image_trim(img), "white")
    plot(img)
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

# get corresponding mz of feature from MS peaklist
getMZFromMSPeakList <- function(featMZ, plist) return(plist$mz[which.min(abs(plist$mz - featMZ))])

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

stripDTRef <- function(dt)
{
    setattr(dt, ".internal.selfref", NULL)
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
    maxn <- brewer.pal.info[name, "maxcolors"]
    if (n > maxn)
        return(colorRampPalette(brewer.pal(maxn, name))(n))

    if (n < 3)
        return(brewer.pal(3, name)[seq_len(n)])

    return(brewer.pal(n, name))
}

#' Conversion of analyses between several open formats
#'
#' This function uses the \command{FileConverter} command from
#' \href{http://www.openms.de/}{OpenMS} to convert analyses files from one open
#' format to another. The supported file formats are (both directions):
#' \file{.mzXML}, \file{.mzML} and \file{.mzData}.
#'
#' @param anaInfo A table with \link[=analysis-information]{analysis
#'   information} that contains the analyses to be converted.
#' @param formatFrom,formatTo A \code{character} that specifies the source and
#'   destination format, respectively. Valid options are: \code{"mzXML"},
#'   \code{"mzML"} and \code{"mzData"}.
#' @param outPath A character vector specifying directories that should be used for the output. Usually
#'   this is the same as the source (otherwise the
#'   \link[=analysis-information]{analysis information} should be changed for
#'   further processing). Will be re-cycled if necessary.
#' @param overWrite Should existing destination file be overwritten
#'   (\code{TRUE}) or not (\code{FALSE})?
#'
#' @references \insertRef{Rst2016}{patRoon}
#' 
#' @export
convertMSFiles <- function(anaInfo, formatFrom, formatTo, outPath = anaInfo$path, overWrite = FALSE)
{
    ac <- checkmate::makeAssertCollection()
    assertAnalysisInfo(anaInfo, add = ac)
    checkmate::assertChoice(formatFrom, c("mzXML", "mzML", "mzData"), add = ac)
    checkmate::assertChoice(formatTo, c("mzXML", "mzML", "mzData"), add = ac)
    checkmate::assertCharacter(outPath, min.chars = 1, min.len = 1, add = ac)
    checkmate::assertDirectoryExists(outPath, "w", add = ac)
    checkmate::assertFlag(overWrite, add = ac)
    checkmate::reportAssertions(ac)
    
    outPath <- rep(outPath, length.out = length(anaInfo$path))

    for (anai in seq_len(nrow(anaInfo)))
    {
        input <- switch(formatFrom,
                        mzXML = getMzXMLAnalysisPath(anaInfo$analysis[anai], anaInfo$path[anai]),
                        mzML = getMzMLAnalysisPath(anaInfo$analysis[anai], anaInfo$path[anai]),
                        mzData = getMzDataAnalysisPath(anaInfo$analysis[anai], anaInfo$path[anai]),
                        stop(sprintf("Unknown input file format: %s!", formatFrom)))

        outExt <- switch(formatTo,
                         mzXML = ".mzXML",
                         mzML = ".mzML",
                         mzData = ".mzData",
                         stop(sprintf("Unknown output file format: %s!", formatTo)))

        output <- file.path(outPath[anai], paste0(anaInfo$analysis[anai], outExt))

        if (!file.exists(input))
            warning(sprintf("Skipping non-existing input analysis %s", input))
        else if (overWrite || !file.exists(output))
        {
            if (executeCommand(getCommandWithOptPath("FileConverter", "OpenMS"),
                               c("-in", input, "-out", output)) != 0)
                stop("Failed to execute OpenMS FileConverter command.")
        }
    }
}

getArgNames <- function(...)
{
    args <- sapply(substitute(list(...))[-1], deparse)
    n <- names(args)
    if (!is.null(n))
        args <- ifelse(n != "", n, args)
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
    refGroups <- unique(anaInfo$ref)
    printf("Analyses: %s (%d total)\n", getStrListWithMax(anaInfo$analysis, 6, ", "), nrow(anaInfo))
    printf("Replicate groups: %s (%d total)\n", getStrListWithMax(rGroups, 8, ", "), length(rGroups))
    printf("Replicate groups used as ref: %s (%d total)\n", getStrListWithMax(refGroups, 8, ", "), length(refGroups))
}

showObjectSize <- function(object) printf("Object size (indication): %s\n", format(object.size(object), "auto", "SI"))

allSame <- function(l)
{
    if (length(l) > 0)
    {
        if (all(is.na(l)))
            return(TRUE)
        if (any(is.na(l)))
            return(FALSE)

        return(all(l[[1]] == l))
    }

    return(TRUE)
}

normalize <- function(x)
{
    xn <- x[!is.na(x)]
    if (length(xn) == 0 || all(xn == 0))
        return(x) # all NA or all values zero

    if (allSame(xn))
        xn <- rep(as(1, typeof(xn[[1]])), length(xn))
    else
    {
        xn <- xn - min(xn)
        xn <- xn / max(xn)
    }

    x[!is.na(x)] <- xn

    return(x)
}

# from https://stackoverflow.com/a/38228840
countCharInStr <- function(str, ch) sum(charToRaw(str) == charToRaw(ch))

makeMSPlot <- function(spec, fragInfo, ..., extraHeightInch = 0)
{
    isMerged <- !is.null(fragInfo[["mergedBy"]]) && nrow(fragInfo) > 0
    oldp <- NULL
    # par(mar = c(par("mar")[1:2], 0, 0))

    if (isMerged)
    {
        allMergedBy <- sapply(fragInfo[["mergedBy"]], function(mb) wrapStr(paste0(unlist(mb), collapse = ", "), 10))
        mbsUnique <- unique(allMergedBy)
        # order from small to big based on number of commas
        mbsUnique <- mbsUnique[order(sapply(mbsUnique, countCharInStr, ch = ",", USE.NAMES = FALSE))]
        mbCombCols <- getBrewerPal(length(mbsUnique), "Paired")

        makeLegend <- function(x, y, ...) legend(x, y, mbsUnique, col = mbCombCols, xpd = NA, bty = "n",
                                                 text.col = mbCombCols, lty = 1, cex = 0.75, ...)

        plot.new()
        leg <- makeLegend(0, 0, plot = FALSE)
        lw <- (grconvertX(leg$rect$w, to = "ndc") - grconvertX(0, to = "ndc"))
        oldp <- par(omd = c(0, 1 - lw, 0, 1), new = TRUE)
    }

    # ym <- (if (is.null(fragInfo)) max(spec$intensity) else max(fragInfo$intensity)) * 1.2
    fragInfoInds <- if (nrow(fragInfo) > 0) match(seq_len(nrow(spec)), fragInfo$PLIndex, 0) else rep(0, nrow(spec))
    
    # see how much extra vertical space is needed by formula labels
    # get character widths (assuming that height of vertically plotted text is the same)
    formWidths <- sapply(fragInfoInds, function(i) if (i) strwidth(fragInfo$formula[i], units = "inches") else 0)
    
    pheight <- par("din")[2] # plot dev height in inches
    relFormHeights <- formWidths / pheight
    
    ym <- max(spec$intensity) # 'regular' plot height in user coordinates, not considering any labels etc
    relIntHeights <- spec$intensity / ym
    maxRelH <- max(relFormHeights + relIntHeights)
    
    ym <- ym * maxRelH * 1.05 # enlarge y limit and add some extra spacing
    
    if (extraHeightInch > 0)
        ym <- ym * (1 + (extraHeightInch / pheight))
    
    plot(0, xlab = "m/z", ylab = "Intensity", xlim = range(spec$mz) * c(0.9, 1.1), ylim = c(0, ym),
         type = "n", bty = "l", ...)

    for (i in seq_len(nrow(spec)))
    {
        infoi <- fragInfoInds[i]

        if (infoi && isMerged)
            specCol <- mbCombCols[match(allMergedBy[infoi], mbsUnique)]
        else if (infoi)
            specCol <- "blue"
        else
            specCol <- "grey"

        segments(spec[i]$mz, 0, spec[i]$mz, spec[i]$intensity, col = specCol, lwd = if (infoi) 2 else 1)

        if (infoi)
            text(spec[i]$mz, spec[i]$intensity + (ym * 0.02), fragInfo$formula[infoi], srt = 90, adj = 0)
    }

    if (isMerged)
        makeLegend(par("usr")[2], par("usr")[4])

    if (!is.null(oldp))
        par(oldp)
}

makeMSPlotGG <- function(spec, fragInfo, ...)
{
    if (!is.null(fragInfo[["mergedBy"]]))
    {
        allMergedBy <- sapply(fragInfo[["mergedBy"]], function(mb) paste0(unlist(mb), collapse = ","))
        mbsUnique <- unique(allMergedBy)
        # order from small to big based on number of commas
        mbsUnique <- mbsUnique[order(sapply(mbsUnique, countCharInStr, ch = ",", USE.NAMES = FALSE))]
        mbCombCols <- getBrewerPal(length(mbsUnique), "Paired")
    }

    if (is.null(fragInfo) || nrow(fragInfo) == 0)
        fragSpecIndices <- NA
    else
        fragSpecIndices <- sapply(seq_len(nrow(spec)), function(r) match(r, fragInfo$PLIndex))

    plotData <- copy(spec)

    plotData[, c("colour", "lab", "lwd", "text") := .("grey", "unassigned", 0.5, "")]
    if (!is.null(fragInfo) && nrow(fragInfo) > 0)
    {
        plotData[, fiInd := sapply(seq_len(nrow(spec)), function(r) match(r, fragInfo$PLIndex))]
        if (!is.null(fragInfo[["mergedBy"]]))
        {
            plotData[!is.na(fiInd), colour := mbCombCols[match(allMergedBy[fiInd], mbsUnique)]]
            plotData[!is.na(fiInd), lab := allMergedBy[fiInd]]
        }
        else
        {
            plotData[!is.na(fiInd), colour := "blue"]
            plotData[!is.na(fiInd), lab := "assigned"]
        }

        plotData[!is.na(fiInd), lwd := 2]
        plotData[!is.na(fiInd), text := fragInfo$formula[fiInd]]
    }

    ggplot(plotData, aes_string(x = "mz", y = 0, label = "text")) + xlim(range(spec$mz) * c(0.9, 1.1)) +
        geom_segment(aes_string(xend = "mz", yend = "intensity",
                                colour = "lab", size = "lwd")) + scale_size(range = c(0.5, 2), guide = FALSE) +
        ggrepel::geom_text_repel(aes_string(y = "intensity", angle = 0), min.segment.length = 0.1,
                                 nudge_y = grid::convertUnit(grid::unit(5, "mm"), "npc", valueOnly = TRUE), size = 3.2) +
        xlab("m/z") + ylab("Intensity") +
        cowplot::theme_cowplot(font_size = 12) + theme(legend.position = "bottom", legend.title = element_blank())
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

isValidMol <- function(mol) !is.null(mol) # && !is.na(mol)
emptyMol <- function() rcdk::parse.smiles("")[[1]]
isEmptyMol <- function(mol) rcdk::get.atom.count(mol) == 0

getMoleculesFromSMILES <- function(SMILES, doTyping = FALSE, emptyIfFails = FALSE)
{
    # vectorization doesn't work if any of the SMILES are not OK
    # mols <- rcdk::parse.smiles(SMILES)
    mols <- lapply(SMILES, function(sm)
    {
        ret <- rcdk::parse.smiles(sm)[[1]]
        if (!isValidMol(ret))
            ret <- rcdk::parse.smiles(sm, kekulise = FALSE)[[1]] # might work w/out kekulization
        if (emptyIfFails && !isValidMol(ret))
            ret <- emptyMol()
        else if (doTyping && isValidMol(ret))
        {
            rcdk::do.typing(ret)
            rcdk::do.aromaticity(ret)
        }
        return(ret)
    })
    
    return(mols)
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
numEQ <- function(x, y, tol = sqrt(.Machine$double.eps)) abs(x - y) < tol
numGTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x > y
numLTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x < y

wrapStr <- function(s, width, sep = "\n") paste0(strwrap(s, width), collapse = sep)

pruneList <- function(l, checkEmptyElements = FALSE)
{
    ret <- l[!sapply(l, is.null)]
    if (checkEmptyElements)
        ret <- ret[lengths(ret) > 0]
    return(ret)
}
