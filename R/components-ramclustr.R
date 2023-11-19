#' @include main.R
#' @include components.R
#' @include feature_groups-set.R
NULL

#' @rdname components-class
componentsRC <- setClass("componentsRC", slots = c(RC = "hclust"),
                         contains = "components")
setMethod("initialize", "componentsRC",
          function(.Object, ...) callNextMethod(.Object, ..., algorithm = "ramclustr"))


#' Componentization of adducts, isotopes etc. with RAMClustR
#'
#' Uses \href{https://github.com/cbroeckl/RAMClustR}{RAMClustR} to generate components from feature groups which follow
#' similar chromatographic retention profiles and annotate their relationships (\emph{e.g.} adducts and isotopes).
#'
#' @templateVar algo RAMClustR
#' @templateVar do generate components
#' @templateVar generic generateComponents
#' @templateVar algoParam ramclustr
#' @template algo_generator
#'
#' @details This method uses the \code{\link[RAMClustR]{ramclustR}} functions for generating the components, whereas
#'   \code{\link[RAMClustR]{do.findmain}} is used for annotation.
#'
#' @param st,sr,maxt,hmax,normalize Arguments to tune the behaviour of feature group clustering. See their documentation
#'   from \code{\link[RAMClustR]{ramclustR}}. When \code{st} is \code{NULL} it will be automatically calculated as the
#'   half of the median for all chromatographic peak widths.
#' @param relMzDev Maximum relative mass deviation (\acronym{ppm}). Sets the \code{ppm.error} argument to
#'   \code{\link[RAMClustR]{do.findmain}}.
#' @param RCExperimentVals A named \code{list} containing two more \code{list}s: \code{design} and \code{instrument}.
#'   These are used to construct the \code{ExpDes} argument passed to \code{\link[RAMClustR]{ramclustR}}.
#' @param extraOptsRC,extraOptsFM Named \code{list} with further arguments to be passed to
#'   \code{\link[RAMClustR]{ramclustR}} and \code{\link[RAMClustR]{do.findmain}}. Set to \code{NULL} to ignore.
#'
#' @templateVar ion TRUE
#' @templateVar RC TRUE
#' @templateVar minSize TRUE
#' @templateVar minReps TRUE
#' @templateVar absMzDev \code{mzabs.error} argument to \code{\link[RAMClustR]{do.findmain}}
#' @template compon_algo-args
#'
#' @inheritParams generateComponents
#'
#' @inherit generateComponents return
#'
#' @templateVar minSize FALSE
#' @template compon_gen-filters
#'
#' @templateVar class componentsSet
#' @template compon_gen-sets-merged
#'
#' @references \insertRef{Broeckling2013}{patRoon} \cr\cr \insertRef{Broeckling2014}{patRoon}
#'
#' @templateVar what generateComponentsRAMClustR
#' @template main-rd-method
#' @export
setMethod("generateComponentsRAMClustR", "featureGroups", function(fGroups, ionization = NULL, st = NULL, sr = NULL,
                                                                   maxt = 12, hmax = 0.3, normalize = "TIC",
                                                                   absMzDev = 0.002, relMzDev = 5,
                                                                   minSize = 2, relMinReplicates = 0.5,
                                                                   RCExperimentVals = list(design = list(platform = "LC-MS"),
                                                                                           instrument = list(ionization = ionization, MSlevs = 1)),
                                                                   extraOptsRC = NULL, extraOptsFM = NULL)
{
    checkPackage("RAMClustR", "cbroeckl/RAMClustR")
    
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    aapply(checkmate::assertNumber, . ~ st + sr + maxt + hmax, lower = 0, null.ok = TRUE, fixed = list(add = ac))
    checkmate::assertString(normalize, add = ac)
    ionization <- checkAndGetIonization(ionization, fGroups, add = ac)
    checkmate::assertNumber(absMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(relMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::assertCount(minSize, positive = TRUE, add = ac)
    checkmate::assertNumber(relMinReplicates, lower = 0, finite = TRUE, add = ac)
    checkmate::assertList(RCExperimentVals, any.missing = FALSE, names = "unique", add = ac)
    checkmate::assertList(extraOptsRC, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::assertList(extraOptsFM, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(componentsRC(componentInfo = data.table(), components = list(),
                            RC = structure(list(), class = "hclust")))

    gTable <- groupTable(fGroups)
    gInfo <- groupInfo(fGroups)
    gNames <- names(fGroups)
    ftinds <- groupFeatIndex(fGroups)
    fTable <- featureTable(fGroups)
    anaInfo <- analysisInfo(fGroups)

    RCMainArgs <- c(list(st = st, sr = sr, maxt = maxt, hmax = hmax, normalize = normalize,
                         minModuleSize = minSize), extraOptsRC)
    FMMainArgs <- c(list(mode = ionization, mzabs.error = absMzDev, ppm.error = relMzDev), extraOptsFM)

    hash <- makeHash(fGroups, minSize, relMinReplicates, RCMainArgs, FMMainArgs)
    cd <- loadCacheData("componentsRC", hash)
    if (!is.null(cd))
        return(cd)

    # RAMClustR needs csv for input
    fGRoupsCSVTable <- as.data.frame(gTable)
    rownames(fGRoupsCSVTable) <- anaInfo$analysis
    colnames(fGRoupsCSVTable) <- paste(gInfo$mzs, gInfo$rts, sep = "_")
    csvFile <- tempfile("fGroups", fileext = ".csv")
    write.table(fGRoupsCSVTable, csvFile, row.names = TRUE, col.names = NA, sep = ",")

    exp <- RAMClustR::defineExperiment(force.skip = TRUE)
    if (!is.null(RCExperimentVals$design))
        exp$design <- modifyList(exp$design, RCExperimentVals$design)
    if (!is.null(RCExperimentVals$instrument))
        exp$instrument <- modifyList(exp$instrument, RCExperimentVals$instrument)

    printf("Generating components with RAMClustR...\n")

    if (is.null(RCMainArgs$st))
    {
        # Get half of median chrom peak width
        # (as is done in ramclustR function when XCMS object is used)

        doAna <- seq_len(nrow(anaInfo))
        doAna <- doAna[sapply(doAna, function(anai) any(ftinds[anai] != 0))] # skip analyses without features

        presentFTIds <- lapply(doAna, function(anai)
        {
            inds <- unlist(ftinds[anai])
            return(inds[inds != 0])
        })
        names(presentFTIds) <- anaInfo$analysis[doAna]

        pWidths <- unlist(sapply(anaInfo$analysis[doAna],
                                 function(ana) fTable[[ana]][presentFTIds[[ana]], retmax - retmin]))

        RCMainArgs$st <- median(pWidths / 2)
    }

    RC <- do.call(RAMClustR::ramclustR, c(list(ms = csvFile, ExpDes = exp, mspout = FALSE), RCMainArgs))

    printf("Annotating components...\n")

    # UNDONE: make optional?
    RC <- do.call(RAMClustR::do.findmain, c(list(ramclustObj = RC, plot.findmain = FALSE, writeMat = FALSE,
                                                 writeMS = FALSE), FMMainArgs))

    comps <- lapply(RC$M.ann, as.data.table)

    # Order of pseudo spectra is changed and can be retrieved from rownames
    psorder <- lapply(RC$M.ann, function(a) as.integer(rownames(a)))
    comps <- lapply(seq_along(comps), function(cmi)
    {
        # label column seems to be the same as adduct column, remove it
        comps[[cmi]][, label := NULL]

        # link feature groups: original order is in RC$xcmsOrd
        comps[[cmi]][, "group" := gNames[RC$xcmsOrd[RC$featclus == cmi][psorder[[cmi]]]]]

        comps[[cmi]][, rt := gInfo[group, "rts"]]
        comps[[cmi]][, intensity := RC$msint[RC$featclus == cmi][psorder[[cmi]]]]
        setnames(comps[[cmi]],
                 c("rt", "int", "isogr", "iso", "adduct"),
                 c("ret", "intensity_rel", "isogroup", "isonr", "adduct_ion"))
        setcolorder(comps[[cmi]], c("ret", "mz", "intensity", "intensity_rel", "group",
                                    "isogroup", "isonr", "charge", "adduct_ion", "ppm"))
        
        # may be converted to logical if all NA
        comps[[cmi]][, c("isogroup", "isonr", "charge") := .(as.integer(isogroup), as.integer(isonr),
                                                             as.integer(charge))]
        comps[[cmi]][, adduct_ion := as.character(adduct_ion)]
        comps[[cmi]][, ppm := as.numeric(ppm)]
    })
    names(comps) <- paste0("CMP", seq_along(comps))
    
    # seems the overall ppm (M.ppm) value was forgotten, make it here
    Mppm <- RC$M.ppm.findmain
    Mppm[!RC$use.findmain] <- RC$M.ppm.ramclustr[!RC$use.findmain]
    
    # UNDONE: include both main+ramclust results?
    cInfo <- data.table(name = names(comps), cmp_ret = RC$clrt, cmp_retsd = RC$clrtsd,
                        neutral_mass = RC$M, cmp_ppm = Mppm, size = sapply(comps, nrow))
    
    # filter components if necessary (do this afterwards as order had to be retained)
    if (relMinReplicates > 0)
    {
        comps <- pruneList(lapply(comps, function(cmp)
        {
            fgCmp <- removeEmptyAnalyses(fGroups[, cmp$group])
            fgCmp <- minReplicatesFilter(fgCmp, relThreshold = relMinReplicates, verbose = FALSE)
            return(cmp[group %in% names(fgCmp)])
        }), checkZeroRows = TRUE)
        
        if (length(comps) != nrow(cInfo))
        {
            # update if components were filtered
            cInfo <- cInfo[name %in% names(comps)]
            if (length(comps) > 0)
            {
                names(comps) <- paste0("CMP", seq_along(comps))
                cInfo[, name := names(comps)]
            }
        }
    }

    ret <- componentsRC(RC = RC, components = comps, componentInfo = cInfo)
    saveCacheData("componentsRC", ret, hash)
    return(ret)
})

#' @rdname generateComponentsRAMClustR
#' @export
setMethod("generateComponentsRAMClustR", "featureGroupsSet", function(fGroups, ionization = NULL, ...)
{
    generateComponentsSet(fGroups, ionization, generateComponentsRAMClustR, setIonization = TRUE, ...)
})
