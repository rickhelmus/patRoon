# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include main.R
#' @include components.R
#' @include feature_groups-set.R
NULL

#' @rdname components-class
#' @export
componentsCamera <- setClass("componentsCamera", slots = c(xsa = "ANY"),
                             contains = "components")
setMethod("initialize", "componentsCamera",
          function(.Object, ...) callNextMethod(.Object, algorithm = "camera", ...))

#' Componentization of adducts, isotopes etc. with CAMERA
#'
#' Interfaces with \href{https://bioconductor.org/packages/release/bioc/html/CAMERA.html}{CAMERA} to generate components
#' from known adducts, isotopes and in-source fragments.
#'
#' @templateVar algo CAMERA
#' @templateVar do generate components
#' @templateVar generic generateComponents
#' @templateVar algoParam camera
#' @template algo_generator
#'
#' @details The specified \code{featureGroups} object is automatically converted to an \code{\link{xcmsSet}} object
#'   using \code{\link{getXCMSSet}}.
#'
#' @param onlyIsotopes Logical value. If \code{TRUE} only isotopes are considered when generating components (faster).
#'   Corresponds to \code{quick} argument of \code{\link[CAMERA:annotate-methods]{CAMERA::annotate}}.
#' @param extraOpts Named character vector with extra arguments directly passed to
#'   \code{\link[CAMERA:annotate-methods]{CAMERA::annotate}}. Set to \code{NULL} to ignore.
#'
#' @templateVar ion TRUE
#' @templateVar minSize TRUE
#' @templateVar minReps TRUE
#' @template compon_algo-args
#'
#' @inheritParams generateComponents
#'
#' @inherit generateComponents return
#' 
#' @templateVar minSize TRUE
#' @template compon_gen-filters
#' 
#' @templateVar class componentsSet
#' @template compon_gen-sets-merged
#'
#' @references \addCitations{CAMERA}
#' 
#' @templateVar what generateComponentsCAMERA
#' @template main-rd-method
#' @export
setMethod("generateComponentsCAMERA", "featureGroups", function(fGroups, ionization = NULL, onlyIsotopes = FALSE,
                                                                minSize = 2, relMinReplicates = 0.5, extraOpts = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    ionization <- checkAndGetIonization(ionization, fGroups, add = ac)
    checkmate::assertFlag(onlyIsotopes, add = ac)
    checkmate::assertCount(minSize, positive = TRUE, add = ac)
    checkmate::assertNumber(relMinReplicates, lower = 0, finite = TRUE, add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(componentsCamera(componentInfo = data.table(), components = list(),
                                xsa = new("xsAnnotate")))

    hash <- makeHash(fGroups, ionization, onlyIsotopes, minSize, relMinReplicates, extraOpts)
    cd <- loadCacheData("componentsCAMERA", hash)
    if (!is.null(cd))
        return(cd)

    gNames <- names(fGroups)
    gTable <- groupTable(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)

    xs <- getXCMSSet(fGroups, loadRawData = TRUE) # UNDONE: handle loadRawData: check if all files are present? Is it necessary?

    # UNDONE: do fillpeaks? What about RC?

    mainArgs <- c(list(polarity = ionization, quick = onlyIsotopes), extraOpts)

    an <- do.call(CAMERA::annotate, c(list(object = xs), mainArgs))

    anPList <- as.data.table(CAMERA::getPeaklist(an))
    anPList <- anPList[, c("rt", "mz", "isotopes", "adduct", "pcgroup")]
    anPList[, pcgroup := as.numeric(pcgroup)] # why is this a character? Convert to numeric for sorting.
    anPList[, group := names(fGroups)]
    setnames(anPList, "rt", "ret")

    # get isotope information
    isoTab <- NULL
    isoGrpInds <- anPList[nzchar(isotopes), which = TRUE]
    if (length(isoGrpInds) > 0)
    {
        isoTab <- rbindlist(an@isotopes[isoGrpInds])
        isoTab[, group := gNames[isoGrpInds]]
        setnames(isoTab, c("y", "iso", "charge"), c("isogroup", "isonr", "charge"))
        isoTab[, val := NULL] # UNDONE: what is this?
    }
    anPList[, isotopes := NULL]

    # get adduct information
    adTable <- NULL
    adGrpInds <- anPList[nzchar(adduct), which = TRUE]
    if (length(adGrpInds) > 0)
    {
        adTable <- rbindlist(lapply(adGrpInds, function(grpi)
        {
            ret <- rbindlist(an@derivativeIons[[grpi]], idcol = "adnr")
            ret[, group := gNames[grpi]]
            setnames(ret, c("rule_id", "charge", "nmol", "name", "mass"),
                     c("adduct_rule", "adduct_charge", "adduct_nmol", "adduct_ion", "M_adduct"))
            return(ret[!duplicated(ret[, -"adnr"])]) # BUG: sometimes there are duplicate entries?
        }))
    }
    anPList[, adduct := NULL]

    # merge isotope/adduct rows: do this afterwards so that uniqueness and order
    # of fGroup rows remain during collection of isotope/adduct information.
    # If unavailable, NA columns are added to get predictable columns
    if (!is.null(isoTab))
        anPList <- merge(anPList, isoTab, all = TRUE, by = "group")
    else
        anPList[, c("isogroup", "isonr", "charge") := NA_real_]
    if (!is.null(adTable))
        anPList <- merge(anPList, adTable, all = TRUE, by = "group")
    else
        anPList[, c("adnr", "adduct_rule", "adduct_charge", "adduct_nmol", "adduct_ion", "M_adduct") := 
                    .(NA_integer_, NA_real_, NA_real_, NA_real_, NA_character_, NA_real_)]
    
    setorderv(anPList, "pcgroup")
    anPList[, pcgroup := paste0("CMP", pcgroup)]
    comps <- split(anPList, by = "pcgroup", keep.by = FALSE)

    # store before filtering/pruning below
    cmpanalyses <- setNames(anaInfo$analysis[unlist(an@psSamples)], names(comps))
    
    # Filter non-abundant feature groups from components and add intensities.
    # NOTE: CAMERA selects for each pcgroup one representative sample
    cnames <- names(comps)
    comps <- pruneList(setNames(sapply(seq_along(comps), function(cmpi)
    {
        cmp <- comps[[cmpi]]
        if (relMinReplicates > 0)
        {
            fgCmp <- removeEmptyAnalyses(fGroups[, cmp$group])
            fgCmp <- minReplicatesFilter(fgCmp, relThreshold = relMinReplicates, verbose = FALSE)
            cmp <- cmp[group %in% names(fgCmp)]
        }
        if (nrow(cmp) < minSize)
            return(cmp[0])
            
        anai <- an@psSamples[[cmpi]]
        cmp[, "intensity" := unlist(gTable[anai, cmp$group, with = FALSE])]
    }, simplify = FALSE), cnames), checkZeroRows = TRUE)
    
    if (length(comps) != length(cmpanalyses))
    {
        # component(s) got filtered out, update
        cmpanalyses <- cmpanalyses[names(comps)] # update
        if (length(comps) > 0)
            names(comps) <- paste0("CMP", seq_along(comps))
    }
    
    rets <- lapply(comps, function(cm) gInfo[cm$group, "rts"])
    if (!is.null(adTable))
    {
        Ms <- sapply(comps, function(cm)
        {
            M <- cm[!is.na(M_adduct), M_adduct]
            if (length(M) > 0)
                return(paste0(unique(round(M, 5)), collapse = "/"))
            return(NA)
        })
    }
    else
        Ms <- NA

    cInfo <- data.table(name = names(comps), cmp_ret = sapply(rets, mean),
                        cmp_retsd = sapply(rets, sd), neutral_mass = Ms,
                        analysis = cmpanalyses, size = sapply(comps, nrow))


    ret <- componentsCamera(xsa = an, components = comps, componentInfo = cInfo)
    saveCacheData("componentsCAMERA", ret, hash)
    return(ret)
})

#' @rdname generateComponentsCAMERA
#' @export
setMethod("generateComponentsCAMERA", "featureGroupsSet", function(fGroups, ionization = NULL, ...)
{
    generateComponentsSet(fGroups, ionization, generateComponentsCAMERA, setIonization = TRUE, ...)
})

