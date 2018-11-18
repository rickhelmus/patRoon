#' @include main.R
#' @include components.R
NULL

componentsCamera <- setClass("componentsCamera", slots = c(xsa = "xsAnnotate"), contains = "components")

#' @details \code{generateComponentsCAMERA} provides an interface to
#'   \href{https://bioconductor.org/packages/release/bioc/html/CAMERA.html}{CAMERA}
#'    which is used to generate components from known adducts, isotopes and
#'   in-source fragments. The specified \code{featureGroups} object is
#'   automatically converted to an \code{\link{xcmsSet}} object using
#'   \code{\link{getXCMSSet}}.
#'
#' @param onlyIsotopes Logical value. If \code{TRUE} only isotopes are
#'   considered when generating components (faster). Corresponds to \code{quick}
#'   argument of \code{\link[CAMERA:annotate-methods]{CAMERA::annotate}}.
#'
#' @references \addCitations{CAMERA}{1}
#'
#' @rdname component-generation
#' @export
generateComponentsCAMERA <- function(fGroups, ionization, onlyIsotopes = FALSE, extraOpts = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertChoice(ionization, c("positive", "negative"), add = ac)
    checkmate::assertFlag(onlyIsotopes, add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(componentsCamera(componentInfo = data.table(), components = list(),
                                algorithm = "CAMERA",
                                xsa = new("xsAnnotate")))

    hash <- makeHash(fGroups, ionization, onlyIsotopes, extraOpts)
    cd <- loadCacheData("componentsCAMERA", hash)
    if (!is.null(cd))
        return(cd)

    gNames <- names(fGroups)
    gTable <- groups(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)

    xs <- getXCMSSet(fGroups, TRUE) # UNDONE: handle exportedData: check if all files are present? Is it necessary?

    # UNDONE: do fillpeaks? What about RC?

    mainArgs <- c(list(polarity = ionization, quick = onlyIsotopes), extraOpts)

    an <- do.call(CAMERA::annotate, c(list(object = xs), mainArgs))

    anPList <- as.data.table(CAMERA::getPeaklist(an))
    anPList <- anPList[, c("rt", "mz", "isotopes", "adduct", "pcgroup")]
    anPList[, pcgroup := as.numeric(pcgroup)] # why is this a character? Convert to numeric for sorting.
    anPList[, group := names(fGroups)]

    # get isotope information
    isoTab <- NULL
    isoGrpInds <- anPList[nzchar(isotopes), which = TRUE]
    if (length(isoGrpInds) > 0)
    {
        isoTab <- rbindlist(an@isotopes[isoGrpInds])
        isoTab[, group := gNames[isoGrpInds]]
        setnames(isoTab, c("y", "iso", "charge"), c("isogroup", "isonr", "charge"))
        isoTab[, val := NULL] # UNDONE: what is this?
        anPList[, isotopes := NULL]
    }

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
        anPList[, adduct := NULL]
    }

    # merge isotope/adduct rows: do this afterwards so that uniqueness and order
    # of fGroup rows remain during collection of isotope/adduct information.
    if (!is.null(isoTab))
        anPList <- merge(anPList, isoTab, all = TRUE, by = "group")
    if (!is.null(adTable))
        anPList <- merge(anPList, adTable, all = TRUE, by = "group")

    setorderv(anPList, "pcgroup")
    anPList[, pcgroup := paste0("CMP", pcgroup)]
    comps <- split(anPList, by = "pcgroup", keep.by = FALSE)

    # add intensities. NOTE: CAMERA selects for each pcgroup one representative
    # sample
    cnames <- names(comps)
    comps <- sapply(seq_along(comps), function(cmpi)
    {
        anai <- an@psSamples[[cmpi]]
        comps[[cmpi]][, "intensity" := unlist(gTable[anai, comps[[cmpi]]$group, with = FALSE])]
    }, simplify = FALSE)
    names(comps) <- cnames

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

    cInfo <- data.table(name = names(comps), ret = sapply(rets, mean),
                        retsd = sapply(rets, sd), neutral_mass = Ms,
                        analysis = anaInfo$analysis[unlist(an@psSamples)],
                        size = sapply(comps, nrow))


    # UNDONE: keep n=1 sized components?

    ret <- componentsCamera(xsa = an, components = comps, componentInfo = cInfo, algorithm = "CAMERA")
    saveCacheData("componentsCAMERA", ret, hash)
    return(ret)
}
