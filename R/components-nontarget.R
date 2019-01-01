#' @include main.R
#' @include components.R
NULL

#' @details \code{generateComponentsNontarget} uses
#'   \href{https://cran.r-project.org/web/packages/nontarget/index.html}{the
#'   nontarget R package} to generate components by unsupervised detection of
#'   homologous series. In the first step the \code{\link{homol.search}}
#'   function is used to detect all homologues within each replicate group
#'   (analyses within each replicate group are averaged prior to detection).
#'   Then, homologous series across replicate groups are merged in case of full
#'   overlap or when merging of partial overlapping series causes no conflicts.
#'
#' @param rtRange A numeric vector containing the minimum and maximum retention
#'   time (in seconds) between homologues. Series are always considered from low
#'   to high \emph{m/z}, thus, a negative minimum retention time allows
#'   detection of homologous series with increasing \emph{m/z} and decreasing
#'   retention times. These values set the \code{minrt} and \code{maxrt}
#'   arguments of \code{\link{homol.search}}.
#' @param mzRange A numeric vector specifying the minimum and maximum \emph{m/z}
#'   increment of a homologous series. Sets the \code{minmz} and \code{maxmz}
#'   arguments of \code{\link{homol.search}}.
#' @param elements A character vector with elements to be considered for
#'   detection of repeating units. Sets the \code{elements} argument of
#'   \code{\link{homol.search}} function.
#' @param maxRTDev,maxMzDev Maximum retention and (absolute) \emph{m/z}
#'   deviation for detection of homologues within series. These arguments set
#'   the \code{rttol} and \code{mztol} arguments of \code{\link{homol.search}}.
#'
#' @references \addCitations{nontarget}{1} \cr\cr
#'   \addCitations{enviPat}{1}
#'
#' @rdname component-generation
#' @export
generateComponentsNontarget <- function(fGroups, ionization, rtRange = c(-120, 120), mzRange = c(5, 120),
                                        elements = c("C", "H", "O"), maxRTDev = 30, maxMzDev = 0.002,
                                        extraOpts = NULL)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertClass(fGroups, "featureGroups", add = ac)
    checkmate::assertChoice(ionization, c("positive", "negative"), add = ac)
    checkmate::assertNumeric(rtRange, finite = TRUE, any.missing = FALSE, len = 2, add = ac)
    checkmate::assertNumeric(mzRange, lower = 0, finite = TRUE, any.missing = FALSE, len = 2, add = ac)
    checkmate::assertCharacter(elements, min.chars = 1, any.missing = FALSE, min.len = 1, add = ac)
    checkmate::assertNumber(maxRTDev, lower = 0, finite = TRUE, add = ac)
    checkmate::assertNumber(maxMzDev, lower = 0, finite = TRUE, add = ac)
    checkmate::assertList(extraOpts, any.missing = FALSE, names = "unique", null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (length(fGroups) == 0)
        return(components(componentInfo = data.table(), components = list(), algorithm = "nontarget"))

    hash <- makeHash(fGroups, ionization, rtRange, mzRange, elements, maxRTDev, maxMzDev, extraOpts)
    cd <- loadCacheData("componentsNontarget", hash)
    if (!is.null(cd))
        return(cd)

    gTable <- groups(fGroups)
    gInfo <- groupInfo(fGroups)
    anaInfo <- analysisInfo(fGroups)

    data(isotopes, package = "enviPat")

    # adduct.search() and homol.search() call stop() when nothing is found (arg!) --> use tryCatch

    # For now just stick with homologues as there is no easy way to do adduct/isotopes over multiple analyses
    # find homologous series for each replicate group
    rGroups <- unique(anaInfo$group)
    groupTablesRG <- sapply(rGroups, function(rg)
    {
        fGrpRep <- replicateGroupFilter(fGroups, rg, verbose = FALSE)
        return(if (length(fGrpRep) == 0) NULL else as.data.table(fGrpRep, average = TRUE))
    }, simplify = FALSE)

    homArgs <- list(isotopes = isotopes, elements = elements, minmz = mzRange[1],
                     maxmz = mzRange[2], minrt = rtRange[1], maxrt = rtRange[2], ppm = FALSE,
                     mztol = maxMzDev, rttol = maxRTDev)
    if (!is.null(extraOpts))
        homArgs <- modifyList(homArgs, extraOpts)

    homList <- sapply(rGroups, function(rg)
    {
        gt <- groupTablesRG[[rg]][, c("mz", rg, "ret"), with = FALSE] # convert to nontarget peaklist format

        if (is.null(gt))
            return(NULL) # rep group doesn't have feature groups

        # UNDONE: find some way to differentiate between actual errors?
        hom <- tryCatch(do.call(nontarget::homol.search, c(list(gt), homArgs)),
                        error = function(e) FALSE)

        if (is.logical(hom))
            return(NULL) # no results

        stopifnot(identical(hom[[3]][["HS IDs"]], hom[[3]][["HS cluster"]])) # what is the difference!?

        homTab <- as.data.table(hom[[3]])

        gNames <- groupTablesRG[[rg]]$group

        # peak IDs are now numeric indices of a subset of the original feature
        # groups. Replace them by group names to allow easy comparison.
        homTab[, groups := list(sapply(`peak IDs`, function(pids)
        {
            ginds <- as.integer(unlist(strsplit(pids, ",")))
            ginds <- ginds[order(gInfo[gNames[ginds], "mzs"])] # ensure they are from low to high m/z
            return(list(gNames[ginds]))
        }, simplify = FALSE))]

        homTab[, c("HS IDs", "peak IDs") := NULL]
        homTab[, (rg) := list(groups)]

        return(homTab)
    }, simplify = FALSE)
    homList <- homList[!sapply(homList, is.null)]

    compTab <- rbindlist(homList, fill = TRUE)

    if (nrow(compTab) == 0)
        return(components(componentInfo = data.table(), components = list(), algorithm = "nontarget"))

    # check if we can merge series
    compTab[, links := list(list())]

    for (r in seq_len(nrow(compTab)))
    {
        series <- compTab[["groups"]][[r]][[1]]
        nseries <- length(series)
        mzRange <- range(gInfo[series, "mzs"])
        links <- integer(0)
        for (ro in seq_len(nrow(compTab)))
        {
            if (r == ro)
                next

            if (abs(compTab[["m/z increment"]][r] - compTab[["m/z increment"]][ro]) > maxMzDev ||
                abs(compTab[["RT increment"]][r] - compTab[["RT increment"]][ro]) > maxRTDev)
                next # different series

            otherSeries <- compTab[["groups"]][[ro]][[1]]

            if (sum(series %in% otherSeries) < 1) # UNDONE: minimum overlap?
                next

            otherMzRange <- range(gInfo[otherSeries, "mzs"])

            # check for conflicts: groups that are not present in both and with
            # m/z values within the m/z range of the series
            missingG <- setdiff(series, otherSeries)
            missingOtherG <- setdiff(otherSeries, series)

            mzWithin <- function(mz1, mz2) numLTE(abs(mz1 - mz2), (compTab[["m/z increment"]][r] - maxMzDev))

            if (any(sapply(missingG, function(mg) any(mzWithin(gInfo[mg, "mzs"], gInfo[otherSeries, "mzs"])))) ||
                any(sapply(missingOtherG, function(mg) any(mzWithin(gInfo[mg, "mzs"], gInfo[series, "mzs"])))))
                next

            links <- c(links, ro)
        }

        set(compTab, r, "links", list(list(links)))
    }

    presentRGroups <- rGroups[sapply(rGroups, function(rg) !is.null(compTab[[rg]]))]

    # merge any pairs of series that are just pointing to each other: they are
    # either exactly the same or otherwise subsets of each other that can be
    # merged without conflicts

    compTab[, keep := TRUE]
    for (r in seq_len(nrow(compTab)))
    {
        if (compTab[["keep"]][r] == FALSE)
            next # already merged

        links <- compTab[["links"]][[r]]
        if (length(links) == 1)
        {
            otherLinks <- compTab[["links"]][[links[1]]]
            if (length(otherLinks) == 1 && otherLinks[1] == r)
            {
                # merge groups
                mGroups <- union(compTab[["groups"]][[r]][[1]], compTab[["groups"]][[links[1]]][[1]])
                mGroups <- mGroups[order(gInfo[mGroups, "mzs"])] # make sure order stays correct
                set(compTab, r, "groups", list(list(list(mGroups))))

                # mark presence
                l <- links[1] # BUG: cannot use "links" name in next line?
                lc <- compTab[l, presentRGroups, with = FALSE]
                for (rg in presentRGroups)
                {
                    if (!is.null(lc[[rg]][[1]]))
                        set(compTab, r, rg, list(list(list(lc[[rg]][[1]])))) # need to rewrap it in a list?
                }

                set(compTab, links[1], "keep", FALSE) # remove other
                set(compTab, r, "links", list(list(integer()))) # unlink
            }
        }
    }

    # clearout merged series
    compTab[, IDs := seq_len(.N)] # store current row order
    compTab <- compTab[keep == TRUE]
    # update links
    for (r in seq_len(nrow(compTab)))
    {
        links <- compTab[["links"]][[r]]
        if (length(links) > 0)
        {
            newLinks <- sapply(links, function(l) compTab[IDs == l, which = TRUE])
            set(compTab, r, "links", list(list(newLinks)))
        }
    }
    compTab[, c("keep", "IDs") := NULL]

    # split all rows in list with tables containing groups per row
    comps <- lapply(seq_len(nrow(compTab)), function(cmpi)
    {
        allGroups <- compTab[["groups"]][[cmpi]][[1]]
        homSeries <- seq_along(allGroups)
        ret <- rbindlist(lapply(presentRGroups, function(rg)
        {
            grp <- compTab[[rg]][[cmpi]][[1]]
            if (!is.null(grp))
                return(data.table(rt = gInfo[grp, "rts"], mz = gInfo[grp, "mzs"], group = grp,
                                  hsnr = match(grp, allGroups), rGroup = rg,
                                  intensity = groupTablesRG[[rg]][group %in% grp, get(rg)]))
            return(NULL)
        }), fill = TRUE)
        setorderv(ret, c("hsnr"))
    })
    names(comps) <- paste0("CMP", seq_along(comps))

    cInfo <- copy(compTab)
    cInfo[, c("groups", "HS cluster") := NULL]
    setnames(cInfo,
             c("m/z increment", "RT increment", "min. RT in series", "max. RT in series", "max.-min. RT"),
             c("mz_increment", "rt_increment", "rt_min", "rt_max", "rt_range"))
    cInfo[, name := names(comps)]
    cInfo[, size := sapply(comps, nrow)]

    # convert from fgroup lists to logical presence
    for (rg in presentRGroups)
        set(cInfo, j = rg, value = !sapply(cInfo[[rg]], is.null))

    ret <- components(componentInfo = cInfo, components = comps, algorithm = "nontarget")
    saveCacheData("componentsNontarget", ret, hash)

    return(ret)
}
