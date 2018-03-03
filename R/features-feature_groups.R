#' @include features.R

if (FALSE){
featuresFromFeatGroups <- setClass("featuresFromFeatGroups", contains = "features")

#' @export
convertFeatureGroupsToFeatures <- function(fgroups, fnames)
{
    flist <- list()
    for (fgi in seq_along(fgroups))
    {
        gt <- groups(fgroups[[fgi]])
        gi <- groupInfo(fgroups[[fgi]])

        # use group info as basis
        ft <- as.data.table(gi)
        setnames(ft, c("rts", "mzs"), c("ret", "mz"))
        ft[, ID := colnames(gt)]

        avgit <- transpose(gt[, lapply(.SD, mean)])
        setnames(avgit, 1, "intensity")
        ft <- cbind(ft, avgit)

        # dummy ranges
        ft[, retmin := ret - 3]
        ft[, retmax := ret + 3]
        ft[, mzmin := mz - 0.0025]
        ft[, mzmax := mz + 0.0025]

        flist[[fnames[[fgi]]]] <- ft
    }

    anaInfo <- data.frame(analysis = fnames, group = fnames, ref = "", path = ".", stringsAsFactors = F)

    return(featuresFromFeatGroups(features = flist, analysisInfo = anaInfo))
}
}
