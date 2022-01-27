#' <%= if (exists("ion")) "@param ionization Which ionization polarity was used to generate the data: should be \\code{\"positive\"}
#'   or \\code{\"negative\"}. If the \\code{featureGroups} object has adduct annotations, and \\code{ionization=NULL}, the
#'   ionization will be detected automatically.
#'
#'   \\setsWF This parameter is not supported for sets workflows, as the ionization will always be detected
#'   automatically." %>
#' <%= if (exists("minSize")) paste("@param minSize The minimum size of a component. Smaller components than this size will be removed. See note below.", if (exists("RC")) "Sets the \\code{minModuleSize} argument to \\code{\\link[RAMClustR]{ramclustR}}." else "") %>
#' <%= if (exists("minReps")) "@param relMinReplicates Feature groups within a component are only kept when they contain data for at least this
#'   (relative) amount of replicate analyses. For instance, \\samp{0.5} means that at least half of the replicates should
#'   contain data for a particular feature group in a component. In this calculation replicates that are fully absent
#'   within a component are not taken in to account. See note below." %>
#' <%= if (exists("absMzDev")) paste("@param absMzDev Maximum absolute \\emph{m/z} deviation. Sets the", absMzDev) %>
#' <%= if (exists("method")) "@param method Clustering method that should be applied (passed to
#'   \\code{\\link[fastcluster:hclust]{fastcluster::hclust}})." %>
#' <%= if (!exists("noDots")) paste("@param \\dots", if (exists("myDots")) paste0(myDots, "\\cr\\cr") else "", "\\setsWF Further arguments passed to the non-sets workflow method.") %>
