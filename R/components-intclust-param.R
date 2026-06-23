# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

#' @export
getComponentsIntClustParamDefs <- paramConfigDefsFact(list(
    method = list(
        default = "complete",
        description = "Clustering method (passed to fastclust::hclust)",
        type = "string"
    ),
    metric = list(
        default = "euclidean",
        description = "Distance metric used for clustering (passed to daisy)",
        type = "string"
    ),
    normalized = list(
        default = TRUE,
        description = "Whether intensities are normalized before clustering",
        type = "flag"
    ),
    average = list(
        default = TRUE,
        description = "Whether to average intensities across replicates before clustering",
        type = "flag"
    ),
    maxTreeHeight = list(
        default = 1,
        description = "Maximum tree height (passed to dynamicTreeCut::cutreeDynamic)",
        type = "number",
        typeCheckArgs = list(finite = TRUE)
    ),
    deepSplit = list(
        default = TRUE,
        description = "Enable deep splitting (passed to dynamicTreeCut::cutreeDynamic)",
        type = "flag"
    ),
    minModuleSize = list(
        default = 1,
        description = "Minimum module (component) size (passed to dynamicTreeCut::cutreeDynamic)",
        type = "count",
        typeCheckArgs = list(positive = TRUE)
    ),
    IMS = list(
        default = "maybe",
        description = "Which feature groups are considered for componentization in IMS workflows ('maybe', 'both', FALSE, TRUE)",
        type = "IMS"
    )
))

#' @export
ComponentsIntClustParam <- setClass("ComponentsIntClustParam", contains = "param")

setMethod("initialize", "ComponentsIntClustParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsIntClustParam", baseName = "ComponentsIntClustParam",
                   description = "Parameters for intensity-clustering component generation", version = "1.0",
                   definitions = getComponentsIntClustParamDefs(), ...)
})
