# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
#' @include utils-param.R
NULL

#' @export
getComponentsSpecClustParamDefs <- paramConfigDefsFact(list(
    method = list(
        default = "complete",
        description = "Clustering method (passed to fastclust::hclust)",
        type = "string"
    ),
    specSimParams = list(
        default = getDefSpecSimParams(),
        description = "Default MS/MS spectral similarity parameters",
        type = "list",
        null.ok = FALSE
    ),
    maxTreeHeight = list(
        default = 1,
        description = "Maximum tree height (passed to dynamicTreeCut::cutreeDynamic)",
        type = "number"
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
        positive = TRUE
    ),
    IMS = list(
        default = "maybe",
        description = "Which feature groups are considered for componentization in IMS workflows ('maybe', 'both', FALSE, TRUE)",
        type = "IMS"
    )
))

#' @export
ComponentsSpecClustParam <- setClass("ComponentsSpecClustParam", contains = "param")

setMethod("initialize", "ComponentsSpecClustParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "ComponentsSpecClustParam", baseName = "ComponentsSpecClustParam",
                   description = "Parameters for MS/MS spectral clustering component generation", version = "1.0",
                   definitions = getComponentsSpecClustParamDefs(), ...)
})
