# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

#' @include param.R
NULL

#' @export
getFilterComponentsParamDefs <- paramConfigDefsFact(list(
    size = list(
        advanced = FALSE,
        default = NULL,
        description = "Component size range",
        type = "range",
        typeCheckArgs = list(null.ok = TRUE)
    ),
    adducts = list(
        advanced = FALSE,
        default = NULL,
        description = "Adducts to keep"
    ),
    isotopes = list(
        advanced = FALSE,
        default = NULL,
        description = "Isotope assignments to keep"
    ),
    rtIncrement = list(
        advanced = FALSE,
        default = NULL,
        description = "Retention increment range for homologous series",
        type = "numeric",
        typeCheckArgs = list(len = 2, any.missing = FALSE, null.ok = TRUE)
    ),
    mzIncrement = list(
        advanced = FALSE,
        default = NULL,
        description = "m/z increment range for homologous series",
        type = "numeric",
        typeCheckArgs = list(len = 2, lower = 0, any.missing = FALSE, null.ok = TRUE)
    ),
    checkComponentsSession = list(
        advanced = FALSE,
        default = NULL,
        description = "Component check session file"
    ),
    sets = list(
        advanced = TRUE,
        default = NULL,
        description = "Sets to keep",
        type = "sets",
        typeCheckArgs = list(multiple = TRUE)
    ),
    negate = list(
        advanced = TRUE,
        default = FALSE,
        description = "Apply filters in the opposite manner",
        type = "flag"
    ),
    verbose = list(
        advanced = TRUE,
        default = TRUE,
        description = "Verbose output",
        type = "flag"
    )
))

FilterComponentsParam <- setClass("FilterComponentsParam", contains = "param")
setMethod("initialize", "FilterComponentsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterComponentsParam", baseName = "filter_components",
                   subDir = "components", description = "Parameters for component filtering", version = "1.0",
                   definitions = getFilterComponentsParamDefs(), ...)
})

setValidity("FilterComponentsParam", function(object)
{
    parsFilled <- paramListFillDefaults(object@data, object@definitions)
    ac <- checkmate::makeAssertCollection()
    checkmate::assert(
        checkmate::checkNull(parsFilled$adducts),
        checkmate::checkFlag(parsFilled$adducts),
        checkmate::checkList(parsFilled$adducts, types = "adduct"),
        checkmate::checkCharacter(parsFilled$adducts, min.chars = 1),
        .var.name = "adducts", add = ac
    )
    checkmate::assert(
        checkmate::checkFlag(parsFilled$isotopes, null.ok = TRUE),
        checkmate::checkIntegerish(parsFilled$isotopes, lower = 0),
        .var.name = "isotopes", add = ac
    )
    if (!is.logical(parsFilled$checkComponentsSession))
        assertCheckSession(parsFilled$checkComponentsSession, mustExist = TRUE, null.ok = TRUE,
                           .var.name = "checkComponentsSession", add = ac)
    tryCatch(checkmate::reportAssertions(ac), error = function(e) e)
})


getFilterComponentsTPsParamDefs <- paramConfigDefsFact(
    list(
        retDirMatch = list(
        advanced = TRUE,
            default = FALSE,
            description = "Require matching retention direction",
            type = "flag"
        ),
        minSpecSim = list(
        advanced = TRUE,
            default = NULL,
            description = "Minimum spectral similarity",
            type = "number",
            typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
        ),
        minSpecSimPrec = list(
        advanced = TRUE,
            default = NULL,
            description = "Minimum precursor-shifted spectral similarity",
            type = "number",
            typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
        ),
        minSpecSimBoth = list(
        advanced = TRUE,
            default = NULL,
            description = "Minimum both-shifted spectral similarity",
            type = "number",
            typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
        ),
        minTotFragMatches = list(
        advanced = TRUE,
            default = NULL,
            description = "Minimum total fragment matches",
            type = "number",
            typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
        ),
        minTotNLMatches = list(
        advanced = TRUE,
            default = NULL,
            description = "Minimum total neutral loss matches",
            type = "number",
            typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
        ),
        minFragMatches = list(
        advanced = TRUE,
            default = NULL,
            description = "Minimum fragment matches",
            type = "number",
            typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
        ),
        minNLMatches = list(
        advanced = TRUE,
            default = NULL,
            description = "Minimum neutral loss matches",
            type = "number",
            typeCheckArgs = list(lower = 0, finite = TRUE, null.ok = TRUE)
        ),
        negate = list(
        advanced = TRUE,
            default = FALSE,
            description = "Apply filters in the opposite manner",
            type = "flag"
        ),
        verbose = list(
        advanced = TRUE,
            default = TRUE,
            description = "Verbose output",
            type = "flag"
        )
    )
)

FilterComponentsTPsParam <- setClass("FilterComponentsTPsParam", contains = "param")
setMethod("initialize", "FilterComponentsTPsParam", function(.Object, ...)
{
    callNextMethod(.Object, name = "FilterComponentsTPsParam", baseName = "filter_components_tps",
                   subDir = "components", description = "Parameters for TP component filtering", version = "1.0",
                   definitions = getFilterComponentsTPsParamDefs(), ...)
})

