#' @include utils-adduct.R
NULL

# Notes on general format:
# SIRIUS: [M+ADDUCT]+/-; [M-ADDUCT]+/-; [M]+/- (only single charge)
# GenForm: predefined list: see msaddion.cpp. Ex: +H, -H, -e, M+3H, 2M+H, M+2K-H
# MetFrag: adduct = addition/subtraction mass, isPositive = T/F --> see http://ipb-halle.github.io/MetFrag/projects/metfragcl/


#' Generic adduct class
#'
#' Objects from this class are used to specify adduct information in an
#' algorithm independent way.
#'
#' @slot add,sub A \code{character} with one or more formulas to add/subtract.
#' @slot molMult How many times the original molecule is present in this
#'   molecule (\emph{e.g.} for a dimer this would be \samp{2}). Default is \samp{1}.
#' @slot charge The final charge of the adduct (default \samp{1}).
#'
#' @usage adduct(...)
#'
#' @param \dots Any of \code{add}, \code{sub}, \code{molMult} and/or \code{charge}. See \verb{Slots}.
#' @param x,object An \code{adduct} object.
#'
#' @examples adduct("H") # [M+H]+
#' adduct(sub = "H", charge = -1) # [M-H]-
#' adduct(add = "K", sub = "H2", charge = -1) # [M+K-H2]+
#' adduct(add = "H3", charge = 3) # [M+H3]3+
#' adduct(add = "H", molMult = 2) # [2M+H]+
#'
#' as.character(adduct("H")) # returns "[M+H]+"
#'
#' @seealso \code{\link{as.adduct}} for easy creation of \code{adduct} objects
#'   and \link[=adduct-utils]{adduct utilities} for other adduct functionality.
#'
#' @export adduct
#' @exportClass adduct
adduct <- setClass("adduct", slots = c(add = "character", sub = "character",
                                       molMult = "numeric", charge = "numeric"))

setMethod("initialize", "adduct", function(.Object, add = character(), sub = character(), molMult = 1, charge = 1)
{
    .Object <- callNextMethod(.Object, add = add, sub = sub, charge = charge, molMult = molMult)

    checkmate::assertCount(.Object@molMult, positive = TRUE)
    checkmate::assertInt(.Object@charge)

    if (length(.Object@add) > 0)
        .Object@add <- sapply(.Object@add, simplifyFormula, USE.NAMES = FALSE)
    if (length(.Object@sub) > 0)
        .Object@sub <- sapply(.Object@sub, simplifyFormula, USE.NAMES = FALSE)

    return(.Object)
})

#' @describeIn adduct Shows summary information for this object.
#' @export
setMethod("show", "adduct", function(object)
{
    printf("An adduct object ('%s')\n", class(object))

    printf("Addition: %s\n", paste0(object@add, collapse = ", "))
    printf("Subtraction: %s\n", paste0(object@sub, collapse = ", "))
    printf("Molecule multiplier: %d\n", object@molMult)
    printf("Charge: %s%d\n", if (object@charge > 0) "+" else "", object@charge)

    printf("Generic textual representation: %s\n", as.character(object))

    showObjectSize(object)
})

#' @describeIn adduct Converts an \code{adduct} object to a specified
#'   \code{character} format.
#' @inheritParams as.adduct
#' @param err If \code{TRUE} then an error will be thrown if conversion fails,
#'   otherwise returns \code{NA}.
#' @export
setMethod("as.character", "adduct", function(x, format = "generic", err = TRUE)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(format, c("generic", "sirius", "genform", "metfrag"), add = ac)
    checkmate::assertFlag(err, add = ac)
    checkmate::reportAssertions(ac)

    getGFMFAdduct <- function()
    {
        simplifyAddSub <- function(x) if (length(x) == 0) "" else simplifyFormula(paste0(x, collapse = ""))
        
        tbl <- if (format == "genform") GenFormAdducts() else MetFragAdducts()
        colName <- if (format == "genform") "adduct" else "adduct_type"
        adds <- tbl[add == simplifyAddSub(x@add) & sub == simplifyAddSub(x@sub) &
                        charge == x@charge & molMult == x@molMult, colName, with = FALSE]

        if (length(adds) == 0)
        {
            if (err)
                stop(sprintf("Invalid adduct for %s! See %s() for valid options.",
                             format, if (format == "GenForm") "GenFormAdducts" else "MetFragAdducts"))
            return(NA_character_)
        }
        
        return(adds[1]) # return first one in case there are multiple hits
    }
    
    if (format == "sirius" || format == "generic")
    {
        if (format == "sirius")
        {
            if (abs(x@charge) != 1)
                stop("SIRIUS only supports a charge of +/- 1")
            if (x@molMult > 1)
                stop("SIRIUS only supports a molecular multiplier of 1")
        }

        adds <- if (length(x@add)) paste0("+", x@add, collapse = "") else ""
        subs <- if (length(x@sub)) paste0("-", x@sub, collapse = "") else ""

        charge <- if (x@charge > 0) "+" else if (x@charge < 0) "-"
        if (abs(x@charge) > 1)
            charge <- paste0(abs(x@charge), charge)
        mult <- if (x@molMult > 1) x@molMult else ""

        return(paste0("[", mult, "M", adds, subs, "]", charge))
    }
    else if (format == "genform")
        return(getGFMFAdduct())
    else if (format == "metfrag")
    {
        if (abs(x@charge) != 1)
        {
            if (err)
                stop("MetFrag only supports a charge of +/- 1")
            return(NA_character_)
        }
        
        if (x@molMult > 1)
        {
            if (err)
                stop("MetFrag only supports a molecular multiplier of 1")
            return(character())
        }

        return(getGFMFAdduct())
    }
})

