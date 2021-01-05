#' Adduct utilities
#'
#' Several utility functions to work with adducts.
#'
#' @param formula A \code{character} vector with formulae to convert.
#'
#' @template adduct-arg
#'
#' @examples as.adduct("[M+H]+")
#' as.adduct("[M+H2]2+")
#' as.adduct("[2M+H]+")
#' as.adduct("[M-H]-")
#' as.adduct("+H", format = "genform")
#' as.adduct(1, isPositive = TRUE, format = "metfrag") # MetFrag adduct ID 1 --> returns [M+H]+
#'
#' calculateIonFormula("C2H4O", "[M+H]+") # C2H5O
#' calculateNeutralFormula("C2H5O", "[M+H]+") # C2H4O
#'
#' @name adduct-utils
NULL

checkAndToAdduct <- function(adduct, fGroups = NULL, .var.name = "adduct")
{
    # if fGroups != NULL and the fGroups are annotated then adducts may be optional
    if (!is.null(fGroups) && is.null(adduct))
    {
        if (nrow(annotations(fGroups)) == 0)
            stop("The input feature groups are not annotated and no adduct is specified. ",
                 "Please either perform annotations (eg with mergeIons()) or specify an adduct.")
        return(NULL) # NULL adduct is fine
    }
    
    checkmate::assert(checkmate::checkString(adduct, min.chars = 1),
                      checkmate::checkClass(adduct, "adduct"),
                      .var.name = .var.name)
    return(as.adduct(adduct))
}

adductMZDelta <- function(adduct)
{
    # NOTE: rcdk::get.formula makes elemental counts absolute
    ret <- 0
    getMZ <- function(form, charge) rcdk::get.formula(form, charge)@mass
    if (length(adduct@add) > 0 && length(adduct@sub) > 0) # NOTE: add charge once
        ret <- sum(sapply(adduct@add, getMZ, charge = adduct@charge)) - sum(sapply(adduct@sub, getMZ, charge = 0))
    else if (length(adduct@add) > 0)
        ret <- sum(sapply(adduct@add, getMZ, charge = adduct@charge))
    else if (length(adduct@sub) > 0)
        ret <- -(sum(sapply(adduct@sub, getMZ, charge = adduct@charge)))
    else # [M]
        ret <- getMZ("H", adduct@charge) - getMZ("H", 0) # electron mass
    
    return(ret)
}

makeAlgoAdducts <- function(adducts, gNames, format)
{
    cat("Converting to algorithm specific adducts... ")
    
    ret <- (mapply(adducts, gNames, FUN = function(add, grp)
    {
        addChr <- as.character(add, format = format, err = FALSE)
        if (is.na(addChr))
        {
            addChr <- if (add@charge < 0) "M-H" else "M+H"
            warning(sprintf("Changing the adduct for group %s from %s to %s because it is not supported by %s",
                            grp, as.character(add), if (add@charge < 0) "[M-H]-" else "[M+H]+", format))
        }
        return(addChr)
    }))
    
    cat("Done!\n")
    return(ret)
}

#' @details \code{GenFormAdducts} returns a table with information on adducts
#'   supported by \command{GenForm}.
#' @rdname adduct-utils
#' @export
GenFormAdducts <- function() copy(adductsGF)

#' @details \code{MetFragAdducts} returns a table with information on adducts
#'   supported by \command{MetFrag}.
#' @rdname adduct-utils
#' @export
MetFragAdducts <- function() copy(adductsMF)

#' @details \code{as.adduct} Converts an object in to an \code{\link{adduct}}
#'   object.
#' @param x The object that should be converted. Should be a \code{character}
#'   string, a \code{numeric} \command{MetFrag} adduct identifier
#'   (\code{adduct_mode} column obtained with \code{MetFragAdducts}) or an
#'   \code{\link{adduct}} object (in which case no conversion occurs).
#' @param format A \code{character} that specifies the source format.
#'
#'   \code{"generic"} is an internally used generic format that supports full
#'   textual conversion. Examples: \code{"[M+H]+"}, \code{"[2M+H]+"},
#'   \code{"[M+3H]3+"}.
#'
#'   \code{"sirius"} Is the format used by \command{SIRIUS}. It is similar to
#'   \code{generic} but does not allow multiple charges/molecules. See the
#'   SIRIUS manual for more details.
#'
#'   \code{"genform"} and \code{"metfrag"} support fixed types of adducts
#'   which can be obtained with the \code{GenFormAdducts} and
#'   \code{MetFragAdducts} functions, respectively.
#' @param isPositive A logical that specifies whether the adduct should be
#'   positive. Should only be set when \code{format="metfrag"} and \code{x} is a
#'   \code{numeric} identifier.
#' @rdname adduct-utils
#' @export
as.adduct <- function(x, format = "generic", isPositive = NULL)
{
    if (is(x, "adduct"))
        return(x)
    if (is.na(x))
        return(x)

    ac <- checkmate::makeAssertCollection()
    checkmate::assertChoice(format, c("generic", "sirius", "genform", "metfrag")) # don't add: this should fail immediately
    if (format == "metfrag")
    {
        checkmate::assert(checkmate::checkInt(x),
                          checkmate::checkString(x, min.chars = 1),
                          .var.name = "x")
        if (is.numeric(x) && x == 0)
            checkmate::assertFlag(isPositive, add = ac)
    }
    else
        checkmate::assertString(x, min.chars = 1, add = ac)
    checkmate::reportAssertions(ac)

    if (format == "generic" || format == "sirius")
    {
        if (format == "generic" && x == "[M]") # special case
        {
            mult <- 1; charge <- 0; adds <- subs <- character()
        }
        else
        {
            if (!grepl("^\\[.+\\].*[\\+\\-]{1}", x))
                stop("Wrong format! (forgot brackets or charge?)")
            
            if (format == "sirius")
                mult <- charge <- 1
            else
            {
                mult <- as.numeric(sub("^\\[([0-9]*)M.*", "\\1", x))
                if (is.na(mult))
                    mult <- 1
                charge <- as.numeric(sub(".*\\]([0-9]*)[\\+\\-]{1}$", "\\1", x))
                if (is.na(charge))
                    charge <- 1
            }
            
            adds <- sub("\\+", "", unlist(regmatches(x, gregexpr("[\\+]{1}[[:alnum:]]+", x))))
            subs <- sub("\\-", "", unlist(regmatches(x, gregexpr("[\\-]{1}[[:alnum:]]+", x))))
            if (endsWith(x, "-"))
                charge <- -charge
        }
    }
    else if (format == "genform")
    {
        gfadds <- getAdductTable("genform")[adduct == x]
        if (nrow(gfadds) == 0)
            stop("Invalid adduct for GenForm! See GenFormAdducts() for valid options.")

        gfadds <- gfadds[1] # in case there are multiple hits
        adds <- gfadds$add; subs <- gfadds$sub; charge <- gfadds$charge; mult <- gfadds$molMult
    }
    else if (format == "metfrag")
    {
        mult <- 1
        if (is.numeric(x))
        {
            if (x == 0)
            {
                # molecular ion is twice in the list (both polarities)
                cha <- if (isPositive) 1 else -1
                mfadds <- getAdductTable("metfrag")[adduct_mode == x & charge == cha]
            }
            else
                mfadds <- getAdductTable("metfrag")[adduct_mode == x]
        }
        else
            mfadds <- getAdductTable("metfrag")[adduct_type == x]

        if (nrow(mfadds) == 0)
            stop("Invalid adduct for MetFrag! See MetFragAdducts() for valid options.")

        mfadds <- mfadds[1] # in case there are multiple hits
        adds <- mfadds$add; subs <- mfadds$sub; charge <- mfadds$charge
    }

    adds <- adds[nzchar(adds)]; subs <- subs[nzchar(subs)]

    return(adduct(add = adds, sub = subs, charge = charge, molMult = mult))
}

#' @details \code{calculateIonFormula} Converts one or more neutral formulae to
#'   adduct ions.
#' @rdname adduct-utils
#' @export
calculateIonFormula <- function(formula, adduct)
{
    checkmate::assertCharacter(formula, min.chars = 1)
    adduct <- checkAndToAdduct(adduct)

    sapply(formula, function(f)
    {
        if (adduct@molMult > 1)
            f <- Reduce(addFormula, rep(f, adduct@molMult-1), f)
        if (length(adduct@add) > 0)
            f <- Reduce(addFormula, adduct@add, f)
        if (length(adduct@sub) > 0)
            f <- Reduce(subtractFormula, adduct@sub, f)
        return(simplifyFormula(f))
    }, USE.NAMES = FALSE)
}

#' @details \code{calculateNeutralFormula} Converts one or more adduct ions to
#'   neutral formulae.
#' @rdname adduct-utils
#' @export
calculateNeutralFormula <- function(formula, adduct)
{
    checkmate::assertCharacter(formula, min.chars = 1)
    adduct <- checkAndToAdduct(adduct)

    sapply(formula, function(f)
    {
        if (length(adduct@add) > 0)
            f <- Reduce(subtractFormula, adduct@add, f)
        if (length(adduct@sub) > 0)
            f <- Reduce(addFormula, adduct@sub, f)
        if (adduct@molMult > 1)
        {
            fl <- splitFormulaToList(f)
            fl <- fl / adduct@molMult
            f <- formulaListToString(fl)
        }
        return(simplifyFormula(f))
    }, USE.NAMES = FALSE)
}
