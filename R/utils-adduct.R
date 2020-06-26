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

checkAndToAdduct <- function(adduct, na.ok = FALSE, .var.name = "adduct")
{
    checkmate::assert(checkmate::checkString(adduct, min.chars = 1, na.ok = na.ok),
                      checkmate::checkClass(adduct, "adduct"),
                      .var.name = .var.name)
    as.adduct(adduct)
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

#' @details \code{GenFormAdducts} returns a table with information on adducts
#'   supported by \command{GenForm}.
#' @rdname adduct-utils
#' @export
GenFormAdducts <- function()
{
    # adapated from GenForm/msaddion.cpp
    ret <- rbindlist(list(
        list("M+3H", "H3", "", "+3", 1),
        list("M+2H+Na", "H2Na", "", "+3", 1),
        list("M+H+2Na", "HNa2", "", "+3", 1),
        list("M+3Na", "Na3", "", "+3", 1),
        list("M+2H", "H2", "", "+2", 1),
        list("M+H+NH4", "H5N", "", "+2", 1),
        list("M+H+Na", " HNa", "", "+2", 1),
        list("M+H+K", "HK", "", "+2", 1),
        list("M+ACN+2H", "C2H5N", "", "+2", 1),# ACN=C2H3N
        list("M+2Na", "Na2", "", "+2", 1),
        list("M+2ACN+2H", "C4H8N2", "", "+2", 1),
        list("M+3ACN+2H", "C6H11N3", "", "+2", 1),
        list("M+H", "H", "", "+1", 1),
        list("M+NH4", "H4N", "", "+1", 1),
        list("M+Na", "Na", "", "+1", 1),
        list("M+CH3OH+H", "CH5O", "", "+1", 1),
        list("M+K", "K", "", "+1", 1),
        list("M+ACN+H", "C2H4N", "", "+1", 1),
        list("M+2Na-H", "Na2", "H", "+1", 1),
        list("M+IsoProp+H", "C3H10O", "", "+1", 1), # IsoProp=C3H9O(?)
        list("M+ACN+Na", "C2H3NNa", "", "+1", 1),
        list("M+2K-H", "K2", "H", "+1", 1),
        list("M+DMSO+H", "C2H7OS", "", "+1", 1), # DMSO=C2H6OS
        list("M+2ACN+H", "C4H7N2", "", "+1", 1),
        list("M+IsoProp+Na+H", "C3H10ONa", "", "+1", 1),
        list("2M+H", "H", "", "+1", 2),
        list("2M+NH4", "H4N", "", "+1", 2),
        list("2M+Na", "Na", "", "+1", 2),
        list("2M+K", "K", "", "+1", 2),
        list("2M+ACN+H", "C2H4N", "", "+1", 2),
        list("2M+ACN+Na", "C2H3NNa", "", "+1", 2),

        list("M-3H", "", "H3", "-3", 1),
        list("M-2H", "", "H2", "-2", 1),
        list("M-H2O-H", "", "H3O", "-1", 1),
        list("M-H", "", "H", "-1", 1),
        list("M+Na-2H", "Na", "H2", "-1", 1),
        list("M+Cl", "Cl", "", "-1", 1),
        list("M+K-2H", "K", "H2", "-1", 1),
        list("M+FA-H", "CH2O2", "H", "-1", 1), # FA=CH2O2
        list("M+Hac-H", "C2H4O2", "H", "-1", 1), # HAc=C2H4O2
        list("M+Br", "Br", "", "-1", 1),
        list("M+TFA-H", "C2HF3O2", "H", "-1", 1), # TFA=C2HF3O2
        list("2M-H", "", "H", "-1", 2),
        list("2M+FA-H", "", "H", "-1", 2),
        list("2M+Hac-H", "C2H4O2", "H", "-1", 2),
        list("3M-H", "", "H", "-1", 3),

        list("M-e", "", "", "+1", 1),
        list("M+e", "", "H", "-1", 1),

        list("-e", "", "", "+1", 1),
        list("+e", "", "", "-1", 1),
        list("+H", "H", "", "+1", 1),
        list("-H", "", "H", "-1", 1),
        list("+Na", "Na", "", "+1", 1)
    ))

    setnames(ret, c("adduct", "add", "sub", "charge", "molMult"))
    ret[, charge := as.numeric(charge)]
    return(ret[])
}

#' @details \code{MetFragAdducts} returns a table with information on adducts
#'   supported by \command{MetFrag}.
#' @rdname adduct-utils
#' @export
MetFragAdducts <- function()
{
    # from http://ipb-halle.github.io/MetFrag/projects/metfragcl/ and Constants.java in MetFragLib
    ret <- rbindlist(list(
        list(1, "[M+H]+", "H", "", 1),
        list(18, "[M+NH4]+", "NH4", "", 1),
        list(23, "[M+Na]+", "Na", "", 1),
        list(39, "[M+K]+", "K", "", 1),
        list(33, "[M+CH3OH+H]+", "CH5O", "", 1), # MeOH+H
        list(42, "[M+ACN+H]+", "C2H4N", "", 1), # ACN+H
        list(64, "[M+ACN+Na]+", "C2H3NNa", "", 1), # ACN+Na
        list(83, "[M+2ACN+H]+", "C4H7N2", "", 1), # 2ACN+H
        list(0, "[M]+", "", "", 1),

        list(-1, "[M-H]-", "", "H", -1),
        list(35, "[M+Cl]-", "Cl", "", -1),
        list(45, "[M+HCOO]-", "CO2H", "", -1),
        list(59, "[M+CH3COO]-", "C2H3O2", "", -1),
        list(0, "[M]-", "", "", -1)
    ))
    setnames(ret, c("adduct_mode", "adduct_type", "add", "sub", "charge"))
    return(ret)
}

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
    else if (format == "genform")
    {
        gfadds <- GenFormAdducts()[adduct == x]
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
                mfadds <- MetFragAdducts()[adduct_mode == x & charge == cha]
            }
            else
                mfadds <- MetFragAdducts()[adduct_mode == x]
        }
        else
            mfadds <- MetFragAdducts()[adduct_type == x]

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
        return(f)
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
        return(f)
    }, USE.NAMES = FALSE)
}
