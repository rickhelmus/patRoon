#' @include utils.R
NULL

#' Adduct utilities
#'
#' Several utility functions to work with adducts.
#'
#' @param formula A \code{character} vector with formulae to convert.
#'
#' @templateVar plain TRUE
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
        if (length(fGroups) > 0 && nrow(annotations(fGroups)) == 0)
            stop("The input feature groups are not annotated and no adduct is specified. ",
                 "Please either perform annotations (eg with selectIons()) or specify an adduct.")
        return(NULL) # NULL adduct is fine
    }
    
    checkmate::assert(checkmate::checkString(adduct, min.chars = 1),
                      checkmate::checkClass(adduct, "adduct"),
                      .var.name = .var.name)
    return(as.adduct(adduct))
}

checkAndGetIonization <- function(ionization, fGroups, .var.name = "ionization", add = NULL)
{
    # if fGroups != NULL and the fGroups are annotated then ioniation arg is optional
    if (is.null(ionization))
    {
        if (length(fGroups) > 0 && nrow(annotations(fGroups)) == 0)
            stop("The input feature groups are not annotated and no ionization is specified. ",
                 "Please either perform annotations (eg with selectIons()) or specify an adduct.")
        return(getIonizationFromAnnTable(annotations(fGroups)))
    }
    
    checkmate::assertChoice(ionization, c("positive", "negative"), .var.name = .var.name, add = add)

    return(ionization)    
}

calculateMasses <- function(masses, adducts, type)
{
    stopifnot(type %in% c("neutral", "mz"))
    
    if (!is.list(adducts))
        adducts <- list(adducts)
    
    if (length(masses) > length(adducts))
        adducts <- rep(adducts, length.out = length(masses))
    else if (length(masses) < length(adducts))
        masses <- rep(masses, length.out = length(adducts))
    
    electronMass <- getFormulaMass("H", -1) - getFormulaMass("H", 0)
    
    return(mapply(masses, adducts, FUN = function(m, add)
    {
        # calculation is split by elemental M multiplier, addition/subtractions, and final charge (z)
        
        em <- if (length(add@add) > 0 && length(add@sub) > 0)
            sum(sapply(add@add, getFormulaMass)) - sum(sapply(add@sub, getFormulaMass))
        else if (length(add@add) > 0)
            sum(sapply(add@add, getFormulaMass))
        else if (length(add@sub) > 0)
            -(sum(sapply(add@sub, getFormulaMass)))
        else # [M]
            0
        
        # NOTE: invert charge (electron is lost in pos and gained in neg)
        em <- em + (electronMass * -add@charge)
        
        if (type == "mz")
            return(((m * add@molMult) + em) / abs(add@charge))
        
        # else neutralize
        return(((m * abs(add@charge)) - em) / add@molMult)
    }))
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

getFGroupAdducts <- function(gNames, annTable, adduct, format)
{
    if (!is.null(adduct))
        grpAdducts <- rep(list(adduct), length(gNames))
    else
        grpAdducts <- lapply(annTable$adduct, as.adduct)
    
    if (format != "generic")
        grpAdductsChr <- makeAlgoAdducts(grpAdducts, gNames, format)
    else
        grpAdductsChr <- sapply(grpAdducts, as.character)
    
    names(grpAdducts) <- names(grpAdductsChr) <- gNames
    
    return(list(grpAdducts = grpAdducts, grpAdductsChr = grpAdductsChr))
}

asAdductVec <- function(v, ...)
{
    uv <- unique(v)
    adds <- setNames(lapply(uv, as.adduct, ...), uv)
    return(adds[v])
}

normalizeAdducts <- function(adds, err = TRUE)
{
    unAdds <- unique(adds)
    nAddsObj <- setNames(lapply(unAdds, as.adduct, err = err), unAdds)
    nAdds <- sapply(nAddsObj, function(x) if (is.null(x)) NA_character_ else as.character(x))
    return(nAdds[adds])
}
    
# NOTE: this and below two functions are split to separate memoised functions to retain proper ref docs
doAsAdduct <- memoise(function(x, format, isPositive, charge, err = TRUE)
{
    if (is(x, "adduct"))
        return(x)
    if (is.na(x))
        return(x)
    
    # check first: this should fail immediately
    checkmate::assertChoice(format, c("generic", "sirius", "genform", "metfrag", "openms", "cliquems"))
    if (format == "openms")
        checkmate::assertInt(charge)
    
    ac <- checkmate::makeAssertCollection()
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
    checkmate::assertFlag(err, add = ac)
    checkmate::reportAssertions(ac)
    
    if (format == "generic" || format == "sirius" || format == "cliquems")
    {
        if (format == "generic" && x == "[M]") # special case
        {
            mult <- 1; charge <- 0; adds <- subs <- character()
        }
        else
        {
            x <- gsub(" ", "", x, fixed = TRUE) # e.g. RAMClustR includes an accidental adduct with space
            
            if (!grepl("^\\[[[:alnum:]\\+\\-]+\\][[:digit:]]*[\\+\\-]{1}$", x))
            {
                if (err)
                    stop("Wrong format! (forgot brackets or charge?)")
                return(NULL)
            }
            
            if (format == "cliquems")
                x <- sub("Cat", "M", x, fixed = TRUE)
            
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
            
            # CliqueMS sometimes adds amount of addition/subtraction elements before the element (eg 2H instead of H2)
            # --> swap
            if (format == "cliquems")
                x <- gsub("([\\+\\-])([[:digit:]]+)([[:alpha:]]+)", "\\1\\3\\2", x)
            
            adds <- sub("\\+", "", unlist(regmatches(x, gregexpr("[\\+]{1}[[:alnum:]]+", x))))
            subs <- sub("\\-", "", unlist(regmatches(x, gregexpr("[\\-]{1}[[:alnum:]]+", x))))
            if (endsWith(x, "-"))
                charge <- -charge
        }
    }
    else if (format == "genform")
    {
        gfadds <- GenFormAdducts()[adduct == x]
        if (nrow(gfadds) == 0)
        {
            if (err)
                stop("Invalid adduct for GenForm! See GenFormAdducts() for valid options.")
            return(NULL)
        }
        
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
        {
            if (err)
                stop("Invalid adduct for MetFrag! See MetFragAdducts() for valid options.")
            return(NULL)
        }
        
        mfadds <- mfadds[1] # in case there are multiple hits
        adds <- mfadds$add; subs <- mfadds$sub; charge <- mfadds$charge
    }
    else if (format == "openms")
    {
        fl <- splitFormulaToList(x)
        adds <- formulaListToString(fl[fl > 0])
        subs <- formulaListToString(abs(fl[fl < 0]))
        mult <- 1 # UNDONE: always one for OpenMS?
    }
    
    adds <- adds[nzchar(adds)]; subs <- subs[nzchar(subs)]
    
    return(adduct(add = adds, sub = subs, charge = charge, molMult = mult))
})

doCalculateIonFormula <- memoise(function(formula, adduct)
{
    checkmate::assertCharacter(formula)
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
})

doCalculateNeutralFormula <- memoise(function(formula, adduct)
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
})

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
#'   
#'   \code{"openms"} is the format used by the \command{MetaboliteAdductDecharger} tool.
#'   
#'   \code{"cliquems"} is the format used by \pkg{cliqueMS}.
#' @param isPositive A logical that specifies whether the adduct should be
#'   positive. Should only be set when \code{format="metfrag"} and \code{x} is a
#'   \code{numeric} identifier.
#' @param charge The final charge. Only needs to be set when \code{format="openms"}.
#' @param err If \code{TRUE} then an error will be thrown if conversion fails,
#'   otherwise returns without data.
#'   
#' @rdname adduct-utils
#' @export
as.adduct <- function(x, format = "generic", isPositive = NULL, charge = NULL,
                      err = TRUE) doAsAdduct(x, format, isPositive, charge, err = err)

#' @details \code{calculateIonFormula} Converts one or more neutral formulae to
#'   adduct ions.
#' @rdname adduct-utils
#' @export
calculateIonFormula <- function(formula, adduct) doCalculateIonFormula(formula, adduct)

#' @details \code{calculateNeutralFormula} Converts one or more adduct ions to
#'   neutral formulae.
#' @rdname adduct-utils
#' @export
calculateNeutralFormula <- function(formula, adduct) doCalculateNeutralFormula(formula, adduct)
