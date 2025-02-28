# SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

defaultPkgOpts <- function(pkgname)
{
    ret <- list(cache.mode = "both",
                checkCentroided = TRUE,
                cache.fileName = "cache.sqlite",
                MS.backends = getMSReadBackends(),
                MS.preferIMS = FALSE,
                threads = parallel::detectCores(logical = FALSE),
                MP.method = "classic",
                # backwards compat
                MP.maxProcs = getOption("patRoon.maxProcAmount", parallel::detectCores(logical = FALSE)),
                MP.futureSched = 1.0,
                MP.logPath = "log",
                path.BrukerTIMS = "",
                path.pwiz = "",
                path.GenForm = "",
                path.MetFragCL = getOption("patRoon.path.metFragCL", ""), # backwards compat
                path.MetFragCompTox = "",
                path.MetFragPubChemLite = "",
                path.SIRIUS = "",
                path.OpenMS = "",
                path.obabel = "",
                path.BioTransformer = "",
                path.limits = "")
    return(setNames(ret, paste0(pkgname, ".", names(ret))))
}

dumpPkgOpts <- function(printFunc)
{
    for (opt in names(defaultPkgOpts(utils::packageName())))
    {
        msg <- sprintf("- %s: %s", opt, paste0('"', getOption(opt), '"', collapse = ", "))
        if (opt == "patRoon.MS.backends")
            msg <- paste(msg, sprintf("(available: %s)", paste0('"', availableBackends(verbose = FALSE), '"', collapse = ", ")))
        printFunc(msg)
    }
}

.onLoad <- function(libname, pkgname)
{
    # intialize any options that are unset
    dopts <- defaultPkgOpts(pkgname)
    missingOpts <- !names(dopts) %in% names(options())
    if (length(missingOpts) > 0)
        options(dopts[missingOpts])
}

.onAttach <- function(libname, pkgname)
{
    packageStartupMessage(sprintf("Welcome to %s %s!", pkgname, utils::packageVersion(pkgname)))
    packageStartupMessage("Configuration:")
    dumpPkgOpts(packageStartupMessage)
}
