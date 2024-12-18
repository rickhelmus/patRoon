# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

defaultPkgOpts <- function(pkgname)
{
    ret <- list(cache.mode = "both",
                checkCentroided = TRUE,
                cache.fileName = "cache.sqlite",
                MP.method = "classic",
                # backwards compat
                MP.maxProcs = getOption("patRoon.maxProcAmount", parallel::detectCores(logical = FALSE)),
                MP.futureSched = 1.0,
                MP.logPath = "log",
                path.pwiz = "",
                path.GenForm = "",
                path.MetFragCL = getOption("patRoon.path.metFragCL", ""), # backwards compat
                path.MetFragCompTox = "",
                path.MetFragPubChemLite = "",
                path.SIRIUS = "",
                path.OpenMS = "",
                path.pngquant = "",
                path.obabel = "",
                path.BioTransformer = "")
    return(setNames(ret, paste0(pkgname, ".", names(ret))))
}

dumpPkgOpts <- function(printFunc)
{
    for (opt in names(defaultPkgOpts(utils::packageName())))
        printFunc(sprintf("- %s: \"%s\"", opt, getOption(opt)))
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
