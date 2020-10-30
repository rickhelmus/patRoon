defaultPkgOpts <- function(pkgname)
{
    ret <- list(cache.mode = "both",
                cache.fileName = "cache.sqlite",
                maxProcAmount = parallel::detectCores(logical = FALSE),
                logPath = "log",
                path.pwiz = "",
                path.GenForm = "",
                path.MetFragCL = getOption("patRoon.path.metFragCL", ""), # backwards compat
                path.MetFragCompTox = "",
                path.MetFragPubChemLite = "",
                path.SIRIUS = "",
                path.OpenMS = "",
                path.pngquant = "",
                path.obabel = "")
    return(setNames(ret, paste0(pkgname, ".", names(ret))))
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
    for (opt in names(defaultPkgOpts(pkgname)))
        packageStartupMessage(sprintf("- %s: \"%s\"", opt, getOption(opt)))
}
