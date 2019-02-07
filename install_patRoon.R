yesNo <- function(...) menu(c("Yes", "No"), ...) == 1
packageInstalled <- function(pkg) requireNamespace(pkg, quietly = TRUE)

checkPackages <- function(pkgs, ask = TRUE, bioc = FALSE)
{
    notInstalled <- pkgs[!sapply(pkgs, packageInstalled)]
    if (length(notInstalled) == 0)
        return(NULL)
    
    if (ask && !yesNo(title = paste0("Install the following required packages?\n",
                                     paste0(notInstalled, collapse = ", "))))
        stop("Aborted. Please install the package(s) manually.")

    if (bioc)
        BiocManager::install(notInstalled, ask = FALSE)
    else
        install.packages(notInstalled)
}

needOptionalPackage <- function(pkg, msg)
{
    return(!packageInstalled(pkg) &&
               yesNo(title = sprintf("The optional package '%s' is not installed. %s\nInstall?", pkg, msg)))
}

addToRProfile <- function(setOpts)
{
    alreadySet <- sapply(names(setOpts), function(o) nzchar(getOption(o, default = "")), USE.NAMES = FALSE)
    
    if (any(alreadySet))
    {
        if (!yesNo(title = sprintf("NOTE: the following options are already set: %s\nDo you still want to add them to your ~/.Rprofile?",
                                   names(setOpts)[alreadySet])))
            setOpts <- setOpts[!alreadySet]
    }
    
    if (length(setOpts) > 0)
    {
        cat(c(paste0("\n\n# Automatically added by install_patRoon R script on ", date()),
              sprintf("options(%s = \"%s\")", names(setOpts), setOpts)),
            file = "~/.Rprofile", append = TRUE, sep = "\n")
        
        # also set in local session
        do.call(options, as.list(setOpts))
    }
}

installPackages <- function()
{
    checkPackages(c("installr", "BiocManager"))
    checkPackages(c("mzR", "xcms", "CAMERA"), bioc = TRUE)
    
    if (needOptionalPackage("RDCOMClient", "This is only required for interfacing with Bruker DataAnalysis."))
        checkPackages("RDCOMClient", ask = FALSE)
    
    if (needOptionalPackage("RAMClustR", paste("This package may be used for componentization (e.g. grouping adducts/isotopes).",
                                               "To install this package R tools is required and will be installed automatically if not yet installed.")))
    {
        checkPackages(c("remotes", "pkgbuild"), ask = FALSE)
        if (!pkgbuild::check_rtools())
            installr::install.rtools(check = FALSE, GUI = FALSE)
        remotes::install_github("cbroeckl/RAMClustR", build_vignettes = TRUE, dependencies = TRUE)
    }

    # UNDONE: check JAVA
}

installExtDeps <- function(extPath)
{
    if (Sys.info()[["machine"]] != "x86-64")
        warning("This script probbaly only works well on a 64 bit system!")
    
    extPath <- normalizePath(extPath, mustWork = FALSE)
    if (!dir.exists(extPath))
    {
        if (!dir.create(extPath))
            stop(sprintf("Failed to create external dependency directory '%s'", extPath))
    }
    
    hasOpenMS <- system2("FeatureFinderMetabo", "--help", stdout = FALSE, stderr = FALSE) == 0
    
    mfBin <- path.expand(getOption("patRoon.path.metFragCL"))
    hasMF <- !is.null(mfBin) && nzchar(mfBin) && file.exists(mfBin)
    
    
    
    instChoices <- c(if (hasOpenMS) "OpenMS (seems already installed)" else "(not installed)",
                     "MetFrag CL", "SIRIUS", "pngquant")
    instWhat <- select.list(instChoices, preselect = if (!hasOpenMS) instChoices[1],
                            title = "Which external dependencies should be installed?",
                            graphics = FALSE, multiple = TRUE)
    
    if (length(instWhat) > 0)
    {
        if (any(grepl("OpenMS", instWhat)))
            installr::install.URL("https://github.com/OpenMS/OpenMS/releases/download/Release2.4.0/OpenMS-2.4.0-Win64.exe", message = FALSE)
        
        setOpts <- character()
        
        if ("MetFrag CL" %in% instWhat)
        {
            url <- "http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-CL.jar"
            dest <- file.path(extPath, basename(url))
            if (download.file(url, dest) != 0)
                warning(paste("Failed to download MetFrag CL from ", url))
            else
                setOpts <- c(setOpts, patRoon.path.metFragCL = dest)
        }
        
        if ("SIRIUS" %in% instWhat)
        {
            url <- "https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.0.1/sirius-4.0.1-win64-headless.zip"
            dest <- file.path(extPath, basename(url))
            if (download.file(url, dest) != 0)
                warning(paste("Failed to download SIRIUS from ", url))
            else
            {
                unzip(dest, exdir = extPath)
                zipdest <- file.path(extPath, "sirius-win64-headless-4.0.1")
                if (!file.exists(zipdest))
                    warning(paste("Failed to extract SIRIUS to ", extPath))
                else
                    setOpts <- c(setOpts, patRoon.path.SIRIUS = zipdest)
                unlink(dest)
            }
        }
        
        if ("pngquant" %in% instWhat)
        {
            url <- "https://pngquant.org/pngquant-windows.zip"
            dest <- file.path(extPath, basename(url))
            if (download.file(url, dest) != 0)
                warning(paste("Failed to download pngquant from ", url))
            else
            {
                unzip(dest, exdir = extPath)
                zipdest <- file.path(extPath, "pngquant")
                if (!file.exists(zipdest))
                    warning(paste("Failed to extract pngquant to ", extPath))
                else
                    setOpts <- c(setOpts, patRoon.path.pngquant = zipdest)
                unlink(dest)
            }
        }
        
        if (length(setOpts) > 0 && yesNo(title = "Do you want to add the location of the downloaded tools to ~/.Rprofile (so you don't need to specify them manually with options())?"))
            addToRProfile(setOpts)
    }
}

installPatRoonPackages <- function(exampleData)
{
    remotes::install_github("rickhelmus/patRoon")
    if (exampleData)
        remotes::install_github("rickhelmus/patRoonData")
    invisible(NULL)
}

installPatRoon <- function(what = c("packages", "external_deps", "patRoon"),
                           extDepPath = "~/patRoon-dependencies", exampleData = TRUE)
{
    validWhat <- c("packages", "external_deps", "patRoon")
    if (!is.character(what) || any(!what %in% validWhat))
        stop(sprintf("what must be a subset of (%s)", paste0(validWhat, collapse = ", ")))
    
    if (!is.character(extDepPath) || length(extDepPath) > 1)
        stop("extDepPath must be valid character string.")

    if ("packages" %in% what)
        installPackages()
    if ("external_deps" %in% what)
        installExtDeps(extDepPath)
    if ("patRoon" %in% what)
        installPatRoonPackages(exampleData)
    
    cat("All done! You may need to restart R.")
    invisible(NULL)
}
