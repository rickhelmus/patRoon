# TODO
# - convert to wizard style?
# - shiny/tcltk GUI?
# - optionally don't unload packages? optionally clear symbols?


# put utility functions in separate environment
utils <- setRefClass("utilsInst", methods = list(

    # force forward slash so we don't have to escape text written to files
    fixPath = function(p) normalizePath(p, winslash = "/", mustWork = FALSE),

    yesNo = function(title, ...) menu(c("Yes", "No"), title = title, ...) == 1,
    
    quoteVariables = function(x) ifelse(sapply(x, function(q) length(q) == 1 && is.character(q)),
                                         paste0("\"", x, "\""), as.character(x)),
    
    packageInstalled = function(pkg)
    {
        # HACK: rJava will be checked afterwards, just see if it's installed
        if (pkg == "rJava")
            return(nzchar(system.file(package = pkg)))
        requireNamespace(pkg, quietly = TRUE)
    },
    
    ensureInstPath = function(instPath)
    {
        if (!dir.exists(instPath))
        {
            if (!dir.create(instPath, recursive = TRUE))
                stop(sprintf("Failed to create external dependency directory '%s'", instPath))
        }
    },
    
    unloadAllPackages = function()
    {
        # adapted from https://stackoverflow.com/a/39235076
        
        loadedPkgs <- names(sessionInfo()$otherPkgs)
        if (is.null(loadedPkgs))
            return(NULL) # nothing to do...
        
        for (pkg in loadedPkgs)
            tryCatch(detach(paste0("package:", pkg), character.only = TRUE, unload = TRUE, force = TRUE),
                     error = function(e) FALSE) # UNDONE: do something with error?
    },
    
    printHeader = function(txt)
    {
        hd <- strrep("-", 25)
        cat(hd, txt, hd, sep = "\n")
    },
    
    checkPackages = function(pkgs, lib, ask = TRUE, type = "installp", ghRepos = NULL, force = FALSE, ...)
    {
        if (force)
            notInstalled <- pkgs
        else
            notInstalled <- pkgs[!sapply(pkgs, utils$packageInstalled)]
        
        if (length(notInstalled) == 0)
            return(NULL)
        
        if (ask && !yesNo(title = paste0("Install the following required packages?\n",
                                         paste0(notInstalled, collapse = ", "))))
            stop("Aborted. Please install the package(s) manually.")
        
        if (type == "bioc")
        {
            cmd <- "BiocManager::install"
            args <- list(notInstalled, ask = TRUE)
        }
        else if (type == "gh")
        {
            cmd <- "remotes::install_github"
            args <- c(list(paste(ghRepos, notInstalled, sep = "/")), list(upgrade_dependencies = FALSE, force = force))
        }
        else
        {
            cmd <- "install.packages"
            args <- list(notInstalled)
        }
        
        args <- c(args, list(...))
        if (!is.null(lib))
            args <- c(args, list(lib = lib))
        
        argNames <- names(args)
        argVals <- quoteVariables(args)
        argsTxt <- paste0(ifelse(nzchar(argNames), paste0(argNames, " = "), ""), argVals, collapse = ", ")
        cat(sprintf("Installing packages: %s\n\nEXECUTE: %s(%s)\n\n", paste0(notInstalled, collapse = ", "), cmd, argsTxt))

        # hopefully to avoid stupid DLL problems etc
        unloadAllPackages()
        gc()
        
        # eval(parse()): make sure we can call namespace prefixed functions
        ret <- tryCatch(do.call(eval(parse(text = cmd)), args), error = function(e) paste("Error:", e))
        
        # update
        notInstalled <- notInstalled[!sapply(notInstalled, utils$packageInstalled)]
        
        # something not installed? (note that package install functions may give errors which are not always that important)
        if (length(notInstalled) > 0)
        {
            if (is.character(ret))
                cat(ret, "\n")
            stop(paste("There were errors/warnings during package installation. Please review the output in the console.",
                       sprintf("Alternatively, try to install the following packages manually with '%s':\n%s", cmd,
                               paste0(notInstalled, collapse = ", ")),
                       sep = "\n"))
        }
    },
    
    needOptionalPackage = function(pkg, msg)
    {
        return(!packageInstalled(pkg) &&
                   yesNo(title = sprintf("The optional package '%s' is not installed. %s\nInstall?", pkg, msg)))
    },
    
    # From utils.R
    getCommandWithOptPath = function(cmd, opt)
    {
        if (Sys.info()[["sysname"]] == "Windows")
            cmd <- paste0(cmd, ".exe") # add file extension for Windows
        
        opt <- paste0("patRoon.path.", opt)
        path <- getOption(opt)
        if (!is.null(path) && nzchar(path))
        {
            ret <- file.path(path.expand(path), cmd)
            if (!file.exists(ret))
                return("")
            return(ret)
        }
        
        # assume command is in PATH --> no need to add path
        if (!nzchar(Sys.which(cmd)))
            return("")
        
        return(cmd)
    },
    
    addToRProfile = function(setOpts, JavaPath)
    {
        # UNDONE: need this?
        if (FALSE && length(setOpts) > 0)
        {
            alreadySet <- sapply(names(setOpts), function(o) nzchar(getOption(o, default = "")), USE.NAMES = FALSE)
            if (any(alreadySet))
            {
                if (!yesNo(title = sprintf("NOTE: the following options are already set: %s\nDo you still want to add them to your ~/.Rprofile?",
                                           paste0(names(setOpts)[alreadySet], collapse = ", "))))
                    setOpts <- setOpts[!alreadySet]
            }
        }
        
        if (length(setOpts) > 0 || nzchar(JavaPath))
        {
            RProfPath <- "~/.Rprofile"
            ownRProfPath <- "~/.Rprofile-patRoon.R"
            
            # use a marker in options to see if everything was loaded        
            optMarker <- "patRoon.Rprof"
            setOpts <- c(setOpts, setNames(list(TRUE), optMarker))
            options(setNames(list(FALSE), optMarker))
            
            RProfFile <- c(paste("# Automatically generated by the install_patRoon R script on", date()),
                           "# This file should not be modified and should be loaded by ~/.Rprofile.",
                           sprintf("options(%s = %s)", names(setOpts), quoteVariables(setOpts)))
            if (nzchar(JavaPath))
                RProfFile <- c(RProfFile, sprintf("Sys.setenv(PATH = paste(Sys.getenv('PATH'), '%s', sep = ';'))",
                                                  fixPath(file.path(JavaPath, "bin"))))

            if (file.exists(ownRProfPath))
            {
                # file already exists, make sure to retain options that are already set
                
                cat("Found existing options file. Merging settings.\n")
                
                f <- readLines(ownRProfPath)
                
                # discard things that we want to write
                f <- f[!grepl("^#", f)] # remove comments
                if (length(setOpts) > 0)
                    f <- f[!grepl(paste0(names(setOpts), collapse = "|"), f)]
                if (nzchar(JavaPath))
                    f <- f[!grepl("Sys.setenv", f)]
                
                if (length(f) > 0)
                    RProfFile <- c(RProfFile, f) # add the rest of the existing contents
            }
            
            cat(RProfFile, file = ownRProfPath, sep = "\n")
            
            # source it: make sure it works and set options in this environment
            sret <- tryCatch(suppressWarnings(source(ownRProfPath, local = TRUE)), error = function(e) FALSE)
            if (is.logical(sret) || !getOption(optMarker, FALSE))
                stop("Failed to generate proper Rprofile settings")

            options(setNames(list(NULL), optMarker)) # reset
            
            cat("Checking your ~/.Rprofile... ")
            
            if (file.exists(RProfPath))
            {
                # there is an ~/.Rprofile: see if it already loads our file (because the installer was already executed before)

                sret <- tryCatch(suppressWarnings(source(RProfPath, local = TRUE)), error = function(e) FALSE)
                if (is.logical(sret))
                    stop("Failed to load current Rprofile!")
            }
            
            # either no current .Rprofile or it doesn't load our settings yet
            if (!getOption(optMarker, FALSE))
            {
                cat(c(paste("\n\n# Automatically added by install_patRoon script on ", date()),
                      sprintf("if (file.exists('%s'))", ownRProfPath),
                      sprintf("    source('%s')", ownRProfPath)),
                    file = RProfPath, sep = "\n", append = TRUE)
                cat("Updated!\n")
            }
            else
                cat("Seems OK, skipping updates\n")
            
            options(setNames(list(NULL), optMarker))
        }
    },
    
    # returns TRUE if Java path needs to be added to Rprofile
    installPrerequisites = function(instPath, lib)
    {
        printHeader("Installing prerequisites...")
        
        ret <- FALSE
        
        checkPackages(c("installr", "BiocManager", "rJava", "remotes", "pkgbuild"), lib)
        
        # see if Java can be executed (e.g. necessary for MetFrag) and can be loaded (e.g. needed for RCDK)
        hasJava <- suppressWarnings(system2("java", "-version", stdout = FALSE, stderr = FALSE) == 0) &&
            tryCatch(rJava::.jinit(), error = function(e) -1) == 0
        if (!hasJava)
        {
            if (yesNo(paste("Could not detect a suitable Java JDK. Do you want automatically install it?",
                            "(The JDK will only be accessible from R and not interfere with the rest of your system)",
                            sep = "\n")))
            {
                ensureInstPath(instPath)
                ret <- installr::install.jdk(path = instPath)
            }
        }
        
        if (!suppressMessages(pkgbuild::has_rtools()))
        {
            if (!yesNo("Rtools doesn't seem to be installed. This is necessary to proceed the installation. Do you want to install Rtools now?"))
                stop("Please install Rtools manually and re-run the installer.")
            # NOTE: set keep_install_file to avoid long delays after installation
            installr::install.rtools(check = FALSE, GUI = FALSE, keep_install_file = TRUE)
        }
        
        return(ret)
    },
    
    # returns which options should be set in Rprofile
    installExtDeps = function(instPath)
    {
        printHeader("Installing external dependencies...")
        
        setOpts <- list()
        
        extDeps <- data.frame(name = c("ProteoWizard", "OpenMS", "SIRIUS", "pngquant"),
                              command = c("msconvert", "FeatureFinderMetabo", "sirius-64", "pngquant"),
                              copt = c("pwiz", "OpenMS", "SIRIUS", "pngquant"),
                              stringsAsFactors = FALSE)
        extDeps$path <- mapply(extDeps$command, extDeps$copt, FUN = utils$getCommandWithOptPath)
        extDeps <- rbind(extDeps, list(name = "MetFrag CL", command = "", copt = "",
                                       path = getOption("patRoon.path.metFragCL", "")))
        
        present <- nzchar(extDeps$path) & file.exists(extDeps$path)
        instChoices <- paste(extDeps$name, ifelse(present, "(seems installed)", "(doesn't seem to be installed)"))
        choices <- c(instChoices, "Missing", "All", "None")
        instWhat <- select.list(choices, instChoices[!present], TRUE, graphics = FALSE,
                                title = "Which external dependencies should be installed?")
        
        if ("All" %in% instWhat)
            instWhat <- instChoices
        else if ("Missing" %in% instWhat)
            instWhat <- instChoices[!present]
        
        if ("None" %in% instWhat)
            instWhat <- character()
        
        # convert back to simple names
        instWhat <- extDeps$name[choices %in% instWhat]
        
        if (length(instWhat) > 0)
        {
            ensureInstPath(instPath)
            
            if ("OpenMS" %in% instWhat)
            {
                # NOTE: set keep_install_file to avoid long delays after installation
                installr::install.URL("https://github.com/OpenMS/OpenMS/releases/download/Release2.4.0/OpenMS-2.4.0-Win64.exe",
                                      message = FALSE, keep_install_file = TRUE)
            }
            
            if ("MetFrag CL" %in% instWhat)
            {
                url <- "http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.4.5-CL.jar"
                dest <- fixPath(file.path(instPath, basename(url)))
                if (download.file(url, dest) != 0)
                    warning(paste("Failed to download MetFrag CL from ", url))
                else
                    setOpts <- c(setOpts, list(patRoon.path.metFragCL = dest))
            }
            
            if ("SIRIUS" %in% instWhat)
            {
                url <- "https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.0.1/sirius-4.0.1-win64-headless.zip"
                dest <- file.path(instPath, basename(url))
                if (download.file(url, dest) != 0)
                    warning(paste("Failed to download SIRIUS from ", url))
                else
                {
                    unzip(dest, exdir = instPath)
                    zipdest <- fixPath(file.path(instPath, "sirius-win64-headless-4.0.1"))
                    if (!file.exists(zipdest))
                        warning(paste("Failed to extract SIRIUS to ", instPath))
                    else
                        setOpts <- c(setOpts, list(patRoon.path.SIRIUS = zipdest))
                    unlink(dest)
                }
            }
            
            if ("pngquant" %in% instWhat)
            {
                url <- "https://pngquant.org/pngquant-windows.zip"
                dest <- file.path(instPath, basename(url))
                if (download.file(url, dest) != 0)
                    warning(paste("Failed to download pngquant from ", url))
                else
                {
                    unzip(dest, exdir = instPath)
                    zipdest <- fixPath(file.path(instPath, "pngquant"))
                    if (!file.exists(zipdest))
                        warning(paste("Failed to extract pngquant to ", instPath))
                    else
                        setOpts <- c(setOpts, list(patRoon.path.pngquant = zipdest))
                    unlink(dest)
                }
            }
            
            if ("ProteoWizard" %in% instWhat &&
                yesNo(title = paste("Due to license agreement restrictions ProteoWizard cannot be installed automatically at this point.",
                                    "Do you want to open the webpage so that you can download and install ProteoWizard manually?",
                                    sep = "\n")))
            {
                browseURL("http://proteowizard.sourceforge.net/download.html")
                while(!yesNo(title = "Did you install ProteoWizard and are ready to continue the patRoon installation?")) {}
            }
        }
        
        return(setOpts)
    },
    
    installRDeps = function(lib)
    {
        printHeader("Pre-Installing R dependencies...")
        
        checkPackages(c("mzR", "xcms", "CAMERA"), lib, type = "bioc")
        
        if (needOptionalPackage("RDCOMClient", "This is only required for interfacing with Bruker DataAnalysis."))
            checkPackages("RDCOMClient", lib, ask = FALSE, repos = "http://www.omegahat.net/R")
        
        if (needOptionalPackage("RAMClustR", "This package may be used for componentization (e.g. grouping adducts/isotopes)."))
            checkPackages("RAMClustR", lib, ask = FALSE, type = "gh", ghRepos = "cbroeckl", build_vignettes = TRUE, dependencies = TRUE)
    },
    
    installPatRoonPackages = function(exampleData, lib, force)
    {
        printHeader("Installing patRoon R package(s)...")
        
        checkPackages("patRoon", lib, ask = FALSE, type = "gh", ghRepos = "rickhelmus", force = force)
        if (exampleData)
            checkPackages("patRoonData", lib, ask = FALSE, type = "gh", ghRepos = "rickhelmus", force = force)
        invisible(NULL)
    }
))()


installPatRoon <- function(what = c("prereq", "external_deps", "packages", "patRoon"),
                           instPath = "~/patRoon-install", exampleData = TRUE,
                           lib = NULL, force = FALSE)
{
    if (Sys.info()[["sysname"]] != "Windows" || Sys.info()[["machine"]] != "x86-64")
        stop("Sorry, this script only works on a 64 bit Windows system at the moment.")
    
    validWhat <- c("prereq", "external_deps", "packages", "patRoon")
    if (!is.character(what) || any(!what %in% validWhat))
        stop(sprintf("what must be a subset of (%s)", paste0(validWhat, collapse = ", ")))
    
    if (!is.character(instPath) || length(instPath) > 1)
        stop("instPath must be valid character string.")
    
    utils$printHeader(paste("This installation script will install the patRoon R package.",
                            "The code uses some functions and other code from the installr package",
                            "(https://github.com/talgalili/installr)",
                            "Credits go to its authors.",
                            sep = "\n"))
    
    cat("Installation configuration:",
        paste("  - what:", paste0(what, collapse = ", ")),
        paste("  - Installation path for external dependencies:", instPath),
        paste("  - Install patRoon example data:", exampleData),
        paste("  - R library location:", if (is.null(lib)) "default" else paste0(lib, collapse = ", ")),
        sep = "\n")
    
    cat("\nNOTE: It is best to run this in a _fresh_ R session!\n")
    
    instPath <- utils$fixPath(instPath)
    
    didJava <- FALSE; setOpts <- list()
    if ("prereq" %in% what)
        didJava <- utils$installPrerequisites(instPath, lib)
    if ("external_deps" %in% what)
        setOpts <- utils$installExtDeps(instPath)
    
    if ((didJava || length(setOpts) > 0) &&
        utils$yesNo(paste("The installer can add code to your ~/.Rprofile file to automatically configure the location of downloaded tools and/or Java.",
                          "An additional file will be created (~/Rprofile-patRoon.R) that will set the necessary options and is sourced from your ~/.Rprofile",
                          "If you do not do this you will have to set the location of downloaded tools (e.g. MetFrag, SIRIUS) and/or Java manually (not recommended)",
                          "Continue?",
                          sep = "\n")))
    {
        jPath <- Sys.getenv("JAVA_HOME") # should be set by installr
        utils$addToRProfile(setOpts, if (didJava) jPath else "")
    }
        
    if ("packages" %in% what)
        utils$installRDeps(lib)
    if ("patRoon" %in% what)
        utils$installPatRoonPackages(exampleData, lib, force)
    
    utils$printHeader("All done! You may need to restart R.")
    
    invisible(NULL)
}

cat("Run installPatRoon() to start the installation\n")
