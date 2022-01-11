# TODO
# - convert to wizard style?
# - shiny/tcltk GUI?
# - optionally don't unload packages? optionally clear symbols?
# - backup ~/.Rprofile
# - allow patRoon GitHub installation of regardless of repos?


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
    
    packagesNotInstalled = function(pkgs) pkgs[!sapply(pkgs, utils$packageInstalled)],
    
    instPathRLib = function(instPath) file.path(instPath, "library"),
    
    ensureInstPath = function(instPath, rlib = FALSE)
    {
        if (rlib)
            instPath <- instPathRLib(instPath)
        if (!dir.exists(instPath))
        {
            if (!dir.create(instPath, recursive = TRUE))
                stop(sprintf("Failed to create external dependency directory '%s'", instPath))
        }
    },
    
    getLibPaths = function(instPath, pkgWhere)
    {
        lpaths <- .libPaths()
        if (pkgWhere == "pDepsIso")
        {
            utils$ensureInstPath(instPath, TRUE)
            lpaths <- c(utils$instPathRLib(instPath), .libPaths())
        }
        return(lpaths)
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
    
    checkPackages = function(pkgs, pkgWhere, ask = TRUE, type = "installp", repos = NULL, force = FALSE,
                             pkgReposNames = pkgs, ...)
    {
        if (force)
            notInstalled <- pkgs
        else
            notInstalled <- packagesNotInstalled(pkgs)
        
        if (length(notInstalled) == 0)
            return(NULL)
        
        if (ask && !yesNo(title = paste0("Install the following required packages?\n",
                                         paste0(notInstalled, collapse = ", "))))
            stop("Aborted. Please install the package(s) manually.")
        
        if (pkgWhere == "normal")
        {
            if (type == "bioc")
            {
                cmd <- "BiocManager::install"
                args <- list(notInstalled, update = FALSE)
            }
            else if (type == "gh")
            {
                # HACK: in some cases the package repos name differs from the installed package name (ie KPIC<-->KPIC2)
                pkgToInstall <- pkgReposNames[notInstalled %in% pkgs]
                cmd <- "remotes::install_github"
                args <- c(list(paste(repos, pkgToInstall, sep = "/")), list(upgrade = FALSE, force = force))
            }
            else
            {
                cmd <- "install.packages"
                args <- list(notInstalled)
                if (!is.null(repos))
                    args <- c(args, list(repos = repos))
            }
        }
        else
        {
            cmd <- "install.packages"
            args <- c(list(repos = "https://rickhelmus.github.io/patRoonDeps/", type = "binary"),
                      list(notInstalled))
            
            # shouldn't be necessary as .libPaths was already set
            # if (pkgWhere == "pDepsIso")
            #     args <- c(args, list(lib = instPathRLib(instPath)))
        }
        args <- c(args, list(...))
        
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
    
    makeOwnRProfile = function(instPath, setOpts, JavaPath, isolatedPackages)
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
        
        ownRProfPath <- fixPath("~/.Rprofile-patRoon.R")
        
        # use a marker in options to see if everything was loaded        
        optMarker <- "patRoon.Rprof"
        setOpts <- c(setOpts, setNames(list(TRUE), optMarker))
        options(setNames(list(FALSE), optMarker))
        
        RProfHeader <- c(paste("# Automatically generated by the install_patRoon R script on", date()),
                         "# This file should not be modified and should be loaded by ~/.Rprofile.")
        RProfSettings <- sprintf("options(%s = %s)", names(setOpts), quoteVariables(setOpts))
        if (!is.null(JavaPath))
        {
            RProfSettings <- c(RProfSettings, sprintf("Sys.setenv(PATH = paste(Sys.getenv('PATH'), '%s', sep = ';'))",
                                                      fixPath(file.path(JavaPath, "bin"))))
            RProfSettings <- c(RProfSettings, sprintf("Sys.setenv(JAVA_HOME = '%s')", fixPath(JavaPath)))
        }
        if (isolatedPackages) # UNDONE: add before or after existing lib paths?
            RProfSettings <- c(RProfSettings, sprintf(".libPaths(c(.libPaths(), \"%s\"))",
                                                      fixPath(instPathRLib(instPath))))
        
        RProfExists <- character()
        if (file.exists(ownRProfPath))
        {
            # file already exists, make sure to retain options that are already set
            
            cat("Found existing options file. Merging settings.\n")
            
            f <- readLines(ownRProfPath)
            
            # discard things that we want to write
            f <- f[!grepl("^#", f)] # remove comments
            if (length(setOpts) > 0)
                f <- f[!grepl(gsub(".", "\\.", paste0(names(setOpts), collapse = "|"), fixed = TRUE), f)]
            if (isolatedPackages)
                f <- f[!grepl("libPaths", f)]
            if (!is.null(JavaPath))
                f <- f[!grepl("^Sys\\.setenv", f)]
            
            if (length(f) > 0)
                RProfExists <- f
        }
        
        cat(c(RProfHeader, RProfExists, RProfSettings), file = ownRProfPath, sep = "\n")
        
        # source it: make sure it works and set options in this environment
        sret <- tryCatch(suppressWarnings(source(ownRProfPath, local = TRUE)), error = function(e) FALSE)
        if (is.logical(sret) || !getOption(optMarker, FALSE))
            stop("Failed to generate proper Rprofile settings")
        
        options(setNames(list(NULL), optMarker)) # reset
        
        return(ownRProfPath)
    },

    editRProfile = function(ownRProfPath)
    {
        RProfPath <- "~/.Rprofile"
        
        # use a marker in options to see if everything was loaded        
        optMarker <- "patRoon.Rprof"
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
            cat("Seems already done, skipping updates\n")
        
        options(setNames(list(NULL), optMarker))
    },
    
    findPWizPath = function()
    {
        # try to find ProteoWizard
        # order: options --> win registry --> PATH
        # the PATH is searched last because OpenMS might have added its own old version.
        
        path <- getOption("patRoon.path.pwiz")
        if (!is.null(path) && nzchar(path))
            return(path)
        
        # Inspired by scan_registry_for_rtools() from pkgload
        key <- "Directory\\shell\\Open with SeeMS\\command"
        reg <- tryCatch(utils::readRegistry(key, "HCR"), error = function(e) NULL)
        
        # not sure if this might occur
        if (is.null(reg))
            reg <- tryCatch(utils::readRegistry(key, "HLM"), error = function(e) NULL)
        
        if (!is.null(reg))
        {
            path <- tryCatch(dirname(sub("\"([^\"]*)\".*", "\\1", reg[[1]])), error = function(e) NULL)
            if (!is.null(path) && file.exists(file.path(path, "msconvert.exe"))) # extra check: see if msconvert is there
                return(path)
        }
        
        # check PATH
        path <- dirname(Sys.which("msconvert.exe"))
        if (nzchar(path))
            return(path)
        
        return(NULL)
    },
    
    downloadFile = function(instPath, what, url, doUnzip, destF = basename(url))
    {
        dest <- file.path(instPath, destF)
        if (download.file(url, dest, mode = "wb") != 0)
        {
            warning(sprintf("Failed to download %s from '%s'", what, url))
            return(NULL)
        }
        else if (doUnzip)
        {
            zipdest <- fixPath(file.path(instPath, tools::file_path_sans_ext(basename(url))))
            unzip(dest, exdir = zipdest)
            if (!file.exists(zipdest))
            {
                warning(paste("Failed to extract %s to '%s'", what, instPath))
                return(NULL)
            }
            unlink(dest)
            return(zipdest)
        }
        return(dest)
    },
    
    installMandatoryRDeps = function(instPath, pkgWhere)
    {
        printHeader("Pre-Installing mandatory R dependencies...")
        
        curLPaths <- .libPaths(); on.exit(.libPaths(curLPaths))
        .libPaths(getLibPaths(instPath, pkgWhere))
        
        packagesCRAN <- packagesNotInstalled(c("installr", "BiocManager", "rJava", "remotes", "pkgbuild"))
        packagesBioC <- packagesNotInstalled(c("mzR", "xcms", "CAMERA"))
        mandatoryPackages <- c(packagesCRAN, packagesBioC)
        
        choices <- character()
        if (length(mandatoryPackages) == 0)
        {
            cat("All mandatory R pre-dependencies already installed\n")
            return(NULL)
        }

        cat("The following mandatory packages will be installed:",
            paste0(mandatoryPackages, collapse = ", "),
            sep = "\n")
        
        checkPackages(packagesCRAN, pkgWhere, ask = FALSE)
        checkPackages(packagesBioC, pkgWhere, ask = FALSE, type = "bioc")
    },
    
    installOptionalRDeps = function(instPath, pkgWhere)
    {
        printHeader("Pre-Installing optional R dependencies...")
        
        curLPaths <- .libPaths(); on.exit(.libPaths(curLPaths))
        .libPaths(getLibPaths(instPath, pkgWhere))
        
        optPackages <- packagesNotInstalled(c("KPIC", "RAMClustR", "cliqueMS", "MetaClean", "RDCOMClient"))
        
        choices <- character()
        if ("KPIC" %in% optPackages)
            choices <- c(choices, KPIC2 = "KPIC2 (algorithm to find/group features)")
        if ("RAMClustR" %in% optPackages)
            choices <- c(choices, RAMClustR = "RAMClustR (algorithm for componentization)")
        if ("cliqueMS" %in% optPackages)
            choices <- c(choices, cliqueMS = "cliqueMS (algorithm for componentization)")
        if ("MetaClean" %in% optPackages)
            choices <- c(choices, MetaClean = "MetaClean (to calculate chromatographic peak qualities)")
        if ("RDCOMClient" %in% optPackages)
            choices <- c(choices, RDCOMClient = "RDCOMClient (required for Bruker DataAnalysis integration)")
        if (length(choices) > 1)
            choices <- c(choices, all = "All")
        else if (length(choices) == 0)
        {
            cat("All optional R pre-dependencies already installed\n")
            return(NULL)
        }

        instWhat <- select.list(choices, multiple = TRUE, graphics = FALSE,
                                title = "Which of the following R packages do you want to install?")
        instWhat <- names(instWhat)

        if (!is.null(instWhat))
        {
            if (any(c("all", "KPIC2") %in% instWhat))
            {
                checkPackages("ropls", pkgWhere, ask = FALSE, type = "bioc")
                checkPackages("KPIC", pkgWhere, ask = FALSE, type = "gh", repos = "rickhelmus", pkgReposNames = "KPIC2")
            }
            if (any(c("all", "RAMClustR") %in% instWhat))
                checkPackages("RAMClustR", pkgWhere, ask = FALSE, type = "gh", repos = "cbroeckl")
            if (any(c("all", "cliqueMS") %in% instWhat))
                checkPackages("cliqueMS", pkgWhere, ask = FALSE, type = "gh", repos = "rickhelmus")
            if (any(c("all", "MetaClean") %in% instWhat))
            {
                checkPackages(c("BiocStyle", "Rgraphviz"), pkgWhere, ask = FALSE, type = "bioc")
                checkPackages("MetaClean", pkgWhere, ask = FALSE, type = "gh", repos = "KelseyChetnik")
            }
            if (any(c("all", "RDCOMClient") %in% instWhat))
                checkPackages("RDCOMClient", pkgWhere, ask = FALSE, type = "gh", repos = "BSchamberger")
        }
    },
    
    # returns JDK path if it was downloaded
    installExtDeps = function(instPath, checkRtools)
    {
        printHeader("Installing external dependencies...")
        
        ret <- NULL
        
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
                url <- "https://download.java.net/java/GA/jdk13/5b8a42f3905b406298b72d750b6919f6/33/GPL/openjdk-13_windows-x64_bin.zip"
                down <- downloadFile(instPath, "Open JDK 13", url, TRUE)
                if (!is.null(down))
                    ret <- paste0(down, "/jdk-13")
                
            }
            else
                cat("NOTE: Please make sure to install a proper JDK before loading patRoon.\n")
        }
        
        if (checkRtools && !suppressMessages(pkgbuild::has_rtools()))
        {
            if (!yesNo("Rtools doesn't seem to be installed. This is necessary to proceed the installation. Do you want to install Rtools now?"))
                stop("Please install Rtools manually and re-run the installer.")
            # NOTE: set keep_install_file to avoid long delays after installation
            installr::install.rtools(check = FALSE, GUI = FALSE, keep_install_file = TRUE)
        }
        
        return(ret)
    },
    
    # returns which options should be set in Rprofile
    installTools = function(instPath)
    {
        printHeader("Installing external tools...")
        
        setOpts <- list()
        
        extDeps <- data.frame(name = c("ProteoWizard", "OpenMS", "SIRIUS", "OpenBabel", "pngquant"),
                              command = c("msconvert", "FeatureFinderMetabo", "sirius", "obabel", "pngquant"),
                              copt = c("pwiz", "OpenMS", "SIRIUS", "obabel", "pngquant"),
                              stringsAsFactors = FALSE)
        extDeps$path <- mapply(extDeps$command, extDeps$copt, FUN = utils$getCommandWithOptPath)
        
        pwiz <- findPWizPath()
        # NOTE: ProteoWizard/msConvert may be installed via OpenMS, which is not the version you typically want to use.
        # Since OpenMS adds their version to PATH, we use findPWizPath() to find the right version
        if (!is.null(pwiz))
            extDeps$path[1] <- pwiz
        
        extDeps <- rbind(extDeps, list(name = "MetFrag CL", command = "", copt = "",
                                       path = getOption("patRoon.path.MetFragCL", "")))
        extDeps <- rbind(extDeps, list(name = "MetFrag CompTox WasteWater DB", command = "", copt = "",
                                       path = getOption("patRoon.path.MetFragCompTox", "")))
        extDeps <- rbind(extDeps, list(name = "MetFrag PubChemLite DB", command = "", copt = "",
                                       path = getOption("patRoon.path.MetFragPubChemLite", "")))
        extDeps <- rbind(extDeps, list(name = "BioTransformer", command = "", copt = "",
                                       path = getOption("patRoon.path.BioTransformer", "")))
        
        # NOTE: only check presence for non-commands, since these should be in PATH if $path isn't empty
        present <- nzchar(extDeps$path) & (file.exists(extDeps$path) | nzchar(extDeps$command))
        instChoices <- paste(extDeps$name, ifelse(present, "(seems installed)", "(doesn't seem to be installed)"))
        choices <- c(instChoices, "Missing", "All", "None")
        instWhat <- select.list(choices, multiple = TRUE, graphics = FALSE,
                                title = "Which external tools should be installed?")
        
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
                installr::install.URL("https://github.com/OpenMS/OpenMS/releases/download/Release2.7.0/OpenMS-2.7.0-Win64.exe",
                                      message = FALSE, keep_install_file = TRUE)
            }
            
            if ("MetFrag CL" %in% instWhat)
            {
                down <- downloadFile(instPath, "MetFrag CL", "https://github.com/rickhelmus/patRoonDeps/raw/master/ext/MetFrag2.4.5-CL.jar", FALSE)
                if (!is.null(down))
                    setOpts <- c(setOpts, list(patRoon.path.MetFragCL = down))
            }
            
            if ("MetFrag CompTox WasteWater DB" %in% instWhat)
            {
                down <- downloadFile(instPath, "MetFrag CompTox database", "https://zenodo.org/record/3472781/files/CompTox_07March19_WWMetaData.csv",
                                     FALSE)
                if (!is.null(down))
                {
                    setOpts <- c(setOpts,
                                 list(patRoon.path.MetFragCompTox = down))
                }
            }

            if ("MetFrag PubChemLite DB" %in% instWhat)
            {
                down <- downloadFile(instPath, "MetFrag PubChemLite database", "https://zenodo.org/record/4432124/files/PubChemLite_01Jan2021_exposomics.csv",
                                     FALSE)
                if (!is.null(down))
                    setOpts <- c(setOpts, list(patRoon.path.MetFragPubChemLite = down))
            }
            
            if ("SIRIUS" %in% instWhat)
            {
                down <- downloadFile(instPath, "SIRIUS", "https://github.com/boecker-lab/sirius/releases/download/v4.9.12/sirius-4.9.12-win64.zip",
                                     TRUE)
                if (!is.null(down))
                    setOpts <- c(setOpts, list(patRoon.path.SIRIUS = fixPath(file.path(down, "sirius-gui"))))
            }
            
            if ("OpenBabel" %in% instWhat)
            {
                installr::install.URL("https://github.com/openbabel/openbabel/releases/download/openbabel-3-1-1/OpenBabel-3.1.1-x64.exe",
                                      message = FALSE)
            }
            
            if ("BioTransformer" %in% instWhat)
            {
                down <- downloadFile(instPath, "BioTransformer", "https://bitbucket.org/djoumbou/biotransformer/get/master.zip",
                                     TRUE, "biotransformer.zip")
                if (!is.null(down))
                {
                    # rename subdirectory with auto generated name...
                    subDir <- list.files(down, pattern = "^djoumbou\\-biotransformer\\-[[:alnum:]]+$", full.names = TRUE)
                    file.rename(subDir, file.path(down, "biotransformer"))
                    
                    # place in jar from patRoonDeps
                    jar <- downloadFile(file.path(down, "biotransformer"), "BioTransformer jar",
                                        "https://github.com/rickhelmus/patRoonDeps/raw/master/ext/biotransformer-3.0.0.jar",
                                        FALSE)
                    if (!is.null(jar))
                        setOpts <- c(setOpts, list(patRoon.path.BioTransformer = jar))
                }
            }
            
            if ("pngquant" %in% instWhat)
            {
                down <- downloadFile(instPath, "pngquant", "https://pngquant.org/pngquant-windows.zip", TRUE)
                if (!is.null(down))
                    setOpts <- c(setOpts, list(patRoon.path.pngquant = fixPath(file.path(down, "pngquant"))))
            }
            
            if ("ProteoWizard" %in% instWhat &&
                yesNo(title = paste("Due to license agreement restrictions ProteoWizard cannot be installed automatically at this point.",
                                    "Do you want to open the webpage so that you can download and install ProteoWizard manually?",
                                    sep = "\n")))
            {
                utils::browseURL("http://proteowizard.sourceforge.net/download.html")
                while(!yesNo(title = "Did you install ProteoWizard and are ready to continue the patRoon installation?")) {}
                
                pwiz <- findPWizPath()
                while(is.null(pwiz) &&
                      !yesNo("Failed to find ProteoWizard. If you continue now you will have to set the patRoon.pwiz.path option manually using options(). Continue?"))
                {
                    pwiz <- findPWizPath()
                }
                
                if (!is.null(pwiz))
                {
                    cat(paste("\nNOTE: Found ProteoWizard at", pwiz, "\n\n"))
                    setOpts <- c(setOpts, list(patRoon.path.pwiz = pwiz))
                }
            }
        }
        
        return(setOpts)
    },
    
    installPatRoonPackages = function(instPath, exampleData, pkgWhere, force)
    {
        printHeader("Installing patRoon R package(s)...")
        
        curLPaths <- .libPaths(); on.exit(.libPaths(curLPaths))
        .libPaths(getLibPaths(instPath, pkgWhere))
        
        checkPackages("patRoon", pkgWhere, ask = FALSE, type = "gh", repos = "rickhelmus", force = force)
        if (exampleData) # NOTE: patRoonData is always installed from GitHub, as it's too big for patRoonDeps
            checkPackages("patRoonData", "normal", ask = FALSE, type = "gh", repos = "rickhelmus", force = force)
        invisible(NULL)
    }
))()


installPatRoon <- function(what = c("mandatory_packages", "optional_packages", "tools", "deps", "patRoon"),
                           instPath = "~/patRoon-install", exampleData = TRUE,
                           force = FALSE)
{
    if (Sys.info()[["sysname"]] != "Windows" || Sys.info()[["machine"]] != "x86-64")
        stop("Sorry, this script only works on a 64 bit Windows system at the moment.")
    
    validWhat <- c("mandatory_packages", "optional_packages", "tools", "deps", "patRoon")
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
        sep = "\n")
    
    cat("\nNOTE: It is best to run this in a _fresh_ R session!\n")
    
    instPath <- utils$fixPath(instPath)
    
    pkgWhere <- "normal"
    if (any(c("mandatory_packages", "optional_packages", "patRoon") %in% what))
    {
        # UNDONE: add option to always install patRoon from GH?
        
        cat("patRoon and its R package dependencies can either be installed from regular repositories (CRAN, BioConductor)",
            "or from the patRoonDeps repository. The former ensures that you will always get the latest package versions.",
            "On the other hand, using the patRoonDeps repository minimizes the risk on installation errors.",
            "Furthermore, when using the patRoonDeps repository you can choose to use a local R library which is only used",
            "if specified explicitly and will therefore not interfere with your current installed packages.",
            "NOTE: patRoonDeps only work with recent versions of R.\n")
        
        choices <- c(pDepsIso = "Install from patRoonDeps repository (inside an isolated R library)",
                     pDeps = "Install from patRoonDeps repository (inside your default R library)",
                     normal = "Install from standard repositories (CRAN, BioConductor, GitHub)")
        pkgWhere <- menu(choices, title = "(From) Where do you want to install patRoon and its dependencies? If unsure, choose 1.")
        
        if (pkgWhere == 0)
            stop("Aborted")
        
        pkgWhere <- names(choices)[pkgWhere]
    }
    
    jPath <- NULL; setOpts <- list()

    if ("mandatory_packages" %in% what)
        utils$installMandatoryRDeps(instPath, pkgWhere)    
    if ("deps" %in% what)
    {
        # UNDONE: enable Rtools if we allow github installation of patRoon regardless of repos
        jPath <- utils$installExtDeps(instPath, "patRoon" %in% what && pkgWhere == "normal")
    }
    if ("optional_packages" %in% what)
        utils$installOptionalRDeps(instPath, pkgWhere)
    if ("tools" %in% what)
        setOpts <- utils$installTools(instPath)
    
    if (!is.null(jPath) || length(setOpts) > 0 || pkgWhere == "pDepsIso")

    {
        ownRProfPath <- utils$makeOwnRProfile(instPath, setOpts, jPath, pkgWhere == "pDepsIso")
        
        cat("An R script that will set all necessary options to use the installed dependencies was automatically generated.",
            paste("File location:", ownRProfPath),
            "This file must be loaded in each R session prior to using patRoon, e.g. with the following command:",
            sprintf("source(\"%s\")", ownRProfPath),
            sep = "\n")
        
        if (utils$yesNo(paste("The installer can add code to your ~/.Rprofile file so that this file is automatically loaded in each R session.",
                              "Do you want to do this?",
                          sep = "\n")))
            utils$editRProfile(ownRProfPath)
    }
    
    if ("patRoon" %in% what)
        utils$installPatRoonPackages(instPath, exampleData, pkgWhere, force) 
    
    utils$printHeader("All done! You may need to restart R.")
    
    # if (pkgWhere == "pDepsIso")
    #     cat("IMPORTANT: Packages where installed in a separate R library.",
    #         "Please make sure that you add the following line to your script before trying to load patRoon:",
    #         sprintf(".libPaths(c(\"%s\", .libPaths()))", utils$instPathRLib(instPath)))
    
    invisible(NULL)
}

cat("Run installPatRoon() to start the installation\n")
