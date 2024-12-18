# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

# make sure reuse is installed: https://reuse.readthedocs.io/en/stable/readme.html

local({
    annotateFiles <- function(paths, exclude = NULL, style = NULL)
    {
        if (!is.null(exclude))
            paths <- setdiff(paths, exclude)
        args <- c("annotate",
                  "--year=2016-2024",
                  "--copyright=Rick Helmus <r.helmus@uva.nl>",
                  "--license=GPL-3.0-only")
        if (!is.null(style))
            args <- c(args, paste0("--style=", style))
        system2("reuse", shQuote(c(args, paths)), stdout = "", stderr = "")
    }
    
    # Rcpp is auto-generated
    annotateFiles(Sys.glob("R/*.R"), "R/RcppExports.R")
    
    # Rcpp is auto-generated, force C style for consistency between .h/.cpp files
    annotateFiles(c(Sys.glob("src/*.h"), Sys.glob("src/*.cpp")), "src/RcppExports.cpp", style = "c")
    
    annotateFiles(c("src/Makevars", "src/Makevars.win"), style = "python") # not exactly Python, but something that uses #
    
    # we ignore the .Rmd files for now, as auto annotation messes up the files
    annotateFiles(c(Sys.glob("inst/js/*.js"), "inst/misc/runSAFD.jl"))
    
    annotateFiles(c("tests/testthat.R", "tests/testthat/*.R"))
    
    annotateFiles(c("tools/update-license-header.R", "data-raw/prepare-data.R"))
    
    # verify --> REUSE.toml lists licensisng for other files
    system2("reuse", "lint", stdout = "", stderr = "")
})
