# executed by Makevars to get compilation flags

MSTKAvailable <- function() requireNamespace("Rmstoolkitlib", quietly = TRUE)

# NOTE: the use of an underscore or dash for the CPU arch seems a bit random so just check for both
OTIMSAvailable <- function() Sys.info()["machine"] %in% c("x86-64", "x86_64") && Sys.info()["sysname"] %in% c("Windows", "Linux")

getFlags <- function(what)
{
    if (MSTKAvailable())
    {
        if (what == "PKG_CXXFLAGS")
            cat("-DWITH_MSTK ")
        cat(Rmstoolkitlib::pkgconfig(what))
    }
    
    if (what %in% c("PKG_CFLAGS", "PKG_CXXFLAGS") && OTIMSAvailable())
        cat(" -DWITH_OTIMS")
}

printConfig <- function()
{
    cat("------------\npatRoon build configuration\n")
    avail <- if (MSTKAvailable())
        paste("yes - version", packageVersion("Rmstoolkitlib"))
    else
        "NOT AVAILABLE! (did you install Rmstoolkitlib?)"
    cat("mstoolkit available: ", avail, "\n")
    avail <- if (OTIMSAvailable())
        "yes"
    else
        "not available! (only on x86-64 Windows/Linux)"
    cat("OpenTIMS available: ", avail, "\n")
    cat("------------\n")
}
