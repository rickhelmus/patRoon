# executed by Makevars to get compilation flags

getFlags <- function(what)
{
    # check if mstoolkit is installed
    if (requireNamespace("Rmstoolkitlib", quietly = TRUE))
    {
        if (what == "PKG_CXXFLAGS")
            cat("-DWITH_MSTK ")
        cat(Rmstoolkitlib::pkgconfig(what))
    }
    # NOTE: the use of an underscore or dash for the CPU arch seems a bit random so just check for both
    if (what %in% c("PKG_CFLAGS", "PKG_CXXFLAGS") && Sys.info()["machine"] %in% c("x86-64", "x86_64") &&
        Sys.info()["sysname"] %in% c("Windows", "Linux"))
        cat(" -DWITH_OTIMS")
}
