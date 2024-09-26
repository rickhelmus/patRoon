# executed by Makevars to get compilation flags

getFlags <- function(what)
{
    # check if mstoolkit is installed
    if (requireNamespace("Rmstoolkitlib", quietly = TRUE))
    {
        if (what == "PKG_CXXFLAGS")
            cat("-DWITH_MSTK ")
        cat(Rmstoolkitlib::pkgconfig(what), sep = " ")
    }
    else
        cat("") # nothing special
}
