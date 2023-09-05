options(Ncpus = 8)

# pre-install some generally needed R packages
install.packages(c("devtools", "BiocManager", "vdiffr", "desc"))

depFile <- tempfile(fileext = ".R")
stopifnot(download.file("https://raw.githubusercontent.com/rickhelmus/patRoonDeps/master/Rdeps.R", depFile) == 0)
source(depFile)
dependencies <- getRDependencies(Sys.getenv("GH_BRANCH", "master"), "linux")
dependencies <- dependencies[names(dependencies) != "patRoon"]

installDeps <- function(deps)
{
    for (pkg in names(deps))
    {
        md <- deps[[pkg]]
        
        if (!is.null(md[["deps"]]))
            installDeps(md$deps)
        
        if (md$type == "cran")
            install.packages(pkg)
        else if (md$type == "bioc")
            BiocManager::install(pkg)
        else # gh
        {
            repos <- if (!is.null(md[["repos"]])) md$repos else pkg
            ref <- if (!is.null(md[["branch"]]))
                md$branch
            else if (!is.null(md[["tag"]]))
                md$tag
            else if (!is.null(md[["commit"]]))
                md$commit
            else
                "master"
            remotes::install_github(paste0(md$user, "/", repos), ref = ref, subdir = md[["pkgroot"]], upgrade = "never")
        }
    }
}

installDeps(dependencies)

getMissingPkgs <- function()
{
    dp <- desc::desc_get_deps()
    dp <- dp[dp$type != "Suggests", "package"]
    dp <- dp[dp != "R"]
    return(dp[!sapply(dp, require, quietly = TRUE, character.only = TRUE)])
}

mp <- getMissingPkgs()
if (length(mp) > 0)
{
    cat(sprintf("Missing packages: %s\n", paste0(mp, collapse = ", ")))
    install.packages(mp)
}

# workaround for sildist not found bug: https://github.com/hemberg-lab/SC3/issues/36
# install.packages("cluster")

# check if anything failed (ie still not available)
mp <- getMissingPkgs()
if (length(mp) > 0)
{
    cat(sprintf("ERROR: couldn't install these packages: %s\n", paste0(mp, collapse = ", ")))
    q(status = 1)
}
