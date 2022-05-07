options(Ncpus = 8)

install.packages(c("devtools", "BiocManager",
                   "dplyr", # missing dep(?)
                   "vdiffr", "desc"))
install.packages(c("circlize", "VennDiagram"),
                 repos="http://cran.r-project.org") # get most recent version (rpm is too old)

BiocManager::install(c("mzR", "xcms", "CAMERA",
                       "Rdisop", # need to install this dep manually (why??)
                       "InterpretMSSpectrum", # for RAMClustR
                       "ropls", # for KPIC2
                       "BiocStyle", "Rgraphviz")) # for MetaClean

remotes::install_github(c("rickhelmus/patRoonData",
                          "thomasp85/farver",
                          "cbroeckl/RAMClustR",
                          "blosloos/enviPick",
                          "blosloos/nontargetData",
                          "blosloos/nontarget",
                          "rickhelmus/KPIC2",
                          "rickhelmus/cliqueMS",
                          "KelseyChetnik/MetaClean", "KelseyChetnik/MetaCleanData"),
                        upgrade = "never")

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
install.packages("cluster")

# check if antyhing failed (ie still not available)
mp <- getMissingPkgs()
if (length(mp) > 0)
{
    cat(sprintf("ERROR: couldn't install these packages: %s\n", paste0(mp, collapse = ", ")))
    q(status = 1)
}
