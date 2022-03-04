options(Ncpus = 6)

install.packages("devtools")
install.packages("BiocManager")

install.packages("dplyr") # missing dep(?)

BiocManager::install(c("mzR", "xcms", "CAMERA"))
BiocManager::install(c("Rdisop", "InterpretMSSpectrum")) # need to install these deps manually (why??)

remotes::install_github("rickhelmus/patRoonData", upgrade = "never")
remotes::install_github("thomasp85/farver", upgrade = "never")
remotes::install_github("cbroeckl/RAMClustR", upgrade = "never")
remotes::install_github("blosloos/nontargetData", upgrade = "never")
remotes::install_github("blosloos/nontarget", upgrade = "never")
install.packages("vdiffr")
install.packages(c("circlize", "VennDiagram"),
                 repos="http://cran.r-project.org") # get most recent version (rpm is too old)

BiocManager::install("ropls") # for KPIC2
remotes::install_github("rickhelmus/KPIC2")
remotes::install_github("rickhelmus/cliqueMS")
BiocManager::install(c("BiocStyle", "Rgraphviz")) # for MetaClean
remotes::install_github("KelseyChetnik/MetaClean")
remotes::install_github("KelseyChetnik/MetaCleanData")

install.packages("desc")

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
