options(Ncpus = 3)

install.packages("devtools")
install.packages("BiocManager")

BiocManager::install(c("mzR", "xcms", "CAMERA"))
BiocManager::install(c("Rdisop", "InterpretMSSpectrum")) # need to install these deps manually (why??)

remotes::install_github("rickhelmus/patRoonData", upgrade = "never")
remotes::install_github("thomasp85/farver", upgrade = "never")
remotes::install_github("cbroeckl/RAMClustR", upgrade = "never")
install.packages("vdiffr")
install.packages("circlize", repos="http://cran.r-project.org") # get most recent version (rpm is too old)

install.packages("desc")
dp <- desc::desc_get_deps()
dp <- dp[dp$type != "Suggests", "package"]
dp <- dp[!sapply(dp, require, quietly = TRUE, character.only = TRUE)]
if (length(dp) > 0)
{
    cat(sprintf("Missing packages: %s\n", paste0(dp, collapse = ", ")))
    install.packages(dp)
}

# workaround for sildist not found bug: https://github.com/hemberg-lab/SC3/issues/36
install.packages("cluster")
