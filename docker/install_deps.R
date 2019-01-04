options(Ncpus = 3)

install.packages("devtools")

source("https://bioconductor.org/biocLite.R")
biocLite(c("mzR", "xcms", "CAMERA"))
biocLite(c("Rdisop", "InterpretMSSpectrum")) # need to install these deps manually (why??)

devtools::install_github("cbroeckl/RAMClustR@0370902", upgrade_dependencies = FALSE)
devtools::install_github("rickhelmus/patRoonData", upgrade_dependencies = FALSE)
# devtools::install_github("yutannihilation/vdiffr", ref = "fix-strip-version", upgrade_dependencies = FALSE)
# devtools::install_github("lionel-/freetypeharfbuzz", upgrade_dependencies = FALSE)
# devtools::install_github("lionel-/vdiffr", upgrade_dependencies = FALSE)
install.packages("vdiffr")
devtools::install_github("mllg/checkmate", upgrade_dependencies = FALSE)

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
