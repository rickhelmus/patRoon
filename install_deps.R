install.packages("rJava")
install.packages("magick")
install.packages("devtools")
# install.packages("vdiffr")

# update
install.packages("DT")

source("https://bioconductor.org/biocLite.R")
biocLite(c("mzR", "xcms", "CAMERA"))
biocLite(c("Rdisop", "InterpretMSSpectrum")) # need to install these deps manually (why??)

devtools::install_github("cbroeckl/RAMClustR@036cdd4")
devtools::install_github("rickhelmus/patRoonData")
devtools::install_github("yutannihilation/vdiffr", ref = "fix-strip-version")
