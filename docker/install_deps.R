options(Ncpus = 8)
options(timeout = 300)

install.packages("remotes")
remotes::install_github("rickhelmus/patRoonInst")

# install everything, including patRoon. This will be removed afterward, but is handy to get all deps.
patRoonInst::install(ask = FALSE, quiet = FALSE)
# UNDONE: pull in other deps
# ---
install.packages("reticulate")
BiocManager::install("ChemmineR")
BiocManager::install('fmcsR')
remotes::install_github("rickhelmus/patRoonData@version30")
remotes::install_github("rickhelmus/patRoonDataIMS")
remotes::install_github("rickhelmus/Rmstoolkitlib")
remotes::install_github("rickhelmus/patRoonExt@version20")
remotes::install_github("rickhelmus/patRoon@version30")
# ---
remove.packages("patRoon")
