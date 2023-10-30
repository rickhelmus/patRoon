options(Ncpus = 8)

install.packages("remotes")
remotes::install_github("rickhelmus/patRoonInst")

# install everything, including patRoon. This will be removed afterward, but is handy to get all deps.
patRoonInst::install(ask = FALSE, quiet = FALSE)
remove.packages("patRoon")
