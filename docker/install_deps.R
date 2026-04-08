options(Ncpus = 8)
options(timeout = 300)

install.packages("remotes")
remotes::install_github("rickhelmus/patRoonInst")

# install everything, including patRoon. This will be removed afterward, but is handy to get all deps.
patRoonInst::install(ask = FALSE, quiet = FALSE)
patRoon::installTIMSCONVERT()
patRoon::installC3SDB()
# remove.packages("patRoon")
