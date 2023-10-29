options(Ncpus = 8)

install.packages("patRoonInst", repos = c('https://rickhelmus.r-universe.dev', 'https://cloud.r-project.org'))

# install everything, including patRoon. This will be removed afterward, but is handy to get all deps.
patRoonInst::install(pkgs="patRoonExt", ask = FALSE, quiet = FALSE)
patRoonInst::install(ask = FALSE, quiet = FALSE)
remove.packages("patRoon")
