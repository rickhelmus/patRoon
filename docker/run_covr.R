# disable flags as otherwise GenForm doesn't compile
options(covr.flags = list(CXXFLAGS = '', LDFLAGS = ''))
options(patRoon.progress.opts = list(style = 1))

install.packages("covr")

withr::with_envvar(c(PATROON_METFRAG = getOption("patRoon.path.MetFragCL"),
                     PATROON_SIRIUS = getOption("patRoon.path.SIRIUS"),
                     PATROON_BIOTRANSFORMER = getOption("patRoon.path.BioTransformer")),
                   covr::codecov(quiet = FALSE, errorsAreFatal = FALSE, clean = FALSE))
