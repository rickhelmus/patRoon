# disable flags as otherwise GenForm doesn't compile
options(covr.flags = list(CXXFLAGS = '', LDFLAGS = ''))
options(patRoon.progress.opts = list(style = 1))

install.packages(c("covr", "testthat", "vdiffr"))

withr::with_envvar(list(NOT_CRAN = "true"), covr::codecov(quiet = FALSE, errorsAreFatal = FALSE, clean = FALSE))
