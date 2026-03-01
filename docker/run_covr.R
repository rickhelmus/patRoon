# disable flags as otherwise GenForm doesn't compile
options(covr.flags = list(CXXFLAGS = '', LDFLAGS = ''))

Sys.setenv(TESTTHAT_CPUS = 3)

install.packages(c("testthat", "vdiffr"))
remotes::install_github("rickhelmus/covr@live-console-update")

withr::with_envvar(list(NOT_CRAN = "true"), covr::codecov(quiet = FALSE, errorsAreFatal = FALSE, clean = FALSE,
                                                          type = "none", code = 'testthat::test_package("patRoon")',
                                                          code_stdout = TRUE))
