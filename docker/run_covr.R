# disable flags as otherwise GenForm doesn't compile
options(covr.flags = list(CXXFLAGS = '', LDFLAGS = ''))

Sys.setenv(TESTTHAT_CPUS = 2)
Sys.setenv(PATROON_MP_MAXPROCS = 2)
Sys.setenv(PKG_BUILD_EXTRA_FLAGS = "false")

# HACK: trigger compile first. It seems that parallel testing triggers multiple compiles, resulting in random compile
# errors.
install.packages("pkgload")
pkgload::load_all()

install.packages(c("testthat", "vdiffr"))
remotes::install_github("rickhelmus/covr@live-console-update")

withr::with_envvar(list(NOT_CRAN = "true"), covr::codecov(quiet = FALSE, errorsAreFatal = FALSE, clean = FALSE,
                                                          type = "none", code = 'testthat::test_package("patRoon")',
                                                          code_stdout = TRUE))
