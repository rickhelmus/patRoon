# disable flags as otherwise GenForm doesn't compile
options(covr.flags = list(CXXFLAGS = '', LDFLAGS = ''))

Sys.setenv(TESTTHAT_CPUS = 2)
Sys.setenv(PATROON_MP_MAXPROCS = 2)
Sys.setenv(PATROON_THREADS = 2)
Sys.setenv(PKG_BUILD_EXTRA_FLAGS = "false")

# HACK: do this until https://github.com/Rdatatable/data.table/issues/7749 is on CRAN
data.table::update_dev_pkg()

# HACK: trigger compile first. It seems that parallel testing triggers multiple compiles, resulting in random compile
# errors.
install.packages("pkgload")
pkgload::load_all()

install.packages(c("testthat", "vdiffr"))
remotes::install_github("rickhelmus/covr@live-console-update")

SIRIUSAPI <- patRoon:::startSIRIUS() # HACK start it now so we can share it between tests and makes things faster.

if (!is.null(Sys.getenv("PATROON_SIRUSER")) && nzchar(Sys.getenv("PATROON_SIRUSER")) &&
    !is.null(Sys.getenv("PATROON_SIRPASS")) && nzchar(Sys.getenv("PATROON_SIRPASS")))
{
    SIRIUSLogin(login = c(username = Sys.getenv("PATROON_SIRUSER"), password = Sys.getenv("PATROON_SIRPASS")),
                SIRIUSAPI = SIRIUSAPI)
}

withr::with_envvar(list(NOT_CRAN = "true"), covr::codecov(quiet = FALSE, errorsAreFatal = FALSE, clean = FALSE,
                                                          type = "none", code = 'testthat::test_package("patRoon")',
                                                          code_stdout = TRUE))
