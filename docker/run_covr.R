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

# disabled: consumes too much RAM on CircleCI
# SIRIUSAPI <- patRoon:::startSIRIUS() # HACK start it now so we can share it between tests and makes things faster.
SIRIUSAPI <- NULL

if (!is.null(Sys.getenv("PATROON_SIRUSER")) && nzchar(Sys.getenv("PATROON_SIRUSER")) &&
    !is.null(Sys.getenv("PATROON_SIRPASS")) && nzchar(Sys.getenv("PATROON_SIRPASS")))
{
    SIRIUSLogin(login = c(username = Sys.getenv("PATROON_SIRUSER"), password = Sys.getenv("PATROON_SIRPASS")),
                SIRIUSAPI = SIRIUSAPI)
}

# --> start SIRIUS manually with limited instances/cores
# HACK: do after login, as that call quits SIRIUS
SIRProc <- processx::process$new(patRoon:::getExtDepPath("sirius"), c("--cores", "1", "--buffer", "1", "REST", "-s", "--headless"))
Sys.sleep(3) # give SIRIUS some time to start up

withr::with_envvar(list(NOT_CRAN = "true"), covr::codecov(quiet = FALSE, errorsAreFatal = FALSE, clean = FALSE,
                                                          type = "none", code = 'testthat::test_package("patRoon")',
                                                          code_stdout = TRUE))

if (SIRProc$is_alive())
{
    SIRProc$interrupt()
    Sys.sleep(3)
    if (SIRProc$is_alive())
        SIRProc$kill()
}
