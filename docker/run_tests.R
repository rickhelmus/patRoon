# unfortunately vdiffr doesn't allow to specify the deps file name.
# file.rename("tests/figs/deps-docker.txt", "tests/figs/deps.txt")

install.packages(c("devtools", "vdiffr"))

# NOTE: these things need to be set as env vars since parallel testthat seems to ignore options() and .Rprofile
Sys.setenv(TESTTHAT_CPUS = 2)
Sys.setenv(PATROON_MP_MAXPROCS = 2)
Sys.setenv(PATROON_THREADS = 2)
Sys.setenv(PKG_BUILD_EXTRA_FLAGS = "false")

# HACK: do this until https://github.com/Rdatatable/data.table/issues/7749 is on CRAN
data.table::update_dev_pkg()

# HACK: trigger compile first. It seems that parallel testing triggers multiple compiles, resulting in random compile
# errors.
devtools::load_all()

# disabled: consumes too much RAM on CircleCI
# SIRIUSAPI <- patRoon:::startSIRIUS() # HACK start it now so we can share it between tests and makes things faster.
SIRIUSAPI <- NULL
# --> start SIRIUS manually with limited instances/cores
SIRProc <- processx::process$new(patRoon:::getExtDepPath("sirius"), c("--cores", "1", "--buffer", "1", "REST", "-s", "--headless"))

if (!is.null(Sys.getenv("PATROON_SIRUSER")) && nzchar(Sys.getenv("PATROON_SIRUSER")) &&
    !is.null(Sys.getenv("PATROON_SIRPASS")) && nzchar(Sys.getenv("PATROON_SIRPASS")))
{
    SIRIUSLogin(login = c(username = Sys.getenv("PATROON_SIRUSER"), password = Sys.getenv("PATROON_SIRPASS")),
                SIRIUSAPI = SIRIUSAPI)
}

# return failure exit code when tests fail: https://github.com/r-lib/testthat/issues/515
tret <- as.data.frame(devtools::test(reporter = testthat::MultiReporter$new(list(testthat::SummaryReporter$new(),
                                                                                 testthat::JunitReporter$new(file = "~/junit.xml")))))
print(tret)

if (SIRProc$is_alive())
{
    SIRProc$interrupt()
    Sys.sleep(3)
    if (SIRProc$is_alive())
        SIRProc$kill()
}

if (sum(tret$failed) > 0 || any(tret$error))
    q(status = 1)
