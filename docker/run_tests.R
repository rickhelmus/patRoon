# unfortunately vdiffr doesn't allow to specify the deps file name.
# file.rename("tests/figs/deps-docker.txt", "tests/figs/deps.txt")

install.packages(c("devtools", "vdiffr"))

# NOTE: these things need to be set as env vars since parallel testthat seems to ignore options() and .Rprofile
Sys.setenv(TESTTHAT_CPUS = 3)
Sys.setenv(PATROON_MP_MAXPROCS = 2)
Sys.setenv(PKG_BUILD_EXTRA_FLAGS = "false")

# return failure exit code when tests fail: https://github.com/r-lib/testthat/issues/515
tret <- as.data.frame(devtools::test(reporter = testthat::MultiReporter$new(list(testthat::SummaryReporter$new(),
                                                                                 testthat::JunitReporter$new(file = "~/junit.xml")))))
print(tret)
if (sum(tret$failed) > 0 || any(tret$error))
    q(status = 1)
