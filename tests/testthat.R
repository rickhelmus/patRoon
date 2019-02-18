library(testthat)
library(patRoon)

envOpts <- Sys.getenv(c("PATROON_METFRAG", "PATROON_SIRIUS"))
if (nzchar(envOpts[["PATROON_METFRAG"]]))
    options(patRoon.path.MetFragCL = envOpts[["PATROON_METFRAG"]])
if (nzchar(envOpts[["PATROON_SIRIUS"]]))
    options(patRoon.path.SIRIUS = envOpts[["PATROON_SIRIUS"]])

options(patRoon.clearCache = TRUE)

# https://github.com/r-lib/devtools/issues/1526
Sys.unsetenv("R_TESTS")

ju <- Sys.getenv("PATROON_JUNIT")
if (nzchar(ju))
{
    test_check("patRoon", reporter = JunitReporter$new(file = ju))
} else
    test_check("patRoon")
