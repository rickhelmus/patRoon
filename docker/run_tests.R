options(patRoon.path.metFragCL = "~/MetFrag2.4.3-CL.jar",
        patRoon.path.SIRIUS = "~/sirius-linux64-headless-4.0/bin")

devtools::install(upgrade_dependencies = FALSE)

# return failure exit code when tests fail: https://github.com/r-lib/testthat/issues/515
tret <- as.data.frame(devtools::test())
if (sum(tret$failed) > 0 || any(tret$error))
    q(status = 1)
