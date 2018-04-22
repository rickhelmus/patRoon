options(patRoon.path.metFragCL = "~/MetFrag2.4.3-CL.jar",
        patRoon.path.SIRIUS = "~/sirius-linux64-headless-4.0/bin",
        patRoon.path.vdiffr = "docker")

devtools::install(upgrade_dependencies = FALSE)

# unfortunately vdiffr doesn't allow to specify the deps file name.
file.rename("tests/figs/deps-docker.txt", "tests/figs/deps.txt")

# return failure exit code when tests fail: https://github.com/r-lib/testthat/issues/515
tret <- as.data.frame(devtools::test(reporter = JunitReporter$new(file = "junit.xml")))
print(tret)
if (sum(tret$failed) > 0 || any(tret$error))
    q(status = 1)
