# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

if (!nzchar(Sys.getenv("PATROON_NOTESTS")))
{
    library(testthat)
    library(patRoon)
    
    envOpts <- Sys.getenv(c("PATROON_OPENMS", "PATROON_OBABEL", "PATROON_METFRAG", "PATROON_SIRIUS",
                            "PATROON_BIOTRANSFORMER", "PATROON_PCLITE"))
    if (nzchar(envOpts[["PATROON_OPENMS"]]))
        options(patRoon.path.OpenMS = envOpts[["PATROON_OPENMS"]])
    if (nzchar(envOpts[["PATROON_OBABEL"]]))
        options(patRoon.path.obabel = envOpts[["PATROON_OBABEL"]])
    if (nzchar(envOpts[["PATROON_METFRAG"]]))
        options(patRoon.path.MetFragCL = envOpts[["PATROON_METFRAG"]])
    if (nzchar(envOpts[["PATROON_SIRIUS"]]))
        options(patRoon.path.SIRIUS = envOpts[["PATROON_SIRIUS"]])
    if (nzchar(envOpts[["PATROON_BIOTRANSFORMER"]]))
        options(patRoon.path.BioTransformer = envOpts[["PATROON_BIOTRANSFORMER"]])
    if (nzchar(envOpts[["PATROON_PCLITE"]]))
        options(patRoon.path.MetFragPubChemLite = envOpts[["PATROON_PCLITE"]])
    
    options(patRoon.clearCache = TRUE)
    options(patRoon.progress.opts = list(style = 1))
    
    # https://github.com/r-lib/devtools/issues/1526
    Sys.unsetenv("R_TESTS")
    
    ju <- Sys.getenv("PATROON_JUNIT")
    if (nzchar(ju))
        test_check("patRoon", reporter = JunitReporter$new(file = ju))
    else
        test_check("patRoon", reporter = "summary")
}
