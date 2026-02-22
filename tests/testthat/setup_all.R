# SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

unlink(getWorkPath(), TRUE)
dir.create(getWorkPath())

if (getOption("patRoon.clearCache", FALSE))
    clearCache("all")

if (!file.exists(getTestDataPath()))
    dir.create(getTestDataPath())

if (doDATests())
{
    revertDAAnalyses(getDAAnaInfo())
    clearCache("featuresBruker")
}

withr::local_options(list(
    datatable.auto.index = FALSE,
    patRoon.MP.logPath = FALSE,
    patRoon.MS.backends = c("mzr", "mstoolkit"),
    pkg.build_extra_flags = FALSE,
    patRoon.progress.opts = list(style = 1),
    patRoon.MP.maxProcs = Sys.getenv("PATROON_MP_MAXPROCS", getOption("patRoon.MP.maxProcs", 2L))
), .local_envir = teardown_env())

if (nzchar(Sys.getenv("PATROON_OPENMS")))
    withr::local_options(list(patRoon.path.OpenMS = Sys.getenv("PATROON_OPENMS")), .local_envir = teardown_env())

# HACK: sometimes plot isn't recognized as an S4 generic
# setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
