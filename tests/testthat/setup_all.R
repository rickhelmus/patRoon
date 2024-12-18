# SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
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

options(datatable.auto.index = FALSE) # should make tests more reproducible
options(patRoon.MP.logPath = FALSE)

# HACK: sometimes plot isn't recognized as an S4 generic
# setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
