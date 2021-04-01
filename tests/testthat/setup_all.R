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
