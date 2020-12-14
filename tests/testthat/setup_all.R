unlink(getWorkPath(), TRUE)
dir.create(getWorkPath())

if (testWithSets())
{
    path <- file.path(getWorkPath(), "set2")
    dir.create(path)
    file.copy(Sys.glob(file.path(patRoonData::exampleDataPath(), "*.mzML")), path)
    anaFiles <- Sys.glob(file.path(path, "*.mzML"))
    newAnaFiles <- file.path(dirname(anaFiles), paste0("set2-", basename(anaFiles)))
    file.rename(anaFiles, newAnaFiles)
}

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
