unlink(getWorkPath(), TRUE)
dir.create(getWorkPath())

if (getOption("patRoon.clearCache", FALSE))
    clearCache("all")

if (!file.exists(getTestDataPath()))
    dir.create(getTestDataPath())
