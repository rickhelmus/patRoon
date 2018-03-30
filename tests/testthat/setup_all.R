unlink(getWorkPath(), TRUE)
dir.create(getWorkPath())

exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
file.copy(exDataFiles, getWorkPath())

testAnaInfo <- generateAnalysisInfo(getWorkPath(),
                                    groups = c(rep("solvent", 3), rep("standard", 3)),
                                    refs = "solvent")

convertMSFiles(testAnaInfo, "mzML", "mzXML")

if (!file.exists(getTestDataPath()))
    dir.create(getTestDataPath())
