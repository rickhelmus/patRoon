unlink(getWorkPath(), TRUE)
dir.create(getWorkPath())

runWithEnviPick <- FALSE

testAnaInfo <- generateAnalysisInfo(patRoonData::exampleDataPath(),
                                    groups = c(rep("solvent", 3), rep("standard", 3)),
                                    refs = "solvent")

if (runWithEnviPick)
{
    exDataFiles <- list.files(patRoonData::exampleDataPath(), "\\.mzML$", full.names = TRUE)
    convertMSFiles(testAnaInfo, "mzML", "mzXML", getWorkPath())
    
    testAnaInfoMzXML <- generateAnalysisInfo(getWorkPath(),
                                             groups = c(rep("solvent", 3), rep("standard", 3)),
                                             refs = "solvent")
}

if (!file.exists(getTestDataPath()))
    dir.create(getTestDataPath())
