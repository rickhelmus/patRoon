getWorkPath <- function() "test_temp"
getTestDataPath <- function() "test_data"
testFile <- function(f, ..., text = FALSE) file.path(getTestDataPath(), paste0(f, ..., if (!text) ".Rds" else ".txt", collapse = ""))