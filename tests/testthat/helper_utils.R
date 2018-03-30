getWorkPath <- function() "test_temp"
getTestDataPath <- function() "test_data"
testFile <- function(f, ..., text = FALSE) file.path(getTestDataPath(), paste0(f, ..., if (!text) ".Rds" else ".txt", collapse = ""))

expect_file <- function(object, file, removeIfExists = TRUE)
{
    if (removeIfExists && file.exists(file))
        file.remove(file)
    
    act <- quasi_label(rlang::enquo(object))
    expect(file.exists(file), sprintf("failed to generate %s", file))
    invisible(act$val)
}

expect_range <- function(object, r)
{
    act <- quasi_label(rlang::enquo(object))
    act$r <- range(act$val)
    expect(act$r[1] >= r[1] && act$r[2] <= r[2],
           sprintf("range of %s is %.1f - %.1f which is outside %.1f - %.1f",
                   act$lab, act$r[1], act$r[2], r[1], r[2]))
    invisible(act$val)
}
