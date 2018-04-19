getWorkPath <- function(file = "", ...) if (nzchar(file)) file.path("test_temp", file, ...) else "test_temp"
getTestDataPath <- function() "test_data"
testFile <- function(f, ..., text = FALSE) file.path(getTestDataPath(), paste0(f, ..., if (!text) ".Rds" else ".txt", collapse = ""))
getTestAnaInfo <- function() generateAnalysisInfo(patRoonData::exampleDataPath(),
                                                  groups = c(rep("solvent", 3), rep("standard", 3)),
                                                  refs = "solvent")
getTestFGroups <- function(anaInfo = getTestAnaInfo()) groupFeatures(findFeatures(anaInfo, "openms", logPath = NULL), "openms")

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

expect_known_show <- function(object, file)
{
    act <- quasi_label(rlang::enquo(object))
    
    text <- capture_output_lines(show(act$val))
    
    # remove last line with object size as it may vary even if object remain the same
    text <- text[seq_len(length(text)-1)]
    
    text <- paste0(text, collapse = "")
    
    # based on simplified expect_known_output
    if (!file.exists(file))
    {
        warning("Creating reference output", call. = FALSE)
        cat(text, file = file)
        succeed()
    }
    else
    {
        ref <- patRoon:::readAllFile(file)
        cat(text, file = file)
        
        cmp <- compare(text, enc2native(ref))
        expect(cmp$equal, sprintf("show reference of %s has changed\n%s", act$lab, cmp$message))
    }
    
    invisible(act$val)
}

expect_plot <- function(object)
{
    tf <- tempfile()
    withr::with_png(tf, act <- quasi_label(rlang::enquo(object)))
    expect(file.exists(tf), "failed to generate plot")
    invisible(act$val)
}

expect_docker <- function(...) vdiffr::expect_doppelganger(..., path = getOption("patRoon.path.vdiffr"))
