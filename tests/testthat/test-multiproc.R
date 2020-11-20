context("multiproc")

# handy tool from processx
px <- paste0(
    system.file(package = "processx", "bin", "px"),
    system.file(package = "processx", "bin", .Platform$r_arch, "px.exe")
)

# test: finish, error and timeout handlers, commands are actually executed, retries work
# finish handler: return TRUE or index and check afterwards if all were written
# error/timeout: as finish
# commands are executed: create file in command queue and check for its existance afterwards
# retries: check if it can succeed after handler is called (e.g. by ensuring in handler command can actually run)
# further: different amount of maxprocs? check output/error (is this actually used?)

performMPTest <- function(len, fail = numeric(), timeout = numeric(), failFirst = FALSE, delay = FALSE,
                          errorHandler = function(...) FALSE, timeoutHandler = function(...) FALSE,
                          maxProcAmount = NULL, ...)
{
    if (!is.null(maxProcAmount))
        withr::local_options(list(patRoon.MP.maxProcs = maxProcAmount))
    
    if (Sys.info()[["sysname"]] == "Windows")
    {
        shFile <- tempfile("makefile", fileext = ".bat")
        cat("if not exist %1 (SET ret=1) else (SET ret=0)",
            "echo y > %1",
            if (delay) paste(paste0("\"", px, "\""), "sleep 2") else "",
            if (failFirst) "exit /b %ret%" else "", # will fail if file didn't exist yet (i.e. first time it's executed)
            sep = "\n", file = shFile)
    }
    else
    {
        shFile <- tempfile("makefile", fileext = ".sh")
        cat("#!/bin/sh",
            "test -f $1",
            "ret=$?",
            "echo y > $1",
            if (delay) paste(paste0("\"", px, "\""), "sleep 2") else "",
            if (failFirst) "exit $ret" else "", # will fail if file didn't exist yet (i.e. first time it's executed)
            sep = "\n", file = shFile)
        Sys.chmod(shFile)
    }

    cmdQueue <- setNames(lapply(seq_len(len), function(i)
    {
        if (i %in% fail)
            return(list(command = px, args = c("return", "2")))
        else if (i %in% timeout)
            return(list(command = px, args = c("sleep", 60)))

        testFile <- tempfile("mptest")
        list(command = shFile, args = testFile, outFile = testFile, index = i)
    }), as.character(seq_len(len)))

    ct <- curTimeMS()
    res <- executeMultiProcess(cmdQueue, function(cmd) { cmd }, procTimeout = 5, showProgress = FALSE,
                               errorHandler = errorHandler, timeoutHandler = timeoutHandler, ...)
    printf("time: %.2f\n", curTimeMS() - ct)

    if (failFirst)
        fail <- timeout <- numeric()

    succeeded <- setdiff(seq_len(len), c(fail, timeout))
    if (length(succeeded) > 0)
    {
        expect_false(any(sapply(res[succeeded], is.null)))
        succeeded <- succeeded[!sapply(res[succeeded], is.null)] # clear out NULL otherwise next tests give errors
        expect_true(all(file.exists(sapply(res[succeeded], "[[", "outFile"))))
        expect_equal(sapply(cmdQueue[succeeded], "[[", "index"), sapply(res[succeeded], "[[", "index"))
    }

    failed <- c(fail, timeout)
    if (length(failed) > 0)
        expect_true(all(sapply(res[failed], is.null)))
}

ehandler <- function(cmd, exitStatus, retries) retries == 0
thandler <- function(cmd, retries) retries == 0

test_that("multi-process functionality", {
    performMPTest(10)
    performMPTest(10, batchSize = 3)
    performMPTest(10, batchSize = 8)
    performMPTest(10, batchSize = 20)
    performMPTest(10, maxProcAmount = 1) # maxProcAmount>1 by default
    performMPTest(3, failFirst = TRUE, errorHandler = ehandler)

    performMPTest(10, fail = c(2, 8))
    performMPTest(10, fail = c(2, 8), batchSize = 3)
    performMPTest(10, fail = c(2, 8), maxProcAmount = 1)
    performMPTest(10, fail = c(2, 8), errorHandler = ehandler)
    performMPTest(10, fail = c(1, 10))
    performMPTest(2, fail = c(1, 2))
    performMPTest(4, fail = c(1, 2), batchSize = 2)

    performMPTest(10, timeout = c(2, 8))
    performMPTest(10, timeout = c(2, 8), batchSize = 3)
    performMPTest(10, timeout = c(2, 8), maxProcAmount = 1)
    performMPTest(10, timeout = c(2, 8), timeoutHandler = thandler)
    performMPTest(10, timeout = c(1, 10))
    performMPTest(2, timeout = c(1, 2))
    performMPTest(4, timeout = c(1, 2), batchSize = 2)
    performMPTest(4, delay = TRUE) # delay is below timeout, shouldn't fail
    performMPTest(4, delay = TRUE, batchSize = 4, maxProcAmount = 1) # batch process will take longer
})

test_that("multi-process future functionality", {
    withr::local_options(list(patRoon.MP.method = "future"))
    future::plan("multisession", workers = 2)
    withr::defer(future::plan("sequential"))
                        
    performMPTest(10)
    performMPTest(3, failFirst = TRUE, errorHandler = ehandler)
    
    performMPTest(10, fail = c(2, 8))
    performMPTest(10, fail = c(2, 8), errorHandler = ehandler)
    performMPTest(10, fail = c(1, 10))
    performMPTest(2, fail = c(1, 2))
    
    performMPTest(10, timeout = c(2, 8))
    performMPTest(10, timeout = c(2, 8), timeoutHandler = thandler)
    performMPTest(10, timeout = c(1, 10))
    performMPTest(2, timeout = c(1, 2))
    performMPTest(4, delay = TRUE) # delay is below timeout, shouldn't fail
})
