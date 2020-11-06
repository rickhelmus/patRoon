executeMultiProcessFuture <- function(commandQueue, finishHandler, timeoutHandler, errorHandler,
                                      prepareHandler, procTimeout, printOutput, printError, logSubDir,
                                      showProgress, batchSize = 1, ...)
{
    results <- future.apply::future_lapply(commandQueue, function(cmd)
    {
        if (!is.null(prepareHandler))
            cmd <- prepareHandler(cmd) # UNDONE
        timeoutRetries <- errorRetries <- 0
        while (TRUE)
        {
            stat <- processx::run(cmd$command, cmd$args, error_on_status = FALSE,
                                  timeout = procTimeout, cleanup_tree = TRUE)
            
            if (stat$timeout)
            {
                if (timeoutHandler(cmd = cmd, retries = errorRetries))
                {
                    timeoutRetries <- timeoutRetries + 1
                    next
                }
            }
            else if (stat$status != 0)
            {
                if (errorHandler(cmd = cmd, exitStatus = stat$status, retries = errorRetries))
                {
                    errorRetries <- errorRetries + 1
                    next
                }
            }
            else # success
            {
                res <- finishHandler(cmd)
                return(list(stdout = stat$stdout, stderr = stat$stderr, result = res))
            }
            
            # failure
            return(list(stdout = stat$stdout, stderr = stat$stderr, result = NULL))
        }
    }, future.scheduling = 1.0)

    logPath <- getOption("patRoon.logPath", FALSE)
    if (!is.null(logSubDir) && !isFALSE(logPath))
    {
        logPath <- file.path(logPath, logSubDir)
        mkdirp(logPath)

        for (i in seq_along(commandQueue))
        {
            if (!is.null(commandQueue[[i]]$logFile))
            {
                tryCatch({
                    fprintf(file.path(logPath, commandQueue[[i]]$logFile),
                            "command: %s\nargs: %s\n", commandQueue[[i]]$command,
                            paste0(commandQueue[[i]]$args, collapse = " "))
                    fprintf(file.path(logPath, commandQueue[[i]]$logFile),
                            "\n---\n\noutput:\n%s\n\nstandard error output:\n%s\n",
                            results[[i]]$stdout, results[[i]]$stderr, append = TRUE)
                }, error = function(e) "")
            }
        }
    }
    
    ret <- sapply(results, "[[", "result")
    return(ret)
}
