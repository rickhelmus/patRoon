executeMultiProcessFuture <- function(commandQueue, finishHandler,
                                      timeoutHandler = function(...) TRUE,
                                      errorHandler = defMultiProcErrorHandler,
                                      prepareHandler = NULL,
                                      procTimeout = NULL, printOutput = FALSE, printError = FALSE,
                                      showProgress = TRUE, waitTimeout = 50,
                                      batchSize = 1, delayBetweenProc = 0)
{
    ret <- future.apply::future_lapply(commandQueue, function(cmd)
    {
        if (!is.null(prepareHandler))
            cmd <- prepareHandler(cmd) # UNDONE
        timeoutRetries <- errorRetries <- 0
        while (TRUE)
        {
            stat <- processx::run(cmd$command, cmd$args, error_on_status = FALSE,
                                  timeout = procTimeout, cleanup_tree = TRUE)
            
            if (!is.null(cmd$logFile))
            {
                tryCatch({
                    fprintf(cmd$logFile, "command: %s\nargs: %s\n", cmd$command, paste0(cmd$args, collapse = " "))
                    fprintf(cmd$logFile, "\n---\n\noutput:\n%s\n\nstandard error output:\n%s\n",
                            stat$stdout, stat$stderr, append = TRUE)
                }, error = function(e) "")
            }
            
            if (stat$timeout)
            {
                if (!timeoutHandler(cmd = cmd, retries = errorRetries))
                    break
                timeoutRetries <- timeoutRetries + 1
            }
            else if (stat$status != 0)
            {
                if (!errorHandler(cmd = cmd, exitStatus = stat$status, retries = errorRetries))
                    break
                errorRetries <- errorRetries + 1
            }
            else # success
                return(finishHandler(cmd))
            
        }
    }, future.scheduling = 1.0)
    return(ret)
}
