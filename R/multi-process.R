defMultiProcErrorHandler <- function(cmd, exitStatus, ...)
{
    stop(sprintf("Failed to run command '%s' with args: '%s'. Exit code: %d",
                 cmd$command, paste0(cmd$args, collapse = " "), exitStatus))
}

executeMultiProcess <- function(...)
{
    # wrap the actual code here so we can act on user interrupts: calling gc()
    # afterwards makes sure processes are cleaned up
    tryCatch(executeMultiProcess2(...), interrupt = function(e)
    {
        invisible(capture.output(gc(verbose = FALSE)))
        # re-throw: Thanks to @gaborcsardi (https://github.com/r-lib/processx/issues/171)
        signalCondition(e)
        invokeRestart("abort")
    })
}

executeMultiProcess2 <- function(commandQueue, finishHandler,
                                timeoutHandler = function(...) TRUE,
                                errorHandler = defMultiProcErrorHandler,
                                procTimeout = NULL, printOutput = FALSE, printError = FALSE,
                                showProgress = TRUE, progressOut = "", waitTimeout = 50,
                                maxProcAmount = getOption("patRoon.maxProcAmount"),
                                maxCmdsPerProc = 1, delayBetweenProc = 0)
{
    if (length(commandQueue) == 0)
        return(list())

    runningProcs <- vector("list", maxProcAmount)
    runningProcInfo <- vector("list", maxProcAmount)
    totCmdCount <- length(commandQueue)

    ret <- vector("list", totCmdCount)
    names(ret) <- names(commandQueue)

    if (showProgress)
        prog <- txtProgressBar(0, totCmdCount, style = 3, file = progressOut)

    nextCommand <- 1
    finishedCommands <- 0
    lastCommandTime <- 0 # at which time (in ms) the last command was started

    while (nextCommand <= totCmdCount || any(sapply(runningProcInfo, function(rp) !is.null(rp) && rp$running)))
    {
        for (pi in seq_along(runningProcs))
        {
            finishedRunning <- !is.null(runningProcs[[pi]]) && !runningProcs[[pi]]$is_alive() && runningProcInfo[[pi]]$running

            if (!is.null(runningProcs[[pi]]))
                cmdInds <- seq(runningProcInfo[[pi]]$startCmdInd, runningProcInfo[[pi]]$endCmdInd)
            
            if (finishedRunning && !runningProcInfo[[pi]]$failed)
            {
                exitStatus = runningProcs[[pi]]$get_exit_status()
                if (exitStatus != 0)
                {
                    if (all(sapply(cmdInds, function(ci) errorHandler(cmd = commandQueue[[ci]],
                                                                      exitStatus = exitStatus,
                                                                      retries = runningProcInfo[[pi]]$noResRetries))))
                    {
                        # Prevents occasional error: "The requested operation cannot be performed on a file with a user-mapped section open"
                        Sys.sleep(1)
                        
                        runningProcInfo[[pi]]$noResRetries <- runningProcInfo[[pi]]$noResRetries + 1
                        runningProcs[[pi]] <- do.call(process$new, runningProcInfo[[pi]]$procArgs)
                        finishedRunning <- FALSE
                    }
                    else
                        runningProcInfo[[pi]]$failed <- TRUE
                }
            }
            
            if (is.null(runningProcs[[pi]]) || finishedRunning)
            {
                finishedProc <- runningProcs[[pi]]
                finishedProcInfo <- runningProcInfo[[pi]]

                # start next command first, so we can process results while next is already running
                if (nextCommand <= totCmdCount)
                {
                    if (delayBetweenProc > 0)
                    {
                        diffTime <- curTimeMS() - lastCommandTime
                        if (diffTime < delayBetweenProc)
                            Sys.sleep((delayBetweenProc - diffTime) / 1000)
                    }
                    
                    nproc <- min(maxCmdsPerProc, (totCmdCount - (nextCommand-1)))

                    procArgs <- list(cleanup_tree = TRUE, supervise = TRUE)
                    if (nproc > 1)
                    {
                        # execute multiple processes at once

                        cmdList <- sapply(commandQueue[seq(nextCommand, nextCommand + (nproc - 1))],
                                          function(cmd) paste(shQuote(cmd$command), paste0(shQuote(cmd$args), collapse = " "),
                                                              sep = " "))

                        if (Sys.info()[["sysname"]] == "Windows")
                        {
                            # on Windows we easily reach the commandline text limit --> execute as batch file
                            shFile <- tempfile(fileext = ".bat")
                            write(paste0(cmdList, collapse = "\n"), shFile)
                            procArgs$command <- shFile
                        }
                        else
                        {
                            # not supported anymore by processx :(
                            # procArgs$commandline <- paste0(cmdList, collapse = " ; ")
                            
                            # --> call sh instead for nix
                            procArgs$command <- "/bin/sh"
                            procArgs$args <- c("-c", paste0(cmdList, collapse = " ; "))
                        }
                    }
                    else
                    {
                        procArgs[c("command", "args")] <- commandQueue[[nextCommand]][c("command", "args")]
                        if (!is.null(commandQueue[[nextCommand]]$stdoutFile))
                            procArgs$stdout <- commandQueue[[nextCommand]]$stdoutFile
                        if (!is.null(commandQueue[[nextCommand]]$stderrFile))
                            procArgs$stderr <- commandQueue[[nextCommand]]$stderrFile
                    }

                    if (printOutput)
                        procArgs <- c(procArgs, stdout = "|")
                    if (printError)
                        procArgs <- c(procArgs, stderr = "|")

                    runningProcs[[pi]] <- do.call(process$new, procArgs)
                    runningProcInfo[[pi]]$procArgs <- procArgs
                    runningProcInfo[[pi]]$startCmdInd <- nextCommand
                    runningProcInfo[[pi]]$endCmdInd <- nextCommand + (nproc - 1)
                    runningProcInfo[[pi]]$failed <- FALSE
                    runningProcInfo[[pi]]$timeOutRetries <- 0
                    runningProcInfo[[pi]]$noResRetries <- 0
                    runningProcInfo[[pi]]$running <- TRUE

                    # printf("started %d-%d on slot %d\n", nextCommand, runningProcInfo[[pi]]$endCmdInd, pi)
                    lastCommandTime <- curTimeMS()

                    if (!is.null(finishedProc))
                    {
                        # printf("prev time: %s\n", (Sys.time() - finishedProc$get_start_time()) * 1000)
                        finishedCommands <- finishedCommands +
                            (runningProcInfo[[pi]]$endCmdInd+1) - runningProcInfo[[pi]]$startCmdInd
                        if (showProgress)
                            setTxtProgressBar(prog, finishedCommands)
                    }

                    nextCommand <- nextCommand + nproc
                }
                else if (!is.null(runningProcs[[pi]]))
                    runningProcInfo[[pi]]$running <- FALSE

                if (!is.null(finishedProc) && !finishedProcInfo$failed)
                    ret[cmdInds] <- lapply(cmdInds, function(ci) finishHandler(cmd = commandQueue[[ci]]))
            }
            else if (!is.null(procTimeout) && runningProcInfo[[pi]]$running &&
                     (Sys.time() - runningProcs[[pi]]$get_start_time()) > procTimeout)
            {
                runningProcs[[pi]]$kill()
                runningProcs[[pi]]$wait()

                if (all(sapply(cmdInds, function(ci) timeoutHandler(cmd = commandQueue[[ci]],
                                                                    retries = runningProcInfo[[pi]]$timeOutRetries))))
                {
                    # As with restarting failed commands short sleep may prevent errors (this situation was never tested)
                    Sys.sleep(1)
                    
                    runningProcInfo[[pi]]$timeOutRetries <- runningProcInfo[[pi]]$timeOutRetries + 1
                    runningProcs[[pi]] <- do.call(process$new, runningProcInfo[[pi]]$procArgs)
                }
                else
                    runningProcInfo[[pi]]$failed <- TRUE
            }
        }

        if (printOutput || printError)
        {
            rp <- pruneList(runningProcs)
            pl <- poll(rp, waitTimeout)

            for (pi in seq_along(rp))
            {
                if (pl[[pi]][["output"]] == "ready")
                    cat(rp[[pi]]$read_output_lines())
                if (pl[[pi]][["error"]] == "ready")
                    cat(rp[[pi]]$read_error_lines())
            }
        }
        else
        {
            # just wait for one of the running processes
            for (pi in seq_along(runningProcs))
            {
                if (!is.null(runningProcInfo[[pi]]) && runningProcInfo[[pi]]$running)
                {
                    runningProcs[[pi]]$wait(waitTimeout)
                    # Sys.sleep(waitTimeout / 1000)
                    break
                }
            }
        }
    }

    if (showProgress)
    {
        setTxtProgressBar(prog, totCmdCount)
        close(prog)
    }

    return(ret)
}
