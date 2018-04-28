executeMultiProcess <- function(commandQueue, finishProc, printOutput = FALSE,
                                printError = FALSE, showProgress = TRUE, progressOut = "", waitTimeout = 50,
                                procTimeout = NULL, timeoutRetries = 3, maxProcAmount = getOption("patRoon.maxProcAmount"),
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
            if (is.null(runningProcs[[pi]]) || (!runningProcs[[pi]]$is_alive() && runningProcInfo[[pi]]$running))
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

                    procArgs <- list()
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
                            procArgs$commandline <- paste0(cmdList, collapse = " ; ")
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
                    runningProcInfo[[pi]]$startCmdInd <- nextCommand
                    runningProcInfo[[pi]]$endCmdInd <- nextCommand + (nproc - 1)
                    runningProcInfo[[pi]]$killed <- FALSE
                    runningProcInfo[[pi]]$retries <- timeoutRetries
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

                if (!is.null(finishedProc))
                {
                    inds <- seq(finishedProcInfo$startCmdInd, finishedProcInfo$endCmdInd)
                    ret[inds] <- lapply(inds, function(ci) finishProc(commandQueue[[ci]],
                                                                      finishedProc$get_exit_status(),
                                                                      finishedProcInfo$retries))
                }
            }
            else if (!is.null(procTimeout) && runningProcInfo[[pi]]$running &&
                     (Sys.time() - runningProcs[[pi]]$get_start_time()) > procTimeout)
            {
                runningProcs[[pi]]$kill()
                runningProcs[[pi]]$wait()

                if (runningProcInfo[[pi]]$retries > 0)
                {
                    runningProcInfo[[pi]]$retries <- runningProcInfo[[pi]]$retries - 1
                    runningProcs[[pi]] <- runningProcs[[pi]]$restart()
                }
                else
                    runningProcInfo[[pi]]$killed <- TRUE
            }
        }

        if (printOutput || printError)
        {
            pl <- poll(runningProcs, waitTimeout)

            for (pi in seq_along(runningProcs))
            {
                if (pl[[pi]]$output == "ready")
                    cat(runningProcs[[pi]]$read_output_lines())
                if (pl[[pi]]$error == "ready")
                    cat(runningProcs[[ci]]$read_error_lines())
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
