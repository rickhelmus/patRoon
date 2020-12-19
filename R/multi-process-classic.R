makeCommandList <- function(commandQueue, cmdInds, sucDir)
{
    ncmd <- length(commandQueue)
    ret <- list()
    
    if (ncmd > 1)
    {
        # execute multiple processes at once
        
        cmdList <- sapply(commandQueue,
                          function(cmd) paste(shQuote(cmd$command), paste0(shQuote(cmd$args), collapse = " "),
                                              sep = " "))
        
        if (Sys.info()[["sysname"]] == "Windows")
        {
            # on Windows we easily reach the commandline text limit --> execute as batch file
            
            # fail with last exit code if a command failed: see https://stackoverflow.com/questions/734598/how-do-i-make-a-batch-file-terminate-upon-encountering-an-error
            # NOTE: need to enable delayed expansion and use exclamation marks for current errorlevel expansion.
            ORDoExit <- "|| exit /b !errorlevel!"
            
            # mark success of a command by creating an empty file named after the command index
            ANDMarkSucceed <- paste("&& type NUL >", file.path(sucDir, cmdInds))
            
            # use call in order to be able to execute batch files
            cmdList <- paste("call", cmdList)
            
            shFile <- tempfile(fileext = ".bat")
            cat("SETLOCAL EnableDelayedExpansion",
                paste(cmdList, ANDMarkSucceed, ORDoExit, collapse = "\n"),
                sep = "\n", file = shFile)
            ret$command <- shFile
        }
        else
        {
            # not supported anymore by processx :( --> call sh instead for nix
            # ret$commandline <- paste0(cmdList, collapse = " ; ")
            
            ORDoExit <- "|| exit $?"
            ANDMarkSucceed <- paste("&& touch", file.path(sucDir, cmdInds))
            
            ret$command <- "/bin/sh"
            ret$args <- c("-c", paste0(cmdList, ANDMarkSucceed, ORDoExit, collapse = " && "))
        }
    }
    else
    {
        ret[c("command", "args")] <- commandQueue[[1]][c("command", "args")]
        if (!is.null(commandQueue[[1]]$logFile))
            ret[c("stdout", "stderr")] <- "|"
    }
    
    return(ret)
}

initCommand <- function(commandQueue, cmdInds, sucDir, printOutput, printError)
{
    procArgs <- makeCommandList(commandQueue, cmdInds, sucDir)
    procArgs <- c(procArgs, list(cleanup_tree = TRUE, supervise = TRUE))
    
    if (printOutput)
        procArgs[["stdout"]] <- "|"
    if (printError)
        procArgs[["stderr"]] <- "|"
    
    ncmd <- length(commandQueue)
    
    ret <- list()
    ret$procArgs <- procArgs
    ret$cmdIndRange <- c(cmdInds[1], cmdInds[ncmd])
    ret$finished <- rep(FALSE, ncmd) # only used for timeout checking
    ret$failed <- rep(FALSE, ncmd)
    ret$timedOut <- rep(FALSE, ncmd)
    ret$timeOutRetries <- rep(0, ncmd)
    ret$startTime <- Sys.time()
    ret$noResRetries <- rep(0, ncmd)
    ret$running <- TRUE
    
    return(ret)
}

maybeRestartCommand <- function(commandQueue, procInfo, sucDir, exitStatus, timeoutHandler, errorHandler)
{
    restart <- FALSE
    cmdInds <- seq(procInfo$cmdIndRange[1], procInfo$cmdIndRange[2])
    ncmd <- length(cmdInds)
    
    starti <- tail(which(procInfo$failed), 1)
    if (length(starti) == 0)
        starti <- 0
    
    for (i in seq(starti + 1, ncmd))
    {
        # succeeded?
        if (ncmd > 1 && file.exists(file.path(sucDir, cmdInds[i])))
            next
        
        # already marked as failed during a previous restart?
        if (procInfo$failed[i])
            next
        
        # timed out?
        if (procInfo$timedOut[i])
            retry <- timeoutHandler(cmd = commandQueue[[i]],
                                    retries = procInfo$timeOutRetries[i])
        else
        {
            # we can assume that only the last process in sequence
            # has failed. Hence the exit status must correspond to
            # this process.
            retry <- errorHandler(cmd = commandQueue[[i]], exitStatus = exitStatus,
                                  retries = procInfo$noResRetries[i])
            if (is.character(retry))
                stop(retry) # special case: stop with error message given by the errorHandler
        }
        
        if (retry)
        {
            # Prevents occasional Windows error: "The requested operation cannot be performed on a file with a user-mapped section open"
            Sys.sleep(1)
            
            if (procInfo$timedOut[i])
            {
                procInfo$timeOutRetries[i] <- procInfo$timeOutRetries[i] + 1
                procInfo$timedOut[i] <- FALSE
            }
            else
                procInfo$noResRetries[i] <- procInfo$noResRetries[i] + 1
            
            reStartInd <- i
        }
        else
        {
            procInfo$failed[i] <- TRUE
            reStartInd <- i + 1
        }
        
        if (reStartInd <= ncmd)
        {
            if (reStartInd != 1)
            {
                # adjust command list if not all commands have to be restarted
                pa <- makeCommandList(commandQueue[seq(reStartInd, ncmd)],
                                      seq(cmdInds[reStartInd], cmdInds[ncmd]), sucDir)
                procInfo$procArgs <- modifyList(procInfo$procArgs, pa)
                procInfo$startTime <- Sys.time()
            }
            
            restart <- TRUE
        }
        
        break
    }
    
    return(list(procInfo = procInfo, restart = restart))
}

executeMultiProcessClassic <- function(commandQueue, finishHandler,
                                       timeoutHandler, errorHandler,
                                       prepareHandler, procTimeout,
                                       printOutput, printError,
                                       logSubDir, showProgress, waitTimeout,
                                       batchSize, delayBetweenProc)
{
    if (length(commandQueue) == 0)
        return(list())
    
    maxProcAmount <- getOption("patRoon.MP.maxProcs")
    runningProcs <- vector("list", maxProcAmount)
    runningProcInfo <- vector("list", maxProcAmount)
    
    if (!is.null(prepareHandler))
        commandQueue <- lapply(commandQueue, prepareHandler)
    
    totCmdCount <- length(commandQueue)
    
    ret <- vector("list", totCmdCount)
    names(ret) <- names(commandQueue)
    
    if (showProgress)
        prog <- openProgBar(0, totCmdCount)
    
    nextCommand <- 1
    finishedCommands <- 0
    lastCommandTime <- 0 # at which time (in ms) the last command was started

    logPath <- getOption("patRoon.MP.logPath", FALSE)
    if (!is.null(logSubDir) && !isFALSE(logPath))
    {
        logPath <- file.path(logPath, logSubDir)
        mkdirp(logPath)
        commandQueue <- lapply(commandQueue, function(cmd) { cmd$logFile <- file.path(logPath, cmd$logFile); return(cmd) })
    }
    else
        logPath <- NULL
    doLog <- !is.null(logPath)
    
    stopifnot(batchSize == 1 || (!doLog && !printOutput && !printError))
    
    # clear up stale processes: see https://github.com/r-lib/processx/issues/171
    on.exit({
        for (pi in seq_along(runningProcs))
        {
            if (!is.null(runningProcs[[pi]]) && runningProcInfo[[pi]]$running)
                runningProcs[[pi]]$kill()
        }
        
        if (doLog)
        {
            for (cmd in commandQueue)
            {
                if (!is.null(cmd$logFile))
                {
                    tryCatch({
                        fprintf(cmd$logFile, "command: %s\nargs: %s\n", cmd$command, paste0(cmd$args, collapse = " "))
                        fprintf(cmd$logFile, "\n---\n\noutput:\n%s\n\nstandard error output:\n%s\n",
                                cmd$stdoutLog, cmd$stderrLog, append = TRUE)
                    }, error = function(e) "")
                }
            }
        }
    }, add = TRUE)
    
    sucDir <- tempfile("suc")
    dir.create(sucDir)
    
    procFinished <- function(pi) !is.null(runningProcs[[pi]]) && !runningProcs[[pi]]$is_alive() && runningProcInfo[[pi]]$running
    
    doProcessOut <- function(txt, print)
    {
        if (print)
            cat(txt)
        if (doLog)
            return(txt)
    }
    
    # reading process output might fail sometimes(?)
    emptyStrOnErr <- function(expr) tryCatch(expr, error = function(e) "")
    
    while (nextCommand <= totCmdCount || any(sapply(runningProcInfo, function(rp) !is.null(rp) && rp$running)))
    {
        for (pi in seq_along(runningProcs))
        {
            finishedRunning <- procFinished(pi)
            
            if (!is.null(runningProcs[[pi]]))
            {
                cmdInds <- seq(runningProcInfo[[pi]]$cmdIndRange[1], runningProcInfo[[pi]]$cmdIndRange[2])
                
                # NOTE: logging/printing currently doesn't work in batch mode
                if (printOutput || printError || doLog)
                {
                    cind <- runningProcInfo[[pi]]$cmdIndRange[1]
                    commandQueue[[cind]]$stdoutLog <- paste0(commandQueue[[cind]]$stdoutLog,
                                                             emptyStrOnErr(doProcessOut(rp[[pi]]$read_output(), printOutput)))
                    commandQueue[[cind]]$stderrLog <- paste0(commandQueue[[cind]]$stderrLog,
                                                             emptyStrOnErr(doProcessOut(rp[[pi]]$read_error(), printError)))
                }
            }
            
            if (finishedRunning)
            {
                ncmd <- length(cmdInds)
                
                # NOTE: as per docs get_exit_status() might return NA, in this
                # case check if a command failed (by checking for missing
                # success marker files)
                # UNDONE: when batchSize=1 we don't create/check success
                # markers. Fix this? So far never had NA exit statuses in that
                # situation.
                exitStatus <- runningProcs[[pi]]$get_exit_status()
                
                if (is.na(exitStatus) || exitStatus != 0) # something (may have) failed?
                {
                    maybe <- maybeRestartCommand(commandQueue[cmdInds], runningProcInfo[[pi]], sucDir, exitStatus,
                                                 timeoutHandler, errorHandler)
                    runningProcInfo[[pi]] <- maybe$procInfo # might have been updated
                    
                    if (maybe$restart)
                    {
                        runningProcs[[pi]] <- do.call(processx::process$new, runningProcInfo[[pi]]$procArgs)
                        finishedRunning <- FALSE
                    }
                }
            }
            
            if (is.null(runningProcs[[pi]]) || finishedRunning)
            {
                finishedProc <- runningProcs[[pi]]
                finishedProcInfo <- runningProcInfo[[pi]]
                
                # start next command first, so we can process results while the next is already running
                if (nextCommand <= totCmdCount)
                {
                    if (delayBetweenProc > 0)
                    {
                        diffTime <- curTimeMS() - lastCommandTime
                        if (diffTime < delayBetweenProc)
                            Sys.sleep((delayBetweenProc - diffTime) / 1000)
                    }
                    
                    if (batchSize == 1)
                        ncmd <- 1
                    else
                    {
                        cmdLeft <- totCmdCount - (nextCommand - 1)
                        freeSlots <- sum(sapply(seq_along(runningProcs), function(i) is.null(runningProcs[[i]]) || procFinished(i)))
                        freeSize <- freeSlots * batchSize
                        if (freeSlots > 1 && cmdLeft < freeSize)
                            ncmd <- ceiling(cmdLeft / freeSlots) # divide between free slots
                        else
                            ncmd <- min(batchSize, cmdLeft)
                        # printf("ncmd %d (slot %d, remain: %d, free: %d)\n", ncmd, pi, cmdLeft, freeSlots)
                    }
                    
                    cs <- seq(nextCommand, nextCommand + (ncmd - 1))
                    runningProcInfo[[pi]] <- initCommand(commandQueue[cs], cs, sucDir, printOutput, printError)
                    runningProcs[[pi]] <- do.call(processx::process$new, runningProcInfo[[pi]]$procArgs)
                    
                    # printf("started %d-%d on slot %d\n", nextCommand, runningProcInfo[[pi]]$cmdIndRange[2], pi)
                    lastCommandTime <- curTimeMS()
                    
                    if (!is.null(finishedProc))
                    {
                        # printf("prev time: %s\n", (Sys.time() - finishedProc$get_start_time()) * 1000)
                        finishedCommands <- finishedCommands +
                            (runningProcInfo[[pi]]$cmdIndRange[2]+1) - runningProcInfo[[pi]]$cmdIndRange[1]
                        if (showProgress)
                            setTxtProgressBar(prog, finishedCommands)
                    }
                    
                    nextCommand <- nextCommand + ncmd
                }
                else if (!is.null(runningProcs[[pi]]))
                    runningProcInfo[[pi]]$running <- FALSE
                
                if (!is.null(finishedProc))
                {
                    inds <- cmdInds[!finishedProcInfo$failed]
                    if (length(inds) > 0)
                        ret[inds] <- lapply(inds, function(ci) finishHandler(cmd = commandQueue[[ci]]))
                }
                
            }
            else if (!is.null(procTimeout) && runningProcInfo[[pi]]$running)
            {
                # check for timeouts
                kill <- FALSE
                
                if (length(cmdInds) > 1)
                {
                    # for batch execution: update start time if a new command was started
                    
                    for (i in seq_along(cmdInds))
                    {
                        if (runningProcInfo[[pi]]$failed[i] || runningProcInfo[[pi]]$finished[i])
                            next
                        
                        if (file.exists(file.path(sucDir, cmdInds[i])))
                        {
                            runningProcInfo[[pi]]$finished[i] <- TRUE # now finished
                            runningProcInfo[[pi]]$startTime <- Sys.time()
                        }
                        else if (difftime(Sys.time(), runningProcInfo[[pi]]$startTime, units = "secs")[[1]] > procTimeout)
                            runningProcInfo[[pi]]$timedOut[i] <- kill <- TRUE
                        
                        break
                    }
                }
                else
                    runningProcInfo[[pi]]$timedOut[1] <- kill <-
                        difftime(Sys.time(), runningProcInfo[[pi]]$startTime, units = "secs") > procTimeout
                
                if (kill)
                {
                    runningProcs[[pi]]$kill()
                    runningProcs[[pi]]$wait()
                }
            }
        }
        
        if (printOutput || printError || doLog)
        {
            rp <- pruneList(runningProcs)
            pl <- processx::poll(rp, waitTimeout)
            
            if (FALSE)
            {
                for (pi in seq_along(rp))
                {
                    if (doLog) # NOTE: logging currently doesn't work in batch mode
                        cind <- runningProcInfo[[pi]]$cmdIndRange[1]
                    
                    if (pl[[pi]][["output"]] == "ready")
                    {
                        txt <- rp[[pi]]$read_output()
                        
                        if (printOutput)
                            cat(txt)
                        if (doLog)
                            commandQueue[[cind]]$stdoutLog <- paste0(commandQueue[[cind]]$stdoutLog, txt)
                        
                    }
                    if (pl[[pi]][["error"]] == "ready")
                    {
                        txt <- rp[[pi]]$read_error()
                        if (printError)
                            cat(txt)
                        if (doLog)
                            commandQueue[[cind]]$stderrLog <- paste0(commandQueue[[cind]]$stderrLog, txt)
                    }
                }
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
                    break
                }
            }
        }
    }
    
    # get rid of potentially large amount of temporary files
    unlink(sucDir, recursive = TRUE)
    
    if (showProgress)
    {
        setTxtProgressBar(prog, totCmdCount)
        close(prog)
    }
    
    return(ret)
}
