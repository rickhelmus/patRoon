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
            # NOTE: need to enable delayed expension and use exclamation marks for current errorlevel expansion.
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

defMultiProcErrorHandler <- function(cmd, exitStatus, ...)
{
    stop(sprintf("Failed to run command '%s' with args: '%s'. Exit code: %d",
                 cmd$command, paste0(cmd$args, collapse = " "), exitStatus))
}

doProcFuture <- function(commandQueue, finishHandler, timeoutHandler, errorHandler,
                         prepareHandler, procTimeout, printOutput, printError, waitTimeout, delayBetweenProc)
{
    # UNDONE: printOutput/printError not here?
    
    {
        sucDir <- withr::local_tempfile()
        dir.create(sucDir)
        
        if (!is.null(prepareHandler))
            commandQueue <- prepareHandler(commandQueue)
        
        lastCommandTime <- 0 # at which time (in ms) the last command was started
        ret <- vector("list", length(commandQueue))

        # UNDONE: or just always wait here? since lastCommandTime is only available in this future
        # if (delayBetweenProc > 0)
        # {
        #     diffTime <- curTimeMS() - lastCommandTime
        #     if (diffTime < delayBetweenProc)
        #         Sys.sleep((delayBetweenProc - diffTime) / 1000)
        # }
            
        ncmd <- length(commandQueue)
        procInfo <- patRoon:::initCommand(commandQueue, seq_along(commandQueue), sucDir, printOutput, printError)
        proc <- do.call(processx::process$new, procInfo$procArgs)
        # lastCommandTime <- curTimeMS()
        
        while (TRUE)
        {
            if (!is.null(proc) && !proc$is_alive() && procInfo$running)
            {
                # NOTE: as per docs get_exit_status() might return NA, in this
                # case check if a command failed (by checking for missing
                # success marker files)
                # UNDONE: when batchSize=1 we don't create/check success
                # markers. Fix this? So far never had NA exit statuses in that
                # situation.
                exitStatus <- proc$get_exit_status()
                
                if (is.na(exitStatus) || exitStatus != 0) # something (may have) failed?
                {
                    maybe <- maybeRestartCommand(commandQueue, procInfo, sucDir, exitStatus,
                                                 timeoutHandler, errorHandler)
                    procInfo <- maybe$procInfo # might have been updated
                    
                    if (maybe$restart)
                        proc <- do.call(processx::process$new, procInfo$procArgs)
                    else
                        break
                }
                else
                    break
            }
            else if (!is.null(procTimeout))
            {
                # check for timeouts
                kill <- FALSE
                
                if (length(commandQueue) > 1)
                {
                    # for batch execution: update start time if a new command was started
                    for (i in seq_along(commandQueue))
                    {
                        if (procInfo$failed[i] || procInfo$finished[i])
                            next
                        
                        if (file.exists(file.path(sucDir, i)))
                        {
                            procInfo$finished[i] <- TRUE # now finished
                            procInfo$startTime <- Sys.time()
                        }
                        else if (difftime(Sys.time(), procInfo$startTime, units = "secs")[[1]] > procTimeout)
                            procInfo$timedOut[i] <- kill <- TRUE
                        
                        break
                    }
                }
                else
                    procInfo$timedOut[1] <- kill <- difftime(Sys.time(), procInfo$startTime, units = "secs") > procTimeout
                
                if (kill)
                {
                    proc$kill()
                    proc$wait()
                }
            }
            else
                proc$wait(waitTimeout)
        }
        
        inds <- which(!procInfo$failed)
        if (length(inds) > 0)
            ret[inds] <- lapply(inds, function(ci) finishHandler(cmd = commandQueue[[ci]]))
        
        ret
    }
}

createProcFuture <- function(expr, globals)
{
    expr <- createProcFutureExpr(commandQueue, finishHandler, timeoutHandler, errorHandler,
                                 procTimeout, printOutput, printError, waitTimeout, delayBetweenProc)
    future::future(eval(expr), globals = globals)
}

executeMultiProcessF <- function(commandQueue, finishHandler,
                                 timeoutHandler = function(...) TRUE,
                                 errorHandler = defMultiProcErrorHandler,
                                 prepareHandler = NULL,
                                 procTimeout = NULL, printOutput = FALSE, printError = FALSE,
                                 showProgress = TRUE, waitTimeout = 50,
                                 maxProcAmount = NULL,
                                 batchSize = 1, delayBetweenProc = 0)
{
    if (FALSE)
    {
        expr <- quote(patRoon:::doProcFuture(commandQueue, finishHandler, timeoutHandler, errorHandler,
                                             procTimeout, printOutput, printError, waitTimeout, delayBetweenProc))
        
        globals <- future::getGlobalsAndPackages(expr)$globals
        futures <- lapply(commandQueue, function(cmd)
        {
            g <- modifyList(globals, list(commandQueue = list(cmd)))
            future::future({ patRoon:::doProcFuture(list(cmd), finishHandler, timeoutHandler, errorHandler,
                                                    procTimeout, printOutput, printError, waitTimeout, delayBetweenProc) }, globals = g)
        })
        # future::resolve(futures, result = TRUE)
        return(unlist(future::values(futures), recursive = FALSE))
    }
    else if (F)
    {
        # executeMultiProcess(commandQueue = commandQueue,
        #                       finishHandler = finishHandler, timeoutHandler = timeoutHandler,
        #                       errorHandler = errorHandler, prepareHandler = prepareHandler, procTimeout = procTimeout,
        #                       printOutput = printOutput, printError = printError, showProgress = showProgress,
        #                       waitTimeout = waitTimeout, maxProcAmount = 4, batchSize = batchSize,
        #                       delayBetweenProc = delayBetweenProc)
        
        n <- future::nbrOfWorkers() # UNDONE: configurable? check for Inf (docs state it may be, then default to 1?)
        if (n == 1)
            chunks <- list(commandQueue)
        else
            chunks <- split(commandQueue, cut(seq_along(commandQueue), n, labels = FALSE))
        args <- list(finishHandler = finishHandler, timeoutHandler = timeoutHandler,
                     errorHandler = errorHandler, prepareHandler = prepareHandler, procTimeout = procTimeout,
                     printOutput = printOutput, printError = printError, showProgress = showProgress,
                     waitTimeout = waitTimeout, batchSize = batchSize,
                     delayBetweenProc = delayBetweenProc)
        if (!is.null(maxProcAmount))
            args <- c(args, maxProcAmount = maxProcAmount) # UNDONE?
        ret <- do.call(future.apply::future_lapply, c(list(chunks, executeMultiProcess), args))
        return(unlist(unname(ret), recursive = FALSE))
    }
    else if (FALSE)
    {
        ret <- future.apply::future_lapply(commandQueue, function(cmd)
        {
            doProcFuture(list(cmd), finishHandler, timeoutHandler, errorHandler, prepareHandler,
                         procTimeout, printOutput, printError, waitTimeout, delayBetweenProc)
        }, future.scheduling = 4.0)
        return(unlist(unname(ret), recursive = FALSE))
    }
    else
    {
        ret <- future.apply::future_lapply(commandQueue, function(cmd)
        {
            if (!is.null(prepareHandler))
                cmd <- prepareHandler(list(cmd))[[1]] # UNDONE
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
}

#' Simultaneous execution of system commands.
#'
#' Execute a queue of system commands in parallel.
#'
#' This function executes a given queue with system commands in parallel to
#' speed up computation. Commands are executed in the background using the
#' \pkg{processx} package. A configurable maximum amount of processes are
#' created to execute multiple commands in parallel.
#'
#' Multiple commands may be executed in sequence that are launched from a single
#' parent process (as part of a batch script on Windows or combined with the
#' shell AND operator otherwise). Note that in this scenario still multiple
#' processes are spawned. Each of these processes will manage a chunk of the
#' command queue (size defined by \code{batchSize} argument). This approach is
#' typically suitable for fast running commands: the overhead of spawning a new
#' process for each command from R would in this case be significant enough to
#' loose most of the speedup otherwise gained with parallel execution. Note that
#' the actual batch size may be adjusted to ensure that a maximum number of
#' processes are running simultaneously.
#'
#' Other functionalities of this function include timeout and error handling.
#'
#' @param commandQueue A list with commands. Should contain \code{command}
#'   (scalar string) and \code{args} (\code{character} vector). More user
#'   defineds fields are allowed and useful to attach command information that
#'   can be used in the finish, timeout and error handlers.
#' @param finishHandler A function that is called when a command has finished.
#'   This function is typically used to process any results generated by the
#'   command. The function is called right after spawning a new process, hence
#'   processing results can occur while the next command is running in the
#'   background. The function signature should be \code{function(cmd)} where
#'   \code{cmd} is the queue data (from \code{commandQueue}) of the command that
#'   has finished.
#' @param timeoutHandler A function that is called whenever a timeout for a
#'   command occurs. Should return \code{TRUE} if execution of the command
#'   should be retried. The function signature should be \code{function(cmd,
#'   retries)} where \code{cmd} is the queue data for that command and
#'   \code{retries} the number of times the command has been retried.
#' @param errorHandler Similar to \code{timeoutHandler}, but called whenever a
#'   command has failed. The signature should be \code{function(cmd, exitStatus,
#'   retries)}. The \code{exitStatus} argument is the exit code of the command
#'   (may be \code{NA} in rare cases this is unknown). Other arguments are as
#'   \code{timeoutHandler}.
#' @param procTimeout The maximum time a process may consume before a timeout
#'   occurs (in seconds). Set to \code{NULL} to disable
#'   timeouts.
#' @param printOutput,printError Set to \code{TRUE} to print stdout/stderr
#'   output to the console. Currently unused and untested.
#' @param showProgress Set to \code{TRUE} to display a progress bar.
#' @param waitTimeout Number of milliseconds to wait before checking if a new
#'   process should be spawned.
#' @param maxProcAmount Maximum number of processes to run simultaneously.
#' @param batchSize Number of commands that should be executed in sequence per
#'   processes. See details.
#' @param delayBetweenProc Minimum number of milliseconds to wait before
#'   spawning a new process. Might be needed to workaround errors.
#'
#' @keywords internal
executeMultiProcess <- function(commandQueue, finishHandler,
                                timeoutHandler = function(...) TRUE,
                                errorHandler = defMultiProcErrorHandler,
                                prepareHandler = NULL,
                                procTimeout = NULL, printOutput = FALSE, printError = FALSE,
                                showProgress = TRUE, waitTimeout = 50,
                                maxProcAmount = getOption("patRoon.maxProcAmount"),
                                batchSize = 1, delayBetweenProc = 0)
{
    if (length(commandQueue) == 0)
        return(list())

    runningProcs <- vector("list", maxProcAmount)
    runningProcInfo <- vector("list", maxProcAmount)
    
    if (!is.null(prepareHandler))
        commandQueue <- prepareHandler(commandQueue)
    
    totCmdCount <- length(commandQueue)

    ret <- vector("list", totCmdCount)
    names(ret) <- names(commandQueue)

    if (showProgress)
        prog <- openProgBar(0, totCmdCount)

    nextCommand <- 1
    finishedCommands <- 0
    lastCommandTime <- 0 # at which time (in ms) the last command was started

    doLog <- any(sapply(commandQueue, function(q) !is.null(q$logFile)))
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
