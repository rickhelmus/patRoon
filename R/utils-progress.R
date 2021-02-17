# This function is heavily inspired by the excellent progressr package by Henrik Bengtsson:
# https://github.com/HenrikBengtsson/progressr
# Some day we might switch (back) to this package, but since it requires an opt-in by the user, we stick with this
# dumbed down version for now...
withProg <- function(end, expr)
{
    prog <- 0
    pb <- openProgBar(0, end, file = stderr())

    # from https://stackoverflow.com/questions/56038299/in-r-how-do-i-evaluate-an-expression-in-a-specific-environment-within-a-functio
    ret <- withCallingHandlers(withVisible(expr), doProg = function(...)
    {
        if (prog < end)
        {
            prog <<- prog + 1
            setTxtProgressBar(pb, prog)
        }
    })

    setTxtProgressBar(pb, end)
    close(pb)
    
    if (ret$visible)
        return(ret$value)
    invisible(ret$value)
}

# helper function in the expression environment to signal progression
doProgress <- function() signalCondition(structure(class = c("doProg", "immediateCondition", "condition"),
                                                   list(message = character(), call = sys.call())))
