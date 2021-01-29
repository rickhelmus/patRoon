#' @include utils-IPO.R
#' @include main.R
NULL

DoEOptimizer <- setRefClass("DoEOptimizer", contains = "VIRTUAL",
                            fields = list(algorithm = "character", maxModelDeviation = "numeric",
                                          parallel = "logical"))

DoEOptimizer$methods(

    # dummy methods that may need to be overloaded
    checkInitialParams = function(params) params,
    defaultParamRanges = function(params) list(),
    convertOptToCallParams = function(params) params,
    fixDesignParam = function(param, value) value,
    fixOptParamBounds = function(params, bounds) bounds,
    fixOptParams = function(params) params,

    # "virtual" methods
    resultIncreased = function(history) stop("VIRTUAL"),
    calculateResponse = function(params, task, keepObject) stop("VIRTUAL"),
    getResponseScores = function(response) stop("VIRTUAL"),
    getFinalScore = function(response) stop("VIRTUAL"),

    getOptSettingRange = function(settingName, params, paramRanges)
    {
        if (!is.null(paramRanges[[settingName]]))
            return(paramRanges[[settingName]])
        return(c(1, Inf)) # UNDONE: for grouping IPO has default of zero... change?
    },

    # Adapted from combineParams() function of IPO
    combineParams = function(params_1, params_2)
    {
        len <- max(unlist(sapply(params_1, length)))

        p_names <- c(names(params_1), names(params_2))
        for (i in seq_along(params_2))
        {
            new_index <- length(params_1) + 1
            fact <- params_2[[i]]
            params_1[[new_index]] <- fact
            params_1[[new_index]][1:len] <- fact
        }
        names(params_1) <- p_names
        return(params_1)
    },

    # Heavily based on xcmsSetExperimentsCluster() from IPO
    performIteration = function(params)
    {
        typParams <- utilsIPO$typeCastParams(params)

        if (length(typParams$to_optimize) > 1)
        {
            design <- utilsIPO$getCcdParameter(typParams$to_optimize)
            designParams <- rsm::decode.data(design)
        }
        else if (length(typParams$to_optimize) == 0)
            stop("No parameters specified for optimization!")
        else
        {
            design <- data.frame(run.order = 1:9, a = seq(-1,1,0.25))
            colnames(design)[2] <- names(typParams$to_optimize)
            designParams <- design
            designParams[,2] <- seq(min(typParams$to_optimize[[1]]),
                                    max(typParams$to_optimize[[1]]),
                                    diff(typParams$to_optimize[[1]]) / 8)
        }
        
        for (col in names(typParams$to_optimize))
            designParams[[col]] <- sapply(designParams[[col]], .self$fixDesignParam, param = col)

        printf("---\nDesign:\n")
        print(designParams)
        printf("---\n")

        designParams <- combineParams(designParams, typParams$no_optimization)
        tasks <- seq_len(nrow(design))

        prog <- progressr::progressor(steps = length(tasks))

        doExp <- function(task)
        {
            # simplified optimizeSlaveCluster() from IPO

            runParams <- as.list(designParams[task, ])
            runParams <- runParams[!names(runParams) %in% c("run.order", "std.order", "Block")]
            runParams <- convertOptToCallParams(runParams)
            result <- calculateResponse(runParams, task, FALSE)
            result$experiment <- task

            prog()

            return(result)
        }
        
        if (parallel)
            response <- rbindlist(future.apply::future_lapply(tasks, doExp))
        else
            response <- rbindlist(lapply(tasks, doExp))

        ret <- list()
        ret$params <- typParams
        ret$design <- design
        ret$response <- response

        return(ret)
    },

    # Heavily based on xcmsSetStatistic() from IPO
    performIterationStat = function(result, maxModelDeviation)
    {
        params <- result$params
        resp <- getResponseScores(result$response)
        result$response$score <- resp

        result$model <- utilsIPO$createModel(result$design, params$to_optimize, resp)
        result$max_settings <- utilsIPO$getMaximumLevels(result$model)

        runParams <- as.list(utilsIPO$decodeAll(result$max_settings[-1], params$to_optimize))
        runParams <- combineParams(runParams, params$no_optimization)

        if (!is.list(runParams))
            runParams <- as.list(runParams)

        runParams <- Map(.self$fixDesignParam, names(runParams), runParams)
        runParams <- convertOptToCallParams(runParams)
        result$finalResult <- calculateResponse(runParams, 1, TRUE)
        result$finalResult$parameters <- runParams
        result$finalResult$score <- getFinalScore(result$response[, -"score"], result$finalResult$response)

        # UNDONE: what about the "..." params in IPO?

        # Sometimes the response from the parameter set returned by the model
        # doesn't actually give an optimum (see: https://github.com/rietho/IPO/issues/61).
        # In this case we simply take the best results from the experiments
        maxExpResp <- result$response[which.max(score)]
        if ((maxExpResp$score * (1 - maxModelDeviation)) > result$finalResult$score)
        {
            printf("Modelled parameter optimum yields significantly lower experimental score: %.1f/%.1f\n",
                   result$finalResult$score, result$max_settings[1])
            printf("Taking data from best experimental value (%.1f) instead...\n", maxExpResp$score)

            # runParams <- as.list(rsm::decode.data(result$design)[maxExpResp$experiment, names(params$to_optimize)])

            if (length(params$to_optimize) == 1)
                vals <- result$design[maxExpResp$experiment, 2, drop = FALSE] # design is just a df, see performIteration()
            else
            {
                # need to use as.data.frame here to avoid strange errors related to subsetting...
                vals <- as.data.frame(result$design)[maxExpResp$experiment,
                                                     paste0("x", seq_along(params$to_optimize)),
                                                     drop = FALSE]
            }

            runParams <- runParamsNoCall <- utilsIPO$decodeAll(vals, params$to_optimize)
            runParams <- Map(.self$fixDesignParam, names(runParams), runParams)
            runParams <- combineParams(runParams, params$no_optimization)

            # re-run as object wasn't stored
            runParams <- convertOptToCallParams(runParams)
            result$finalResult <- calculateResponse(runParams, 1, TRUE)
            result$finalResult$parameters <- runParams
            result$finalResult$score <- getFinalScore(result$response[, -"score"], result$finalResult$response)

            # update max_settings from experimental results
            result$max_settings <- sapply(seq_along(params$to_optimize),
                                          function(i) utilsIPO$encode(runParamsNoCall[[i]], params$to_optimize[[i]]))
            result$max_settings <- c(result$finalResult$score, result$max_settings)
            names(result$max_settings) <- c("response", paste0("x", seq_len(length(result$max_settings) - 1)))
            result$max_settings <- t(result$max_settings) # convert to usual data format...
        }

        printf("\n---\nResponse:\n")
        print(result$response)
        printf("---\n")

        cat("Best params: "); printf("%s: %s; ", names(runParams), runParams); cat("\n")
        br <- result$finalResult$response
        cat("Best results: "); printf("%s: %s; ", names(br), br); cat("\n")
        printf("---\n")

        return(result)
    },

    # heavily based on optimizeXcmsSet() from IPO
    optimize = function(allParams, templateParams, paramRanges, maxIterations, maxModelDeviation)
    {
        pSets <- lapply(seq_along(allParams), function(pi)
        {
            params <- modifyList(templateParams, allParams[[pi]])

            params <- startParams <- checkInitialParams(params)
            paramRanges <- modifyList(defaultParamRanges(params), paramRanges)

            history <- list()
            bestRange <- 0.25

            for (iter in seq_len(maxIterations))
            {
                if (length(allParams) > 1)
                    printf("Starting new DoE (parameter set %d/%d, iteration %d):\n", pi, length(allParams), iter)
                else
                    printf("Starting new DoE (iteration %d):\n", iter)
                cat(paste0(rbind(paste0(names(params), ": "),
                                 paste0(params, "\n"))), sep = "")

                result <- performIteration(params)
                result <- performIterationStat(result, maxModelDeviation)

                history[[iter]] <- result
                params <- result$params

                lastRun <- iter == maxIterations
                if (!resultIncreased(history) || lastRun)
                {
                    if (lastRun)
                        warning(sprintf("Reached maximum number of iterations (maxIterations=%d). Returning possibly suboptimal result.",
                                        maxIterations))

                    break
                }

                for (i in seq_len(length(params$to_optimize)))
                {
                    setting <- result$max_settings[i+1]
                    bounds <- params$to_optimize[[i]]
                    settingName <- names(params$to_optimize)[i]

                    settingRange <- getOptSettingRange(settingName, params, paramRanges)

                    # - if the parameter is NA, we increase the range by 20%,
                    # - if it was within the inner 25% of the previous range or
                    #   at the minimum value we decrease the range by 20%
                    if (is.na(setting))
                        stepFactor <- 1.2
                    else if (abs(setting) < bestRange ||
                             (setting == -1 && utilsIPO$decode(-1, params$to_optimize[[i]]) == settingRange[1]))
                        stepFactor <- 0.8
                    else
                        stepFactor <- 1

                    step <- (diff(bounds) / 2) * stepFactor

                    # CHANGED: check min-max for step
                    if (all(is.finite(settingRange)))
                    {
                        if ((diff(settingRange) / 2) < step)
                            printf("Changed step size to make sure its within range: %f --> %f\n", step, diff(settingRange) / 2)
                        step <- min(step, diff(settingRange) / 2)
                    }

                    if (is.na(setting))
                        setting <- 0

                    newCenter <- utilsIPO$decode(setting, bounds)

                    # CHANGED: also take max into account
                    # if ((newCenter-settingRange[1]) > step)
                    #     newBounds <- c(newCenter - step, newCenter + step)
                    # else
                    #     newBounds <- c(settingRange[1], 2*step+settingRange[1])

                    newBounds <- c(newCenter - step, newCenter + step)
                    if (is.finite(settingRange[1]) && newBounds[1] < settingRange[1])
                        newBounds <- c(settingRange[1], 2 * step + settingRange[1])
                    if (is.finite(settingRange[2]) && newBounds[2] > settingRange[2])
                    {
                        oldnb <- newBounds
                        newBounds <- c(settingRange[2] - 2 * step, settingRange[2])
                        printf("changed max range from (%f, %f) to (%f, %f)\n", oldnb[1], oldnb[2], newBounds[1], newBounds[2])
                    }

                    names(newBounds) <- NULL

                    params$to_optimize[[i]] <- fixOptParamBounds(names(params$to_optimize)[i], newBounds)
                }

                params <- fixOptParams(params)
                params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
            }

            #params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)

            maxima <- 0
            maxIndex <- 1
            for (i in seq_len(length(history)))
            {
                if (history[[i]]$max_settings[1] > maxima)
                {
                    maxima <- history[[i]]$max_settings[1]
                    maxIndex <- i
                }
            }

            finalParams <- as.list(utilsIPO$decodeAll(history[[maxIndex]]$max_settings[-1],
                                                      history[[maxIndex]]$params$to_optimize))
            finalParams <- combineParams(finalParams, params$no_optimization)

            if (!is.list(finalParams))
                finalParams <- as.list(finalParams)

            bestResults <- list()
            bestResults$parameters <- convertOptToCallParams(finalParams)

            bestResults$object <- history[[maxIndex]]$finalResult$object
            bestResults$DoEIteration <- maxIndex
            bestResults$score <- history[[maxIndex]]$finalResult$score

            return(list(startParams = startParams, bestResults = bestResults, iterations = history))
        })

        bestPS <- which.max(sapply(pSets, function(ps) ps$bestResults$score))

        printf("===\nDONE!\nBest parameter set: %d\nBest DoE iteration: %d\n",
               bestPS, pSets[[bestPS]]$bestResults$DoEIteration)

        return(list(paramSets = pSets, bestParamSet = bestPS))
    }
)

#' Class containing optimization results.
#'
#' Objects from this class contain optimization results resulting from design of
#' experiment (DoE).
#'
#' Objects from this class are returned by \code{\link{optimizeFeatureFinding}} and
#' \code{\link{optimizeFeatureGrouping}}.
#'
#' @slot algorithm A character specifying the algorithm that was optimized.
#' @slot paramSets A \code{list} with detailed results from each parameter set
#'   that was tested.
#' @slot bestParamSet Numeric index of the parameter set yielding the best
#'   response.
#'
#' @param obj,x,object An \code{optimizationResult} object.
#' @param paramSet Numeric index of the parameter set (\emph{i.e.} the first
#'   parameter set gets index \samp{1}). For some methods optional: if
#'   \code{NULL} the best will be selected.
#' @param DoEIteration Numeric index specifying the DoE iteration within the
#'   specified \code{paramSet}. For some methods optional: if \code{NULL} the
#'   best will be selected.
#'
#' @export
optimizationResult <- setClass("optimizationResult",
                               slots = c(algorithm = "character",
                                         paramSets = "list", bestParamSet = "numeric"))

#' @describeIn optimizationResult Returns the algorithm that was used for finding features.
setMethod("algorithm", "optimizationResult", function(obj) obj@algorithm)

#' @describeIn optimizationResult Obtain total number of experimental design iterations performed.
#' @export
setMethod("length", "optimizationResult", function(x) sum(lengths(x)))

#' @describeIn optimizationResult Obtain number of experimental design iterations performed for each parameter set.
#' @param use.names Ignored.
#' @export
setMethod("lengths", "optimizationResult", function(x, use.names = FALSE) sapply(x@paramSets, function(ps) length(ps$iterations),
                                                                                 USE.NAMES = use.names))
#' @describeIn optimizationResult Shows summary information for this object.
#' @export
setMethod("show", "optimizationResult", function(object)
{
    printf("An optimization result object ('%s')\n", class(object))
    printf("Algorithm: %s\n", algorithm(object))

    for (pi in seq_along(object@paramSets))
    {
        ps <- object@paramSets[[pi]]
        printf("Parameter set %d/%d%s:\n", pi, length(object@paramSets),
               if (pi == object@bestParamSet) " (BEST)" else "")

        printf("    Experimental designs tested: %d\n", length(ps$iterations))
        printf("    Starting params:\n"); printf("    - %s: %s\n", names(ps$startParams), ps$startParams)
        printf("    Optimized params:\n"); printf("    - %s: %s\n", names(ps$bestResults$parameters), ps$bestResults$parameters)

        bexp <- ps$iterations[[ps$bestResults$DoEIteration]]
        br <- bexp$finalResult$response
        printf("    Best results: "); cat(paste(names(br), br, sep = ": ", collapse = "; ")); cat("\n")

    }

    printf("Best parameter set: %d\n", object@bestParamSet)

    # printf("\nOptimized object:\n---\n"); show(object@bestResults$object); cat("---\n")

    showObjectSize(object)
})

#' @describeIn optimizationResult Generates response plots for all or a selected
#'   set of parameters.
#'
#' @param paramsToPlot Which parameters relations should be plot. If \code{NULL}
#'   all will be plot. Alternatively, a \code{list} containing one or more
#'   \code{character} vectors specifying each two parameters that should be
#'   plotted. Finally, if only one pair should be plotted, can be a
#'   \code{character} vector specifying both parameters.
#' @param maxCols Multiple parameter pairs are plotted in a grid. The maximum
#'   number of columns can be set with this argument. Set to \code{NULL} for no
#'   limit.
#' @param type The type of plots to be generated: \code{"contour"},
#'   \code{"image"} or \code{"persp"}. The equally named functions will be
#'   called for plotting.
#' @param image Passed to \code{\link{contour}} (if \code{type="contour"}).
#' @param contours Passed to \code{\link{persp}} (if \code{type="persp"}).
#' @param \dots Further arguments passed to \code{\link{contour}},
#'   \code{\link{image}} or \code{\link{persp}} (depending on \code{type}).
#'
#' @examples \dontrun{
#' # ftOpt is an optimization object.
#'
#' # plot contour of all parameter pairs from the first parameter set/iteration.
#' plot(ftOpt, paramSet = 1, DoEIteration = 1)
#'
#' # as above, but only plot two parameter pairs
#' plot(ftOpt, paramSet = 1, DoEIteration = 1,
#'      paramsToPlot = list(c("mzPPM", "chromFWHM"), c("chromFWHM", "chromSNR")))
#'
#' # plot 3d perspective plots
#' plot(ftOpt, paramSet = 1, DoEIteration = 1, type = "persp")
#' }
#'
#' @export
setMethod("plot", "optimizationResult", function(x, paramSet, DoEIteration, paramsToPlot = NULL, maxCols = NULL, type = "contour",
                                                 image = TRUE, contours = "colors", ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertInt(paramSet, lower = 1, upper = length(x@paramSets)) # don't add, this should fail before the next line
    checkmate::assertInt(DoEIteration, lower = 1, upper = length(x@paramSets[[paramSet]]$iterations), add = ac)
    checkmate::assert(checkmate::checkList(paramsToPlot, types = "character", any.missing = FALSE, min.len = 1, null.ok = TRUE),
                      checkmate::checkCharacter(paramsToPlot, min.chars = 1, len = 2),
                      checkmate::checkNull(paramsToPlot),
                      .var.name = "paramsToPlot")
    checkmate::assertCount(maxCols, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::assertChoice(type, c("contour", "image", "persp"), add = ac)
    checkmate::assertFlag(image, add = ac)
    checkmate::assert(checkmate::checkFlag(contours),
                      checkmate::checkCharacter(contours),
                      checkmate::checkList(contours),
                      .var.name = "contours")
    checkmate::reportAssertions(ac)

    ex <- x@paramSets[[paramSet]]$iterations[[DoEIteration]]

    if (length(ex$params$to_optimize) < 2)
        stop("Need at least two optimized parameters for plotting.")

    if (is.null(paramsToPlot))
    {
        paramsToPlot <- list()
        optNames <- names(ex$params$to_optimize)
        for (i in seq(1, length(optNames)-1))
        {
            for (j in seq(i+1, length(optNames)))
                paramsToPlot <- c(paramsToPlot, list(c(optNames[i], optNames[j])))
        }
    }
    else if (is.character(paramsToPlot))
        paramsToPlot <- list(paramsToPlot)

    codedNames <- names(ex$design)
    decodedNames <- rsm::truenames(ex$design)

    forms <- lapply(paramsToPlot, function(pn)
    {
        # change to coded names
        pn <- sapply(pn, function(n) codedNames[decodedNames == n])
        return(as.formula(paste(pn[2], "~", pn[1])))
    })

    maxSlice <- ex$max_settings[1, -1]
    maxSlice[is.na(maxSlice)] <- 1

    formsLen <- length(forms)
    if (formsLen > 1) # multiple plots?
    {
        if (is.null(maxCols))
            maxCols <- ceiling(sqrt(formsLen))

        if (formsLen <= maxCols)
        {
            cols <- formsLen
            rows <- 1
        }
        else
        {
            cols <- maxCols
            rows <- ceiling(formsLen / cols)
        }

        withr::local_par(list(mfrow = c(rows, cols)))
    }

    switch(type,
           contour = contour(ex$model, forms, image = image, at = maxSlice, ...),
           image = image(ex$model, forms, at = maxSlice, ...),
           persp = persp(ex$model, forms, contours = contours, at = maxSlice, ...))
})

#' @describeIn optimizationResult Returns parameter set yielding optimal
#'   results. The \code{paramSet} and \code{DoEIteration} arguments can be
#'   \code{NULL}.
#' @aliases optimizedParameters
#' @export
setMethod("optimizedParameters", "optimizationResult", function(object, paramSet, DoEIteration)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertInt(paramSet, lower = 1, upper = length(object@paramSets),
                         null.ok = is.null(DoEIteration), add = ac)
    checkmate::assertCount(DoEIteration, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (is.null(paramSet))
        paramSet <- object@bestParamSet

    if (!is.null(DoEIteration))
        return(object@paramSets[[paramSet]]$iterations[[DoEIteration]]$finalResult$parameters)

    return(object@paramSets[[paramSet]]$bestResults$parameters)
})

#' @describeIn optimizationResult Returns the object (\emph{i.e.} a
#'   \code{\link{features}} or \code{\link{featureGroups}} object) that was
#'   generated with optimized parameters. The \code{paramSet} argument can be
#'   \code{NULL}.
#' @aliases optimizedObject
#' @export
setMethod("optimizedObject", "optimizationResult", function(object, paramSet)
{
    checkmate::assertInt(paramSet, lower = 1, upper = length(object@paramSets), null.ok = TRUE)

    if (is.null(paramSet))
        paramSet <- object@bestParamSet

    return(object@paramSets[[paramSet]]$bestResults$object)
})


#' @describeIn optimizationResult Returns optimization scores. The
#'   \code{paramSet} and \code{DoEIteration} arguments can be \code{NULL}.
#' @aliases scores
#' @export
setMethod("scores", "optimizationResult", function(object, paramSet, DoEIteration)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertInt(paramSet, lower = 1, upper = length(object@paramSets),
                         null.ok = is.null(DoEIteration), add = ac)
    checkmate::assertCount(DoEIteration, positive = TRUE, null.ok = TRUE, add = ac)
    checkmate::reportAssertions(ac)

    if (is.null(paramSet))
        paramSet <- object@bestParamSet

    if (is.null(DoEIteration))
        DoEIteration <- object@paramSets[[paramSet]]$bestResults$DoEIteration

    return(object@paramSets[[paramSet]]$iterations[[DoEIteration]]$finalResult$response)
})

#' @describeIn optimizationResult Returns a \code{list} with optimization
#'   information from an DoE iteration.
#' @aliases experimentInfo
#' @export
setMethod("experimentInfo", "optimizationResult", function(object, paramSet, DoEIteration)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertInt(paramSet, lower = 1, upper = length(object@paramSets)) # don't add, this should fail before the next line
    checkmate::assertInt(DoEIteration, lower = 1, upper = length(object@paramSets[[paramSet]]$iterations), add = ac)
    checkmate::reportAssertions(ac)

    return(object@paramSets[[paramSet]]$iterations[[DoEIteration]])
})
