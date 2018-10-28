#' @include utils-IPO.R
#' @include main.R
NULL

DoEOptimizer <- setRefClass("DoEOptimizer", contains = "VIRTUAL",
                            fields = list(algorithm = "character",
                                          maxModelDeviation = "numeric"))

DoEOptimizer$methods(

    # dummy methods that may need to be overloaded
    checkInitialParams = function(params) params,
    getMinOptSetting = function(settingName, params) 0,
    fixParamBounds = function(params, bounds) bounds,
    fixOptParams = function(params) params,

    # "virtual" methods
    resultIncreased = function(history) stop("VIRTUAL"),
    calculateResponse = function(params, task, final = FALSE) stop("VIRTUAL"),
    getResponseScores = function(response) stop("VIRTUAL"),
    getFinalScore = function(response) stop("VIRTUAL"),


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
        else
        {
            design <- data.frame(run.order = 1:9, a = seq(-1,1,0.25))
            colnames(design)[2] <- names(typParams$to_optimize)
            designParams <- design
            designParams[,2] <- seq(min(typParams$to_optimize[[1]]),
                                    max(typParams$to_optimize[[1]]),
                                    diff(typParams$to_optimize[[1]]) / 8)
        }

        designParams <- combineParams(designParams, typParams$no_optimization)
        tasks <- seq_len(nrow(design))

        printf("---\nDesign:\n")
        print(rsm::decode.data(design))
        printf("---\n")

        prog <- txtProgressBar(0, length(tasks), style = 3)

        response <- rbindlist(lapply(tasks, function(task)
        {
            # simplified optimizeSlaveCluster() from IPO

            runParams <- as.list(designParams[task, ])
            runParams <- runParams[!names(runParams) %in% c("run.order", "std.order", "Block")]
            result <- calculateResponse(runParams, task)
            result$experiment <- task

            setTxtProgressBar(prog, task)

            return(result)
        }))

        setTxtProgressBar(prog, length(tasks))
        close(prog)

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

        result$finalResponse <- calculateResponse(runParams, 1, final = TRUE)
        result$finalResponse$score <- getFinalScore(result$response[, -"score"], result$finalResponse)

        # UNDONE: what about the "..." params in IPO?

        # Sometimes the response from the parameter set returned by the model
        # doesn't actually give an optimum (see: https://github.com/rietho/IPO/issues/61).
        # In this case we simply take the best results from the experiments
        maxExpResp <- result$response[which.max(score)]
        if ((maxExpResp$score * (1 - maxModelDeviation)) > result$finalResponse$score)
        {
            printf("Modelled parameter optimum yields significantly lower experimental score: %.1f/%.1f\n",
                   result$finalResponse$score, result$max_settings[1])
            printf("Taking data from best experimental value (%.1f) instead...\n", maxExpResp$score)

            runParams <- as.list(rsm::decode.data(result$design)[maxExpResp$experiment, names(params$to_optimize)])

            # re-run as object wasn't stored
            result$finalResponse <- calculateResponse(combineParams(runParams, params$no_optimization), 1, final = TRUE)
            result$finalResponse$score <- getFinalScore(result$response[, -"score"], result$finalResponse)

            # update max_settings from experimental results
            result$max_settings <- sapply(seq_along(params$to_optimize),
                                          function(i) utilsIPO$encode(runParams[[i]], params$to_optimize[[i]]))
            result$max_settings <- c(result$finalResponse$score, result$max_settings)
            names(result$max_settings) <- c("response", names(rsm::codings(result$design)))
            result$max_settings <- t(result$max_settings) # convert to usual data format...
        }

        printf("\n---\nResponse:\n")
        print(result$response)
        printf("---\n")

        cat("Best params: "); printf("%s: %s; ", names(runParams), runParams); cat("\n")
        br <- result$finalResponse[names(result$finalResponse) != "object"]
        cat("Best results: "); printf("%s: %s; ", names(br), br); cat("\n")
        printf("---\n")

        return(result)
    },

    # heavily based on optimizeXcmsSet() from IPO
    optimize = function(params, maxIterations, maxModelDeviation)
    {
        params <- startParams <- checkInitialParams(params)

        history <- list()
        bestRange <- 0.25

        for (iter in seq_len(maxIterations))
        {
            printf("Starting new DoE with (#%d):\n", iter)
            printf(paste0(rbind(paste0(names(params), ": "),
                                paste0(params, "\n"))))

            result <- performIteration(params)
            result <- performIterationStat(result, maxModelDeviation)

            history[[iter]] <- result
            params <- result$params

            if (!resultIncreased(history))
            {
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

                bestSettings <- list()
                bestSettings$parameters <- finalParams

                bestSettings$object <- history[[maxIndex]]$finalResponse$object
                bestSettings$result <- history[[maxIndex]]$finalResponse[names(history[[maxIndex]]$finalResponse) != "object"]
                bestSettings$score <- history[[maxIndex]]$finalResponse$score
                history$bestSettings <- bestSettings

                break
            }

            for (i in seq_len(length(params$to_optimize)))
            {
                setting <- result$max_settings[i+1]
                bounds <- params$to_optimize[[i]]
                settingName <- names(params$to_optimize)[i]

                minSetting <- getMinOptSetting(settingName, params)

                # - if the parameter is NA, we increase the range by 20%,
                # - if it was within the inner 25% of the previous range or
                #   at the minimum value we decrease the range by 20%
                if (is.na(setting))
                    stepFactor <- 1.2
                else if (abs(setting) < bestRange ||
                         (setting == -1 && utilsIPO$decode(-1, params$to_optimize[[i]]) == minSetting))
                    stepFactor <- 0.8
                else
                    stepFactor <- 1

                step <- (diff(bounds) / 2) * stepFactor

                if (is.na(setting))
                    setting <- 0

                newCenter <- utilsIPO$decode(setting, bounds)

                if ((newCenter-minSetting) > step)
                    newBounds <- c(newCenter - step, newCenter + step)
                else
                    newBounds <- c(minSetting, 2*step+minSetting)

                names(newBounds) <- NULL

                params$to_optimize[[i]] <- fixParamBounds(names(params$to_optimize)[i], newBounds)
            }

            params <- fixOptParams(params)
            params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)
        }

        #params <- utilsIPO$attachList(params$to_optimize, params$no_optimization)

        return(list(startParams = startParams, finalResults = history$bestSettings, experiments = history[seq_len(iter)]))
    }
)

# based on part of optimizeXcmsSet() function from IPO
fixOptParamRange <- function(params, paramPairs)
{
    for (pp in paramPairs)
    {
        if (!is.null(params$to_optimize[[pp[1]]]) || !is.null(params$to_optimize[[pp[2]]]))
        {
            if (is.null(params$to_optimize[[pp[1]]]))
                pmin <- params$no_optimization[[pp[1]]]
            else
                pmin <- max(params$to_optimize[[pp[1]]])

            if (is.null(params$to_optimize[[pp[2]]]))
                pmax <- params$no_optimization[[pp[2]]]
            else
                pmax <- min(params$to_optimize[[pp[2]]])

            if (pmin >= pmax)
            {
                additional <- abs(pmin-pmax) + 1
                if (!is.null(params$to_optimize[[pp[2]]]))
                    params$to_optimize[[pp[2]]] <- params$to_optimize[[pp[2]]] + additional
                else
                    params$no_optimization[[pp[2]]] <- params$no_optimization[[pp[2]]] + additional
            }
        }
    }

    return(params)
}


#' @export
optimizationResult <- setClass("optimizationResult",
                               slots = c(algorithm = "character",
                                         startParams = "list", finalResults = "list",
                                         experiments = "list"))

#' @describeIn optimizationResult Returns the algorithm that was used for finding features.
setMethod("algorithm", "optimizationResult", function(obj) obj@algorithm)

#' @describeIn optimizationResult Obtain total number of experimental design iteratations performed.
#' @export
setMethod("length", "optimizationResult", function(x) length(x@experiments))

#' @describeIn optimizationResult Shows summary information for this object.
#' @export
setMethod("show", "optimizationResult", function(object)
{
    printf("An optimization result object ('%s')\n", class(object))
    printf("Algorithm: %s\n", algorithm(object))
    printf("Experimental designs performed: %d\n", length(object))
    printf("Starting params:\n"); printf("- %s: %s\n", names(object@startParams), object@startParams)
    printf("Optimized params:\n"); printf("- %s: %s\n", names(object@finalResults$parameters), object@finalResults$parameters)

    br <- object@finalResults$result
    br <- br[!names(br) %in% "ExpId"]
    printf("Best results: "); cat(paste(names(br), br, sep = ": ", collapse = "; ")); cat("\n")

    printf("\nOptimized object:\n---\n"); show(object@finalResults$object); cat("---\n")

    showObjectSize(object)
})

#' @export
setMethod("plot", "optimizationResult", function(x, index, paramsToPlot = NULL, maxCols = NULL, type = "contour",
                                                 image = TRUE, contours = "colors", ...)
{
    ac <- checkmate::makeAssertCollection()
    checkmate::assertInt(index, lower = 1, upper = length(x))
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

    ex <- x@experiments[[index]]

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

