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

            # CHANGED: for now ignore fixing the range if the min/max is not known (unlike IPO default params can be omitted)
            if (is.null(pmin) || is.null(pmax))
                next
            
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

assertOptimArgs <- function(params, templateParams, paramRanges, maxIterations, maxModelDeviation, parallel, ac)
{
    aapply(checkmate::assertList, ~ templateParams + paramRanges, names = "unique",
           fixed = list(add = ac))
    checkmate::assertList(params, min.len = 1, add = ac)
    checkmate::assertCount(maxIterations, positive = TRUE, add = ac)
    checkmate::assertNumber(maxModelDeviation, finite = TRUE, add = ac)
    checkmate::assertFlag(parallel, add = ac)

    checkmate::qassertr(params, "l+", "...")
}
