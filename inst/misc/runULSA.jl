using ULSA
using JLD
using CSV
using DecisionTree

# based on example: https://bitbucket.org/SSamanipour/ulsa.jl/src/master/examples/USLA_run.jl

lastArg = 1
function popArg()
    global lastArg
    ret = ARGS[lastArg]
    lastArg += 1
    return ret
end

mode = popArg()
source = popArg()
weightF = [parse(Float64, popArg()), parse(Float64, popArg()), parse(Float64, popArg()), parse(Float64, popArg()), parse(Float64, popArg()),
    parse(Float64, popArg()), parse(Float64, popArg())]
inputFile = popArg()

featureID_comp(mode, source, inputFile, weightF)
