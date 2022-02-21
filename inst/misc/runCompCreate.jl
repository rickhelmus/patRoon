using CompCreate
using MS_Import
using JLD

# based on example: https://bitbucket.org/SSamanipour/compcreate.jl/src/master/examples/CompCreateSingle.jl

lastArg = 1
function popArg()
    global lastArg
    ret = ARGS[lastArg]
    lastArg += 1
    return ret
end

inPath = popArg()
filename = popArg()
mz_thresh = [parse(Float64, popArg()), parse(Float64, popArg())]

onlyMS = popArg() == "TRUE"
mass_win_per = parse(Float64, popArg())
ret_win_per = parse(Float64, popArg())
r_thresh = parse(Float64, popArg())
delta_mass = parse(Float64, popArg())
min_int = parse(Int64, popArg()) # Not needed for comp_ms1()

featureInput = popArg()

GC.gc()

chrom = import_files(inPath, [ filename ], mz_thresh)

if onlyMS
    comp_ms1(chrom, featureInput, mass_win_per, ret_win_per, r_thresh, delta_mass)
else
    compcreate(chrom, featureInput, mass_win_per, ret_win_per, r_thresh, delta_mass, min_int)
end
