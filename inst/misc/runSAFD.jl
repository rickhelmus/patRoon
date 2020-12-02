using SAFD
using MS_Import
using Printf

# based on example script: https://bitbucket.org/SSamanipour/safd.jl/src/master/examples/SAFD_Single_run.jl

lastArg = 1
function popArg()
    global lastArg
    ret = ARGS[lastArg]
    lastArg += 1
    return ret
end

inPath = popArg()
outPath = popArg()
filename = popArg()

mz_thresh = [parse(Float64, popArg()), parse(Float64, popArg())]

max_numb_iter = parse(Int64, popArg())
max_t_peak_w = parse(Int64, popArg())
res = parse(Int64, popArg())
min_ms_w = parse(Float64, popArg())
r_thresh = parse(Float64, popArg())
min_int = parse(Int64, popArg())
sig_inc_thresh = parse(Int64, popArg())
S2N = parse(Int64, popArg())
min_peak_w_s = parse(Int64, popArg())

GC.gc()

mz_vals, mz_int, t0, t_end, m, pathin = import_files_MS1(inPath, filename, mz_thresh)

GC.gc()

rep_table, final_table = safd_s3D(mz_vals, mz_int, t0, t_end, filename, outPath, max_numb_iter,
    max_t_peak_w, res, min_ms_w, r_thresh, min_int, sig_inc_thresh, S2N, min_peak_w_s)
