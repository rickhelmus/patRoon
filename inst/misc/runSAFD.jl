# SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
#
# SPDX-License-Identifier: GPL-3.0-only

using SAFD
using MS_Import
using Printf
using DelimitedFiles

# based on example script: https://bitbucket.org/SSamanipour/safd.jl/src/master/examples/SAFD_Single_run.jl

lastArg = 1
function popArg()
    global lastArg
    ret = ARGS[lastArg]
    lastArg += 1
    return ret
end

inPath = popArg()
filename = popArg()

doCent = popArg() == "TRUE"

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

centMethod = popArg()
centDM = parse(Float64, popArg())

mzCSV = popArg()
intCSV = popArg()
rtCSV = popArg()

@printf("CSV Files: %s / %s / %s\n", mzCSV, intCSV, rtCSV)

outPath = popArg()

GC.gc()

#mz_vals, mz_int, t0, t_end, m, pathin, msModel, msIonisation, msManufacturer, polarity, Rt = import_files_MS1(inPath, filename, mz_thresh)
mz_vals = readdlm(mzCSV, ',', Float32)
mz_int = readdlm(intCSV, ',', Float32)
Rt = vec(readdlm(rtCSV, ',', Float32))

GC.gc()

if (doCent)
    safd_s3d_cent(mz_vals, mz_int, Rt, filename, outPath, max_numb_iter, max_t_peak_w, res, min_ms_w, r_thresh, min_int,
                  sig_inc_thresh, S2N, min_peak_w_s, centMethod, (centMethod == "CT") ? centDM : [])
else
    safd_s3D(mz_vals, mz_int, Rt, filename, outPath, max_numb_iter, max_t_peak_w, res, min_ms_w, r_thresh, min_int,
             sig_inc_thresh, S2N, min_peak_w_s)
end
