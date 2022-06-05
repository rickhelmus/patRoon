/*
 *   OpenTIMS: a fully open-source library for opening Bruker's TimsTOF data files.
 *   Copyright (C) 2020-2021 Michał Startek and Mateusz Łącki
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License, version 3 only,
 *   as published by the Free Software Foundation.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <cassert>
#include <cstdint>
#include <string>
#include <cstring>
#include <atomic>
#include <vector>
#include <iostream>
#include <locale>
#include <memory>
#include <limits>
#include <thread>
#include <unordered_map>


#include "platform.h"

#ifdef OPENTIMS_BUILDING_R
#include "mio.h"
#include "opentimsr_types.h"
#else
#include "mio.hpp"
#include "opentims.h"
#endif

#include "tof2mz_converter.h"
#include "scan2inv_ion_mobility_converter.h"
#include "thread_mgr.h"

TimsFrame::TimsFrame(uint32_t _id,
                     uint32_t _num_scans,
                     uint32_t _num_peaks,
                     uint32_t _msms_type,
                     double _intensity_correction,
                     double _time,
                     const char* frame_ptr,
                     TimsDataHandle& parent_hndl
                     )
:
    bytes0(nullptr),
    tims_bin_frame(frame_ptr),
    parent_tdh(parent_hndl),
    id(_id),
    num_scans(_num_scans),
    num_peaks(_num_peaks),
    msms_type(_msms_type),
    intensity_correction(_intensity_correction),
    time(_time)
{}

TimsFrame TimsFrame::TimsFrameFromSql(char** sql_row, TimsDataHandle& parent_handle)
{
    assert(sql_row != nullptr);
    assert(sql_row[0] != nullptr);
    assert(sql_row[1] != nullptr);
    assert(sql_row[2] != nullptr);
    assert(sql_row[3] != nullptr);
    assert(sql_row[4] != nullptr);
    assert(sql_row[5] != nullptr);
    assert(sql_row[6] != nullptr);

    return TimsFrame(
            atol(sql_row[0]),
            atol(sql_row[1]),
            atol(sql_row[2]),
            atol(sql_row[3]),
            100.0 / atof(sql_row[4]),
            atof(sql_row[5]),
            std::strtoul(sql_row[6], nullptr, 10) + parent_handle.tims_data_bin.data(),
            parent_handle
    );
}

void TimsFrame::print() const
{
#ifndef OPENTIMS_BUILDING_R
    std::cout << "Frame description: id: " << id << ", num_scans: " << num_scans << ", num_peaks: " << num_peaks << std::endl;
#endif
}


void TimsFrame::decompress(char* decompression_buffer, ZSTD_DCtx* decomp_ctx)
{
    uint32_t tims_packet_size = *reinterpret_cast<const uint32_t*>(tims_bin_frame);
    assert(num_scans == *(reinterpret_cast<const uint32_t*>(tims_bin_frame)+1));

    size_t dsbytes = data_size_bytes();

    if(decompression_buffer == nullptr)
    {
        decompression_buffer = parent_tdh.decompression_buffer.get();
//        back_buffer = std::make_unique<char[]>(dsbytes);
//        decompression_buffer = reinterpret_cast<char*>(back_buffer.get());
    }

    if(decomp_ctx == nullptr)
        decomp_ctx = parent_tdh.zstd_dctx;

    size_t dec_result = ZSTD_decompressDCtx(decomp_ctx, decompression_buffer, dsbytes, tims_bin_frame + 8, tims_packet_size - 8);
    if(ZSTD_isError(dec_result))
    {
        std::string err = "Error uncompressing frame, error code: ";
        err += std::to_string(dec_result);
        err += ". File is either corrupted, or in a (yet) unsupported variant of the format.";
        throw std::runtime_error(err);
    }

    size_t dsints = data_size_ints();
    bytes0 = decompression_buffer;
    bytes1 = bytes0 + dsints;
    bytes2 = bytes1 + dsints;
    bytes3 = bytes2 + dsints;
}

void TimsFrame::close()
{
    bytes0 = nullptr;
    back_buffer.reset(nullptr);
}

void TimsFrame::save_to_buffs(uint32_t* frame_ids,
                              uint32_t* scan_ids,
                              uint32_t* tofs,
                              uint32_t* intensities,
                              double* mzs,
                              double* inv_ion_mobilities,
                              double* retention_times,
                              ZSTD_DCtx* decomp_ctx)
{
    if(num_peaks == 0)
        return;

    std::unique_ptr<uint32_t[]> scan_ids_hndl;
    std::unique_ptr<uint32_t[]> tofs_hndl;
    std::unique_ptr<uint32_t[]> intensities_hndl;

    if(scan_ids == nullptr && inv_ion_mobilities != nullptr)
    {
        scan_ids_hndl = std::make_unique<uint32_t[]>(num_peaks);
        scan_ids = scan_ids_hndl.get();
    }
    if(tofs == nullptr)
    {
        tofs_hndl = std::make_unique<uint32_t[]>(num_peaks);
        tofs = tofs_hndl.get();
    }
    if(intensities == nullptr)
    {
        intensities_hndl = std::make_unique<uint32_t[]>(num_peaks);
        intensities = intensities_hndl.get();
    }

    bool needs_closure = false;
    if(bytes0 == nullptr)
    {
        decompress(nullptr, decomp_ctx);
        needs_closure = true;
    }

    uint32_t peaks_processed = 0;
    size_t read_offset = num_scans;
    uint32_t accum_tofs;

    uint32_t num_scans_m1 = num_scans - 1;

    for(uint32_t scan_idx = 0; scan_idx < num_scans_m1; scan_idx++)
    {
        accum_tofs = -1;

        const uint32_t no_peaks = back_data(scan_idx+1) / 2;

        const uint32_t for_loop_end = no_peaks + peaks_processed;

        if(scan_ids != nullptr)
            for(uint32_t ii = peaks_processed; ii < for_loop_end; ii++)
                scan_ids[ii] = scan_idx;

        for(uint32_t ii = 0; ii < no_peaks; ii++)
        {
            accum_tofs += back_data(read_offset);
            tofs[peaks_processed] = accum_tofs;
            read_offset++;
            intensities[peaks_processed] = back_data(read_offset);
            read_offset++;
            peaks_processed++;
        }
    }

    accum_tofs = -1;

    const uint32_t nnum_peaks = num_peaks;

    if(scan_ids != nullptr)
        for(uint32_t ii = peaks_processed; ii < nnum_peaks; ii++)
            scan_ids[ii] = num_scans_m1;

    while(peaks_processed < nnum_peaks)
    {
        accum_tofs += back_data(read_offset);
        tofs[peaks_processed] = accum_tofs;
        read_offset++;
        intensities[peaks_processed] = back_data(read_offset);
        read_offset++;
        peaks_processed++;
    }


    for(size_t idx = 0; idx < nnum_peaks; idx++)
        intensities[idx] = static_cast<double>(intensities[idx]) * intensity_correction + 0.5;

    if(mzs != nullptr)
        parent_tdh.tof2mz_converter->convert(id, mzs, tofs, nnum_peaks);

    if(frame_ids != nullptr)
        for(size_t idx = 0; idx < nnum_peaks; idx++)
            frame_ids[idx] = id;

    if(retention_times != nullptr)
        for(size_t idx = 0; idx < nnum_peaks; idx++)
            retention_times[idx] = time;

    if(inv_ion_mobilities != nullptr)
        parent_tdh.scan2inv_ion_mobility_converter->convert(id, inv_ion_mobilities, scan_ids, nnum_peaks);

    if(needs_closure)
        close();
}

int tims_sql_callback(void* out, [[maybe_unused]] int cols, char** row, char**)
{
    assert(cols == 7);
    assert(row != NULL);
    assert(row[0] != NULL);
    uint32_t frame_id = atol(row[0]);
    TimsDataHandle* hndl = reinterpret_cast<TimsDataHandle*>(out);
    hndl->frame_descs.emplace(frame_id, TimsFrame::TimsFrameFromSql(row, *hndl));
    return 0;
}

int check_compression(void*, [[maybe_unused]] int cols, char** row, char**)
{
    assert(cols == 1);
    assert(row != NULL);
    assert(row[0] != NULL);
    if(atoi(row[0]) != 2)
    {
        std::string error_msg = "Compression algorithm used in your TDF dataset: ";
        error_msg += row[0];
        error_msg += " is not (yet) supported by OpenTIMS. Right now only algorithm 2 (zstd) is supported.";
        throw std::runtime_error(error_msg);
    }
    return 0;
}

#ifndef OPENTIMS_BUILDING_R
namespace{
class RAIILocaleHelper
{
    const std::locale previous_locale;
 public:
    RAIILocaleHelper() : previous_locale(std::locale::global(std::locale("C"))) {};
    ~RAIILocaleHelper() { std::locale::global(previous_locale); };
};

class RAIISqlite
{
    sqlite3* db_conn;
 public:
    RAIISqlite(const std::string& tims_tdf_path) : db_conn(nullptr)
    {
        if(sqlite3_open_v2(tims_tdf_path.c_str(), &db_conn, SQLITE_OPEN_READONLY, NULL))
            throw std::runtime_error(std::string("ERROR opening database: " + tims_tdf_path + " SQLite error msg: ") + sqlite3_errmsg(db_conn));
    }
    ~RAIISqlite()
    {
        if(db_conn != nullptr)
            sqlite3_close(db_conn);
    }
    sqlite3* release_connection() { sqlite3* ret = db_conn; db_conn = nullptr; return ret; }
    void query(const std::string& sql, int (*callback)(void*,int,char**,char**), void* arg)
    {
        char* error = NULL;

        if(sqlite3_exec(db_conn, sql.c_str(), callback, arg, &error) != SQLITE_OK)
        {
	    std::string err_msg(std::string("ERROR performing SQL query. SQLite error msg: ") + error);
	    sqlite3_free(error);
	    throw std::runtime_error(err_msg);
        }
    }

};
}
#endif

void TimsDataHandle::read_sql(const std::string& tims_tdf_path)
{
#ifndef OPENTIMS_BUILDING_R
    RAIILocaleHelper locale_guard;
    RAIISqlite DB(tims_tdf_path);

    const std::string sql = "SELECT Id, NumScans, NumPeaks, MsMsType, AccumulationTime, Time, TimsId from Frames;";

    DB.query(sql, tims_sql_callback, this);
    DB.query("SELECT Value FROM GlobalMetadata WHERE Key == \"TimsCompressionType\";", check_compression, nullptr);

    db_conn = DB.release_connection();

#endif
}


void TimsDataHandle::init()
{
    _min_frame_id = (std::numeric_limits<uint32_t>::max)();
    _max_frame_id = (std::numeric_limits<uint32_t>::min)();
    decomp_buffer_size = 0;
    for(auto it = frame_descs.begin(); it != frame_descs.end(); it++)
    {
        _min_frame_id = (std::min)(_min_frame_id, it->first);
        _max_frame_id = (std::max)(_max_frame_id, it->first);
        decomp_buffer_size = (std::max)(decomp_buffer_size, it->second.data_size_bytes());
    }
    decompression_buffer = std::make_unique<char[]>(decomp_buffer_size);

    zstd_dctx = ZSTD_createDCtx();

    tof2mz_converter = DefaultTof2MzConverterFactory::produceDefaultConverterInstance(*this);
    scan2inv_ion_mobility_converter = DefaultScan2InvIonMobilityConverterFactory::produceDefaultConverterInstance(*this);
}

TimsDataHandle::TimsDataHandle(const std::string& tims_tdf_bin_path, const std::string& tims_tdf_path, const std::string& tims_data_dir)
: tims_dir_path(tims_data_dir), tims_data_bin(tims_tdf_bin_path), zstd_dctx(nullptr)
#ifndef OPENTIMS_BUILDING_R
, db_conn(nullptr)
#endif
{
#ifndef OPENTIMS_BUILDING_R
    read_sql(tims_tdf_path);
#endif
    init();
}

TimsDataHandle::TimsDataHandle(const std::string& tims_data_dir)
: TimsDataHandle(tims_data_dir + "/analysis.tdf_bin", tims_data_dir + "/analysis.tdf", tims_data_dir)
{}

#ifdef OPENTIMS_BUILDING_R

union braindead_r
{
    double as_dbl;
    int64_t as_int;
};

template<typename T> std::vector<T> braindead_r_extract_as_int(const SEXP& vec)
{
    /* Quick and dirty workaround for R being braindead.
     * R holds integers as either IntegerVector or as NumericVector depending on size.
     * Converting a NumericVector to IntegerVector fails silently, giving zeros.
     * This function forces it to a sensible C++ int type.
     */

    std::vector<T> res;

    if(Rf_isInteger(vec))
    {
        Rcpp::IntegerVector vint(vec);
        res.reserve(vint.size());
        for(int ii=0; ii<vint.size(); ii++)
            res.push_back(vint[ii]);
        return res;
    }
    else
    {
        braindead_r converter;
        Rcpp::NumericVector vnum(vec);
        res.reserve(vnum.size());
        for(int ii=0; ii<vnum.size(); ii++)
        {
            converter.as_dbl = vnum[ii];
            res.push_back(converter.as_int);
        }
        return res;
    }
}

TimsDataHandle::TimsDataHandle(const std::string& tims_data_dir, const Rcpp::List& analysis_tdf) :
TimsDataHandle(tims_data_dir)
{

    std::vector<uint32_t> ids = braindead_r_extract_as_int<uint32_t>(analysis_tdf("Id"));
    std::vector<uint32_t> num_scans = braindead_r_extract_as_int<uint32_t>(analysis_tdf("NumScans"));
    std::vector<uint32_t> num_peaks = braindead_r_extract_as_int<uint32_t>(analysis_tdf("NumPeaks"));
    std::vector<uint32_t> msms_type = braindead_r_extract_as_int<uint32_t>(analysis_tdf("MsMsType"));
    Rcpp::NumericVector accum_time = analysis_tdf("AccumulationTime");
    Rcpp::NumericVector time = analysis_tdf("Time");
    std::vector<uint64_t> tims_id = braindead_r_extract_as_int<uint64_t>(analysis_tdf("TimsId"));

    for(size_t ii = 0; ii < ids.size(); ii++)
    {
        frame_descs.emplace(ids[ii], TimsFrame(
                ids[ii],
                num_scans[ii],
                num_peaks[ii],
                msms_type[ii],
                100.0 / accum_time[ii],
                time[ii],
                tims_id[ii] + tims_data_bin.data(),
                *this));
    }

    init();
}
#endif

TimsDataHandle::~TimsDataHandle()
{
    if(zstd_dctx != nullptr)
        ZSTD_freeDCtx(zstd_dctx);
#ifndef OPENTIMS_BUILDING_R
    if(db_conn != nullptr)
        sqlite3_close(db_conn);
#endif
    // std::cout << "KABOOM!!!" << std::endl; // JUST A TEST: this can be triggered by Python GC with pybind11.
}


TimsFrame& TimsDataHandle::get_frame(uint32_t frame_no)
{ 
    return frame_descs.at(frame_no); 
}

std::unordered_map<uint32_t, TimsFrame>& TimsDataHandle::get_frame_descs()
{
    return frame_descs;
}

size_t TimsDataHandle::no_peaks_in_frames(const uint32_t* indexes, size_t no_indexes)
{
    size_t ret = 0;
    for(size_t ii = 0; ii < no_indexes; ii++)
        ret += frame_descs.at(indexes[ii]).num_peaks;
    return ret;

}

size_t TimsDataHandle::no_peaks_in_slice(uint32_t start, uint32_t end, uint32_t step)
{
    size_t ret = 0;
    for(uint32_t ii = start; ii < end; ii += step)
        ret += frame_descs.at(ii).num_peaks;
    return ret;
}

size_t TimsDataHandle::no_peaks_total() const
{
    size_t ret = 0;
    for(auto it = frame_descs.begin(); it != frame_descs.end(); it++)
        ret += it->second.num_peaks;
    return ret;
}

void TimsDataHandle::set_converter(std::unique_ptr<Tof2MzConverter>&& converter)
{
    if(converter)
        tof2mz_converter = std::move(converter);
    else
        tof2mz_converter = DefaultTof2MzConverterFactory::produceDefaultConverterInstance(*this);
}

void TimsDataHandle::set_converter(std::unique_ptr<Scan2InvIonMobilityConverter>&& converter)
{
    if(converter)
        scan2inv_ion_mobility_converter = std::move(converter);
    else
        scan2inv_ion_mobility_converter = DefaultScan2InvIonMobilityConverterFactory::produceDefaultConverterInstance(*this);
}

void TimsDataHandle::extract_frames(const uint32_t* indexes,
                                    size_t no_indexes,
                                    uint32_t* result)
{
    size_t no_peaks = no_peaks_in_frames(indexes, no_indexes);

    uint32_t* offset0 = result;
    uint32_t* offset1 = offset0 + no_peaks;
    uint32_t* offset2 = offset1 + no_peaks;
    uint32_t* offset3 = offset2 + no_peaks;

    for(size_t ii = 0; ii < no_indexes; ii++)
    {
        TimsFrame& frame = frame_descs.at(indexes[ii]);
        frame.save_to_buffs(offset0, offset1, offset2, offset3, nullptr, nullptr, nullptr, zstd_dctx);
        offset0 += frame.num_peaks;
        offset1 += frame.num_peaks;
        offset2 += frame.num_peaks;
        offset3 += frame.num_peaks;
    }
}


void TimsDataHandle::extract_frames_slice(uint32_t start,
                                          uint32_t end,
                                          uint32_t step,
                                          uint32_t* result)
{
    size_t no_peaks = no_peaks_in_slice(start, end, step);

    uint32_t* offset0 = result;
    uint32_t* offset1 = offset0 + no_peaks;
    uint32_t* offset2 = offset1 + no_peaks;
    uint32_t* offset3 = offset2 + no_peaks;

    for(uint32_t ii = start; ii < end; ii += step)
    {
        TimsFrame& frame = frame_descs.at(ii);
        frame.save_to_buffs(offset0, offset1, offset2, offset3, nullptr, nullptr, nullptr, zstd_dctx);
        offset0 += frame.num_peaks;
        offset1 += frame.num_peaks;
        offset2 += frame.num_peaks;
        offset3 += frame.num_peaks;
    }
}

#define move_ptr(ptr) if(ptr) ptr += n;

void TimsDataHandle::extract_frames(const uint32_t* indexes,
                                    size_t no_indexes,
                                    uint32_t* frame_ids,
                                    uint32_t* scan_ids,
                                    uint32_t* tofs,
                                    uint32_t* intensities,
                                    double* mzs,
                                    double* inv_ion_mobilities,
                                    double* retention_times)
{
    for(size_t ii = 0; ii < no_indexes; ii++)
    {
        TimsFrame& frame = frame_descs.at(indexes[ii]);
        const size_t n = frame.num_peaks;
        frame_descs.at(indexes[ii]).save_to_buffs(frame_ids, scan_ids, tofs, intensities, mzs, inv_ion_mobilities, retention_times, zstd_dctx);
        move_ptr(frame_ids);
        move_ptr(scan_ids);
        move_ptr(tofs);
        move_ptr(intensities);
        move_ptr(mzs);
        move_ptr(inv_ion_mobilities);
        move_ptr(retention_times);
    }
}

void TimsDataHandle::extract_frames_slice(uint32_t start,
                                          uint32_t end,
                                          uint32_t step, 
                                          uint32_t* frame_ids,
                                          uint32_t* scan_ids,
                                          uint32_t* tofs,
                                          uint32_t* intensities,
                                          double* mzs,
                                          double* inv_ion_mobilities,
                                          double* retention_times)
{
    for(uint32_t ii = start; ii < end; ii += step)
    {
        TimsFrame& frame = frame_descs.at(ii);
        const size_t n = frame.num_peaks;
        frame_descs.at(ii).save_to_buffs(frame_ids, scan_ids, tofs, intensities, mzs, inv_ion_mobilities, retention_times, zstd_dctx);
        move_ptr(frame_ids);
        move_ptr(scan_ids);
        move_ptr(tofs);
        move_ptr(intensities);
        move_ptr(mzs);
        move_ptr(inv_ion_mobilities);
        move_ptr(retention_times);
    }
}



size_t TimsDataHandle::max_peaks_in_frame()
{
    size_t ret = 0;
    for(auto it = frame_descs.begin(); it != frame_descs.end(); it++)
        if(it->second.num_peaks > ret)
            ret = it->second.num_peaks;
    return ret;
}

void TimsDataHandle::allocate_buffers()
{
    size_t size = max_peaks_in_frame();
    _scan_ids_buffer = std::make_unique<uint32_t[]>(size);
    _tofs_buffer = std::make_unique<uint32_t[]>(size);
    _intensities_buffer = std::make_unique<uint32_t[]>(size);
}

void TimsDataHandle::free_buffers()
{
    _scan_ids_buffer = nullptr;
    _tofs_buffer = nullptr;
    _intensities_buffer = nullptr;
}

size_t TimsDataHandle::expose_frame(size_t frame_no)
{
    ensure_buffers_allocated();
    TimsFrame& frame = get_frame(frame_no);
    frame.save_to_buffs(nullptr, _scan_ids_buffer.get(), _tofs_buffer.get(), _intensities_buffer.get(), nullptr, nullptr, nullptr, zstd_dctx);
    return frame.num_peaks;
}

void TimsDataHandle::extract_frames(const std::vector<uint32_t>& indexes,
                                    uint32_t* const * frame_ids,
                                    uint32_t* const * scan_ids,
                                    uint32_t* const * tofs,
                                    uint32_t* const * intensities,
                                    double* const * mzs,
                                    double* const * inv_ion_mobilities,
                                    double* const * retention_times)
{
    std::atomic<size_t> current_task(0);

    ThreadingManager::get_instance().set_shared_threading();
    size_t n_threads = ThreadingManager::get_instance().get_no_opentims_threads();

    std::vector<std::thread> threads;
    for(size_t ii=0; ii<n_threads; ii++)
        threads.emplace_back([&](){
            std::unique_ptr<ZSTD_DCtx, decltype(&ZSTD_freeDCtx)> zstd(ZSTD_createDCtx(), &ZSTD_freeDCtx);
            std::unique_ptr<char[]> decomp_buffer = std::make_unique<char[]>(decomp_buffer_size);
            while(true)
            {
                size_t my_task = current_task.fetch_add(1);
                if(my_task < indexes.size())
                {
                    TimsFrame& frame = get_frame(indexes[my_task]);
                    frame.decompress(decomp_buffer.get(), zstd.get());
                    frame.save_to_buffs(frame_ids[my_task], scan_ids[my_task], tofs[my_task], intensities[my_task], mzs[my_task], inv_ion_mobilities[my_task], retention_times[my_task]);
                    frame.close();
                }
                else
                    break;
            }
        });
    for (auto& th : threads) th.join();
    ThreadingManager::get_instance().set_converter_threading();
}


void TimsDataHandle::per_frame_TIC(uint32_t* result)
{
    std::unique_ptr<uint32_t[]> intensities = std::make_unique<uint32_t[]>(max_peaks_in_frame());

    for(auto it = frame_descs.begin(); it != frame_descs.end(); it++)
    {
        it->second.save_to_buffs(nullptr, nullptr, nullptr, intensities.get(), nullptr, nullptr, nullptr, zstd_dctx);
        uint32_t acc = 0;
        const size_t n_peaks = it->second.num_peaks;
        for(size_t ii = 0; ii < n_peaks; ii++)
            acc += intensities[ii];
        result[it->first - 1] = acc;
    }
}
