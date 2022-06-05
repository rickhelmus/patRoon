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

#pragma once
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "platform.h"

#ifndef OPENTIMS_BUILDING_R
#include "sqlite/sqlite3.h"
#endif

#include "zstd/zstd.h"


#ifdef OPENTIMS_BUILDING_R
#define STRICT_R_HEADERS
#include "mio.h"
#include <Rcpp.h>
#else
#include "mio.hpp"
#endif



class TimsDataHandle;

int tims_sql_callback(void* out, int cols, char** row, char** colnames);

class TimsFrame
{
    std::unique_ptr<char[]> back_buffer;

    char* bytes0;
    char* bytes1;
    char* bytes2;
    char* bytes3;

    inline uint32_t back_data(size_t index)
    {
        uint32_t ret;
        char* bytes = reinterpret_cast<char*>(&ret);

        bytes[0] = bytes0[index];
        bytes[1] = bytes1[index];
        bytes[2] = bytes2[index];
        bytes[3] = bytes3[index];

        return ret;
    }

    const char * const tims_bin_frame;

    friend class TimsDataHandle;
    friend int tims_sql_callback(void* out, int cols, char** row, char** colnames);

    TimsDataHandle& parent_tdh;

    /* TODO: implement the below one day.
    template void save_to_buffs_impl<bool frame_ids_present,
                                     bool scan_ids_present,
                                     bool intensities_present,
                                     bool mzs_present,
                                     bool inv_ion_mobilities_present,
                                     bool retention_times_present,
                                     > (uint32_t* frame_ids,
                                        uint32_t* scan_ids,
                                        uint32_t* tofs,
                                        uint32_t* intensities,
                                        double* mzs,
                                        double* inv_ion_mobilities,
                                        double* retention_times,
                                        ZSTD_DCtx* decomp_ctx = nullptr);
    */

    TimsFrame(uint32_t _id,
              uint32_t _num_scans,
              uint32_t _num_peaks,
              uint32_t _msms_type,
              double _intensity_correction,
              double _time,
              const char* frame_ptr,
              TimsDataHandle& parent_hndl
            );

    static TimsFrame TimsFrameFromSql(char** sql_row, 
                                      TimsDataHandle& parent_handle);

    inline size_t data_size_ints() const { return num_scans + num_peaks + num_peaks; };


public:

    const uint32_t id;          ///< ID of the frame
    const uint32_t num_scans;   ///< Number of scans this frame contains
    const uint32_t num_peaks;   ///< Number of peaks this frame contains (summed across all scans)
    const uint32_t msms_type;   ///< The MS/MS type of this frame
    const double intensity_correction;
    const double time;

     //! \brief Prints out to stdout a short summary of this frame.
    void print() const;

    //! Return the size of back buffer needed to store raw TIMS data.
    inline size_t data_size_bytes() const { return data_size_ints() * 4; };

    //! Precalculates and memorizes all the information contained within this frame.
    /**
     * This method extracts all the data associated with the frame and stores it in memory.
     * Without it, the data (mz values, intensities, etc.) will be re-calculated every time they
     * are accessed. After this method is called, they will be just retrieved from RAM.
     *
     * This is a purely performance-enchancing method - it is recommended to call it if
     * repeated access to the frame is necessary, before the access, and to call close()
     * afterward.
     *
     * @param decompression_buffer optional, a pre-allocated buffer which will be used for 
     *        the decompression.
     *        If null, a pre-allocated buffer from parent TimsDataFrame will be used,
     *        making this method non-thread-safe. If TimsFrame is to be used in multithreaded
     *        context, a thread-local buffer must be passed here (or a mutex must be used).
     *        The buffer must be at least data_size_bytes() large.
     * @param decomp_ctx optonal, a decompression buffer to be used by the method. If null,
     *        a ctx from parent handle will be used, possibly making this non-thread-safe.
     *        To ensure thread safety a thread-local context must be passed here.
     */
    void decompress(char* decompression_buffer = nullptr, ZSTD_DCtx* decomp_ctx = nullptr);

    //! Releases the storage taken by decompress() method, without destroying the frame. Will be called in destructor.
    void close();

    //! Retrieve the MS peak data held by the frame.
    /**
     * The data is saved to the passed buffers - if some of this data is unnecessary, then nullptr
     * may be passed as the corresponding pointer.
     * The buffers are passed and returned by columns. Each row corresponds to one MS peak.
     * Each buffer must be albe to hold at least this->num_peaks values.
     *
     * @param frame_ids     The repeated ID of this frame.
     * @param scan_ids      IDs of the scan a peak comes from.
     * @param tofs          Times of Flight of peaks.
     * @param intensities   Signal intensities of peaks.
     * @param mzs           M/Z ratios of peaks.
     * @param inv_ion_mobilities    Inverse ion mobilities (in seconds).
     * @param retention_times       Retention times (in seconds).
     */
    void save_to_buffs(uint32_t* frame_ids,
                       uint32_t* scan_ids,
                       uint32_t* tofs,
                       uint32_t* intensities,
                       double* mzs,
                       double* inv_ion_mobilities,
                       double* retention_times,
                       ZSTD_DCtx* decomp_ctx = nullptr);

    //! This function is deprecated and intentionally undocumented; do not use.
    void save_to_matrix_buffer(uint32_t* buf,
                               ZSTD_DCtx* decomp_ctx = nullptr)
    { save_to_buffs(buf, buf+num_peaks, buf+2*num_peaks, buf+3*num_peaks, nullptr, nullptr, nullptr, decomp_ctx); };

    friend class TimsDataHandle;
};

class BrukerTof2MzConverter;
class Tof2MzConverter;
class BrukerScan2InvIonMobilityConverter;
class Scan2InvIonMobilityConverter;

class TimsDataHandle
{
friend class BrukerTof2MzConverter;
friend class BrukerScan2InvIonMobilityConverter;

private:
    const std::string tims_dir_path;
    mio::mmap_source tims_data_bin;
    std::unordered_map<uint32_t, TimsFrame> frame_descs;
    void read_sql(const std::string& tims_tdf_path);
    uint32_t _min_frame_id;
    uint32_t _max_frame_id;

    std::unique_ptr<char[]> decompression_buffer;
    size_t decomp_buffer_size;

    std::unique_ptr<uint32_t[]> _scan_ids_buffer;
    std::unique_ptr<uint32_t[]> _tofs_buffer;
    std::unique_ptr<uint32_t[]> _intensities_buffer;

    ZSTD_DCtx* zstd_dctx;

#ifndef OPENTIMS_BUILDING_R
    sqlite3* db_conn;
#endif

public:
    size_t get_decomp_buffer_size() const { return decomp_buffer_size; };
    std::unique_ptr<Tof2MzConverter> tof2mz_converter;
    std::unique_ptr<Scan2InvIonMobilityConverter> scan2inv_ion_mobility_converter;

private:
    void init();

#ifdef OPENTIMS_BUILDING_R
    void* setupFromAnalysisList(const Rcpp::List& analysis_tdf);
#endif /* OPENTIMS_BUILDING_R */

    TimsDataHandle(const std::string& tims_tdf_bin_path,
                   const std::string& tims_tdf_path,
                   const std::string& tims_data_dir);

    void set_converter(std::unique_ptr<Tof2MzConverter>&& converter);
    void set_converter(std::unique_ptr<Scan2InvIonMobilityConverter>&& converter);

public:

    //! Open TimsTOF dataset.
    /**
     * Open a TimsTOF dataset (read-only access supported). Return a handle using which
     * Tims data may be retrieved.
     *
     * @param Path to the dataset (typically a *.d directory, containing a analysis.tdf
     * and analysis.tdf_bin files).
     */
    TimsDataHandle(const std::string& tims_data_dir);

#ifdef OPENTIMS_BUILDING_R
    //! Internal use only.
    TimsDataHandle(const std::string& tims_data_dir, const Rcpp::List& analysis_tdf);
#endif /* OPENTIMS_BUILDING_R */

    //! Close and deallocate the TimsTOF data handle (destructor).
    ~TimsDataHandle();

    //! Access a single frame by its ID.
    TimsFrame& get_frame(uint32_t frame_no);

    //! Access a dictionary containing all the frames from this dataset, keyed by ID.
    std::unordered_map<uint32_t, TimsFrame>& get_frame_descs();

    //! Returns the total number of MS peaks in this handle.
    size_t no_peaks_total() const;

    //! Count the peaks in a subset of frames
    /**
     * Returns the total number of peaks in the frames with given indexes.
     *
     * @param indexes    Indexes of frames to count
     * @param no_indexes Number of indexes (and length of the indexes[] table).
     */
    size_t no_peaks_in_frames(const uint32_t indexes[],
                              size_t no_indexes);

    //! Count the peaks in a subset of frames
    /**
     * Returns the total number of peaks in the frames with given indexes.
     *
     * @param indexes    Indexes of frames to count
     */
    size_t no_peaks_in_frames(const std::vector<uint32_t>& indexes)
    { return no_peaks_in_frames(indexes.data(), indexes.size()); };

    //! Count the peaks in a subset of frames, selected by a slice
    /**
     * Returns the total number of peaks in frames with IDs contained
     * in the start:stop:end slice
     */
    size_t no_peaks_in_slice(uint32_t start,
                             uint32_t end,
                             uint32_t step);

    //! Access the lowest id of a valid frame from this dataset.
    uint32_t min_frame_id() const { return _min_frame_id; };

    //! Access the highest id of a valid frame from this dataset.
    uint32_t max_frame_id() const { return _max_frame_id; };

    //! Check whether a frame with provided ID exists in the dataset.
    bool has_frame(uint32_t frame_id) const { return frame_descs.count(frame_id) > 0; };

    //! This function is deprecated, and left deliberately undocumented; do not use.
    void extract_frames(const uint32_t* indexes,
                        size_t no_indexes,
                        uint32_t* result);

    //! This function is deprecated, and left deliberately undocumented; do not use.
    void extract_frames(const std::vector<uint32_t>& indexes,
                        uint32_t* result)
    { extract_frames(indexes.data(), indexes.size(), result); };

    //! This function is deprecated, and left deliberately undocumented; do not use.
    void extract_frames_slice(uint32_t start,
                              uint32_t end,
                              uint32_t step,
                              uint32_t* result);

    //! This function is deprecated, and left deliberately undocumented; do not use.
    void extract_frames(const uint32_t* indexes,
                        size_t no_indexes,
                        uint32_t* frame_ids,
                        uint32_t* scan_ids,
                        uint32_t* tofs,
                        uint32_t* intensities,
                        double* mzs,
                        double* inv_ion_mobilities,
                        double* retention_times);

    //! Extract a subset of frames, selected by indexes, filling provided buffers with MS peak data.
    /**
     * The data is saved to the passed buffers - if some of this data is unnecessary, then nullptr
     * may be passed as the corresponding pointer.
     * The buffers are passed and returned by columns. Each row corresponds to one MS peak.
     * Each buffer must be able to hold at least no_peaks_in_frames(indexes) values.
     *
     * @param indexes       Set of indexes of frames for which data is to be obtained.
     * @param frame_ids     The IDs of frames containing the associated peaks.
     * @param scan_ids      IDs of the scan a peak comes from.
     * @param tofs          Times of Flight of peaks.
     * @param intensities   Signal intensities of peaks.
     * @param mzs           M/Z ratios of peaks.
     * @param inv_ion_mobilities    Inverse ion mobilities (in seconds).
     * @param retention_times       Retention times (in seconds).
     */
    void extract_frames(const std::vector<uint32_t>& indexes,
                        uint32_t* frame_ids,
                        uint32_t* scan_ids,
                        uint32_t* tofs,
                        uint32_t* intensities,
                        double* mzs,
                        double* inv_ion_mobilities,
                        double* retention_times)
    { extract_frames(indexes.data(),
                     indexes.size(),
                     frame_ids,
                     scan_ids,
                     tofs,
                     intensities,
                     mzs,
                     inv_ion_mobilities,
                     retention_times); };

    //! Extract a subset of frames, selected by a slice, filling provided buffers with MS peak data.
    /**
     * The data is saved to the passed buffers - if some of this data is unnecessary, then nullptr
     * may be passed as the corresponding pointer.
     * The buffers are passed and returned by columns. Each row corresponds to one MS peak.
     * Each buffer must be able to hold at least no_peaks_in_slice(start, stop, end) values.
     *
     * IDs of the returned frames come from start:stop:step slice.
     *
     * @param indexes       Start of the slice.
     * @param end           End of the slice.
     * @param step          Step of the slice.
     * @param frame_ids     The IDs of frames containing the associated peaks.
     * @param scan_ids      IDs of the scan a peak comes from.
     * @param tofs          Times of Flight of peaks.
     * @param intensities   Signal intensities of peaks.
     * @param mzs           M/Z ratios of peaks.
     * @param inv_ion_mobilities    Inverse ion mobilities (in seconds).
     * @param retention_times       Retention times (in seconds).
     */
    void extract_frames_slice(uint32_t start,
                              uint32_t end, 
                              uint32_t step,
                              uint32_t* frame_ids,
                              uint32_t* scan_ids,
                              uint32_t* tofs,
                              uint32_t* intensities,
                              double* mzs,
                              double* inv_ion_mobilities,
                              double* retention_times);

    //! Extract a subset of frames, selected by indexes, filling provided buffers with MS peak data.
    /**
     * The data is saved to the passed arrays of buffers - if some of this data is unnecessary, then
     * nullptr-filled array may be passed as the corresponding argument.
     *
     * Each array (for example scan_ids array) must hold indexes.size() pointers.
     * Each of these pointers must be a nullptr, or must point to a block of memory
     * able to hold enough values: for example array pointed to by tofs[13] must be able
     * to hold this->get_frame(indexes[13]).num_peaks values.
     *
     * @param indexes       Set of indexes of frames for which data is to be obtained.
     * @param frame_ids     The IDs of frames containing the associated peaks.
     * @param scan_ids      IDs of the scan a peak comes from.
     * @param tofs          Times of Flight of peaks.
     * @param intensities   Signal intensities of peaks.
     * @param mzs           M/Z ratios of peaks.
     * @param inv_ion_mobilities    Inverse ion mobilities (in seconds).
     * @param retention_times       Retention times (in seconds).
     */
    void extract_frames(const std::vector<uint32_t>& indexes,
                        uint32_t* const * frame_ids,
                        uint32_t* const * scan_ids,
                        uint32_t* const * tofs,
                        uint32_t* const * intensities,
                        double* const * mzs,
                        double* const * inv_ion_mobilities,
                        double* const * retention_times);

    void allocate_buffers();

    inline void ensure_buffers_allocated() { if(_scan_ids_buffer) return; allocate_buffers(); };

    void free_buffers();

    //! Return the maximal number of peaks in the biggest frame in this dataset.
    size_t max_peaks_in_frame();

    //! Expermental API - use discouraged.
    size_t expose_frame(size_t frame_id);

    //! Expermental API - use discouraged.
    const std::unique_ptr<uint32_t[]>& scan_ids_buffer() { return _scan_ids_buffer; };

    //! Expermental API - use discouraged.
    const std::unique_ptr<uint32_t[]>& tofs_buffer() { return _tofs_buffer; };

    //! Expermental API - use discouraged.
    const std::unique_ptr<uint32_t[]>& intensities_buffer() { return _intensities_buffer; };

//    const sqlite3* db_connection() { return db_conn; };

    //! Obtain the Total Ionic Current for each frame present in the spectrum
    /** The data is saved to the argument buffer - which must be able to hold at least
     * max_frame_id()-1 values. The number at nth index corresponds to n+1st frame (as
     * frames are numbered starting at 1).
     */
    void per_frame_TIC(uint32_t* result);

    friend int tims_sql_callback(void* out, int cols, char** row, char** colnames);

    friend class TimsFrame;
};
