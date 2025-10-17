/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <Rcpp.h>

#include "msdata-otims.h"
#include "spectrum-raw.h"

#ifdef WITH_OTIMS

#include "opentims++/opentims_all.h"

namespace {

auto getTIMSDecompBuffers(size_t size)
{
    // get buffers for decompression (see TimsDataHandle::extract_frames())
    return std::make_pair(std::unique_ptr<ZSTD_DCtx, decltype(&ZSTD_freeDCtx)>(ZSTD_createDCtx(), &ZSTD_freeDCtx),
                          std::make_unique<char[]>(size));
}

// From http://mochan.info/c++/2019/06/12/returntype-deduction.html
using TIMSDecompBufferPair = decltype(getTIMSDecompBuffers(std::declval<size_t>()));

struct OTIMSThreadData
{
    TIMSDecompBufferPair buffers;
    
    size_t curFrame;
    bool initialized = false;
    std::vector<uint32_t> IDs, intensities;
    std::vector<double> mzs, mobilities;
    
    OTIMSThreadData(TIMSDecompBufferPair &&b) : buffers(std::move(b)) { }
};

// utility class to load and use the TIMS SDK
class BrukerLibrary
{
    std::string path;
    bool loaded = false;
    // NOTE: we conveniently use LoadedLibraryHandle to handle loading the SDK library
    std::unique_ptr<LoadedLibraryHandle> libHandle = nullptr;

    // CCS <--> mobility conversion from TIMS-SDK
    typedef double tims_oneoverk0_to_ccs_for_mz_fun_t(const double ook0, const int charge, const double mz);
    tims_oneoverk0_to_ccs_for_mz_fun_t *tims_oneoverk0_to_ccs_for_mz = nullptr;
    typedef double tims_ccs_to_oneoverk0_for_mz_fun_t(const double ccs, const int charge, const double mz);
    tims_ccs_to_oneoverk0_for_mz_fun_t *tims_ccs_to_oneoverk0_for_mz = nullptr;
    
    void verifyIfLoaded(void) const
    {
        if (!loaded)
            throw(std::runtime_error("Cannot convert mobility/CCS data: Bruker SDK library not loaded!"));
    }
    
public:
    void load(const std::string &p)
    {
        loaded = false;
        
        setup_bruker(p);
        
        // NOTE: this disables threading from TIMS-SDK, which otherwise seems to cause issues with the threading applied
        // in patRoon
        ThreadingManager::get_instance().set_opentims_threading();
        
        libHandle = std::make_unique<LoadedLibraryHandle>(p);
        
        tims_oneoverk0_to_ccs_for_mz = libHandle->symbol_lookup<tims_oneoverk0_to_ccs_for_mz_fun_t>("tims_oneoverk0_to_ccs_for_mz");
        tims_ccs_to_oneoverk0_for_mz = libHandle->symbol_lookup<tims_ccs_to_oneoverk0_for_mz_fun_t>("tims_ccs_to_oneoverk0_for_mz");
        
        loaded = true; // No exception thrown if we are here, so should be fine
        path = p;
    }
    
    bool isLoaded(void) const { return loaded; }
    const std::string &currentPath(void) { return path; }
    
    double getCCS(double mob, int charge, double mz) const
    {
        verifyIfLoaded();
        return tims_oneoverk0_to_ccs_for_mz(mob, charge, mz);
    }
    double getMob(double ccs, int charge, double mz) const
    {
        verifyIfLoaded();
        return tims_ccs_to_oneoverk0_for_mz(ccs, charge, mz);
    }
    
};
BrukerLibrary brukerLibrary;

}

// [[Rcpp::interfaces(r, cpp)]]

void MSReadBackendOTIMS::doOpen(const std::string &file)
{
    handle = std::make_unique<TimsDataHandle>(file);
    
    if (handle->tof2mz_converter->description() != "BrukerTof2MzConverter" ||
        handle->scan2inv_ion_mobility_converter->description() != "BrukerScan2InvIonMobilityConverter")
        Rcpp::stop("Could not properly initialize OpenTIMS conversion functions. You may need to reload patRoon or restart R.");
}

void MSReadBackendOTIMS::doClose(void)
{
    handle.reset();
}

MSReadBackend::ThreadDataType MSReadBackendOTIMS::doGetThreadData(void) const
{
    return std::make_shared<OTIMSThreadData>(getTIMSDecompBuffers(handle->get_decomp_buffer_size()));
}

SpectrumRaw MSReadBackendOTIMS::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                               const SpectrumRawSelection &scanSel,
                                               const SpectrumRawTypes::MobilityRange &mobRange,
                                               SpectrumRawTypes::Intensity minIntensityIMS) const
{
    const auto &meta = getSpecMetadata();
    const auto &metaMS = (MSLevel == SpectrumRawTypes::MSLevel::MS1) ? meta.first : meta.second;

    // UNDONE: It seems that OpenTIMS and/or TIMS-SDK have some odd multithreading issues. It should work, but doesn't.
    // 
    // what doesn't work (consistently...!)
    // * trying to disable threading in TIMS-SDK via setting threads via OpenMP or OpenTIMS or OMP_NUM_THREADS
    // * using thread local handles and/or buffers
    // * save_to_buffs(): not converting mzs and mobilities (ie not using TIMS-SDK) seemed to sometimes imrove things,
    //   but inconsistently. Putting only that function in a critical section is not sufficient.
    // All or some of these issues are sometimes already triggered if called from a block with just one thread.
    // Naming the critical sections seem important too.
    //
    // so far the 2nd approach below seems best, ie still sharing handles/buffers but using them in a critical section.
    // But this needs more testing!
    
#if 0
    // can be used if not called from OpenMP threading block
    
    
    auto &tframe = handle->get_frame(metaMS.scans[scanSel.index]);
    
    std::vector<uint32_t> IDs(tframe.num_peaks), intensities(tframe.num_peaks);
    std::vector<double> mzs(tframe.num_peaks), mobilities(tframe.num_peaks);
    
    auto *bufs = reinterpret_cast<TIMSDecompBufferPair *>(tdata.get());
    tframe.decompress(bufs->second.get(), bufs->first.get());
    tframe.save_to_buffs(nullptr, IDs.data(), nullptr, intensities.data(), mzs.data(), mobilities.data(), nullptr);
    tframe.close();
#elif 1
    auto *otd = reinterpret_cast<OTIMSThreadData *>(tdata.get());
    if (!otd->initialized || otd->curFrame != scanSel.index)
    {
        #pragma omp critical (VeryCritical1)
        {
            auto &tframe = handle->get_frame(metaMS.scans[scanSel.index]);
            
            otd->IDs.assign(tframe.num_peaks, 0); otd->intensities.assign(tframe.num_peaks, 0.0);
            otd->mzs.assign(tframe.num_peaks, 0.0); otd->mobilities.assign(tframe.num_peaks, 0.0);
            
            auto *bufs = reinterpret_cast<TIMSDecompBufferPair *>(tdata.get());
            tframe.decompress(otd->buffers.second.get(), otd->buffers.first.get());
            tframe.save_to_buffs(nullptr, otd->IDs.data(), nullptr, otd->intensities.data(), otd->mzs.data(),
                                 otd->mobilities.data(), nullptr);
            tframe.close();
        }
        otd->curFrame = scanSel.index;
        otd->initialized = true;
    }

    
#else
    std::vector<uint32_t> IDs, intensities;
    std::vector<double> mzs, mobilities;
    
    #pragma omp critical (VeryCritical1)
    {
        //omp_set_num_threads(1);
        
        
        TimsDataHandle TDH(getCurrentFile());
        auto &tframe = TDH.get_frame(metaMS.scans[scanSel.index]);
        IDs.assign(tframe.num_peaks, 0); intensities.assign(tframe.num_peaks, 0.0);
        mzs.assign(tframe.num_peaks, 0.0); mobilities.assign(tframe.num_peaks, 0.0);
        
        //auto &tframe = handle->get_frame(metaMS.scans[scanSel.index]);
        /*if (tframe.num_peaks == 0)
            return SpectrumRaw();*/
        
        //auto *bufs = reinterpret_cast<TIMSDecompBufferPair *>(tdata.get());
        auto bufs = getTIMSDecompBuffers(TDH.get_decomp_buffer_size());
        
        tframe.decompress(bufs.second.get(), bufs.first.get());
        //#pragma omp critical (VeryCritical1)
        tframe.save_to_buffs(nullptr, IDs.data(), nullptr, intensities.data(), mzs.data(), mobilities.data(), nullptr);
        //tframe.save_to_buffs(nullptr, IDs.data(), nullptr, intensities.data(), nullptr, nullptr, nullptr);
        tframe.close();
    }
#endif
    
    SpectrumRaw ret;
    for (size_t i=0; i<otd->IDs.size(); ++i)
    {
        if (mobRange.isSet())
        {
            // NOTE: mobilities are sorted from high to low
            if (otd->mobilities[i] < mobRange.start)
                continue;
            if (otd->mobilities[i] > mobRange.end)
                break;
        }
        
        if (!scanSel.MSMSFrameIndices.empty())
        {
            bool inRange = false;
            const auto curScanID = otd->IDs[i];
            for (size_t j=0; j<scanSel.MSMSFrameIndices.size() && !inRange; ++j)
            {
                const auto fri = scanSel.MSMSFrameIndices[j];
                const auto ss = meta.second.MSMSFrames[scanSel.index].subScans[fri];
                const auto se = meta.second.MSMSFrames[scanSel.index].subScanEnds[fri];
                // NOTE: scan ranges are zero for isCID/bbCID
                // NOTE: scan ends are exclusive
                inRange = (ss == 0 && se == 0) || (curScanID >= ss && curScanID < se);
            }
            if (!inRange)
            {
                // try again with next scan: increment until last element (i will be incremented again in main for loop)
                for (; i < (otd->IDs.size()-1) && otd->IDs[i+1] == curScanID; ++i)
                    ;
                continue;
            }
        }
        
        if (otd->intensities[i] >= minIntensityIMS)
            ret.append(otd->mzs[i], otd->intensities[i], otd->mobilities[i]);
    }
    
    return ret;
}

#endif // WITH_OTIMS

// [[Rcpp::export]]
bool initBrukerLibrary(const std::string &path, bool force = false)
{
#ifdef WITH_OTIMS
    const std::string lastPath = brukerLibrary.currentPath();
    
    if (!force && brukerLibrary.isLoaded() && path == brukerLibrary.currentPath())
        return true; // already loaded
    
    try
    {
        brukerLibrary.load(path);
    }
    catch(const std::exception &e)
    {
        if (lastPath != path) // only warn the first time
            Rcpp::warning("Failed to load Bruker TIMS library ('%s'): %s", path.c_str(), e.what());
    }
    
    return brukerLibrary.isLoaded();
#else
    return false;
#endif
}

// [[Rcpp::export]]
Rcpp::NumericVector getBrukerCCS(Rcpp::NumericVector mobs, Rcpp::IntegerVector charges, Rcpp::NumericVector mzs)
{
#ifdef WITH_OTIMS
    Rcpp::NumericVector ret(mobs.size());
    
    for (int i=0; i<ret.size(); ++i)
    {
        if (mobs[i] == NA_REAL || charges[i] == NA_INTEGER || mzs[i] == NA_REAL)
            ret[i] = NA_REAL;
        else
            ret[i] = brukerLibrary.getCCS(mobs[i], charges[i], mzs[i]);
    }
    
    return ret;
#else
    Rcpp::stop("Cannot calculate CCS values with Bruker library: backend unavailable!");
    return Rcpp::NumericVector();
#endif
}

// [[Rcpp::export]]
Rcpp::NumericVector getBrukerMob(Rcpp::NumericVector ccss, Rcpp::IntegerVector charges, Rcpp::NumericVector mzs)
{
#ifdef WITH_OTIMS
    Rcpp::NumericVector ret(ccss.size());
    
    for (int i=0; i<ret.size(); ++i)
    {
        if (ccss[i] == NA_REAL || charges[i] == NA_INTEGER || mzs[i] == NA_REAL)
            ret[i] = NA_REAL;
        else
            ret[i] = brukerLibrary.getMob(ccss[i], charges[i], mzs[i]);
    }
    
    return ret;
#else
    Rcpp::stop("Cannot calculate mobility values with Bruker library: backend unavailable!");
    return Rcpp::NumericVector();
#endif
}
