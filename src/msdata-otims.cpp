
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
    
}

// [[Rcpp::interfaces(r, cpp)]]

void MSReadBackendOTIMS::doOpen(const std::string &file)
{
    handle = std::make_unique<TimsDataHandle>(file);
}

void MSReadBackendOTIMS::doClose(void)
{
    handle.reset();
}

MSReadBackend::ThreadDataType MSReadBackendOTIMS::doGetThreadData(void) const
{
    return std::make_shared<TIMSDecompBufferPair>(getTIMSDecompBuffers(handle->get_decomp_buffer_size()));
}

SpectrumRaw MSReadBackendOTIMS::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                               const SpectrumRawSelection &scanSel,
                                               const SpectrumRawTypes::MobilityRange &mobRange) const
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
    std::vector<uint32_t> IDs, intensities;
    std::vector<double> mzs, mobilities;

    #pragma omp critical (VeryCritical1)
    {
        auto &tframe = handle->get_frame(metaMS.scans[scanSel.index]);
        
        IDs.assign(tframe.num_peaks, 0); intensities.assign(tframe.num_peaks, 0.0);
        mzs.assign(tframe.num_peaks, 0.0); mobilities.assign(tframe.num_peaks, 0.0);
        
        auto *bufs = reinterpret_cast<TIMSDecompBufferPair *>(tdata.get());
        tframe.decompress(bufs->second.get(), bufs->first.get());
        tframe.save_to_buffs(nullptr, IDs.data(), nullptr, intensities.data(), mzs.data(), mobilities.data(), nullptr);
        tframe.close();
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
    for (size_t i=0; i<IDs.size(); ++i)
    {
        if (mobRange.isSet())
        {
            // NOTE: mobilities are sorted from high to low
            if (mobilities[i] < mobRange.start)
                continue;
            if (mobilities[i] > mobRange.end)
                break;
        }
        
        if (!scanSel.MSMSFrameIndices.empty())
        {
            bool inRange = false;
            const auto curScanID = IDs[i];
            for (size_t j=0; j<scanSel.MSMSFrameIndices.size() && !inRange; ++j)
            {
                const auto ss = meta.second.MSMSFrames[scanSel.index].subScans[j];
                const auto se = meta.second.MSMSFrames[scanSel.index].subScanEnds[j];
                // NOTE: scan ranges are zero for isCID/bbCID
                // NOTE: scan ends are exclusive
                inRange = (ss == 0 && se == 0) || (curScanID >= ss && curScanID < se);
            }
            if (!inRange)
            {
                // try again with next scan: increment until last element (i will be incremented again in main for loop)
                for (; i < (IDs.size()-1) && IDs[i+1] == curScanID; ++i)
                    ;
                continue;
            }
        }
        ret.append(mzs[i], intensities[i], mobilities[i]);
    }
    
    return ret;
}

#endif // WITH_OTIMS

// [[Rcpp::export]]
bool initBrukerLibrary(const std::string &path, bool force = false)
{
#ifdef WITH_OTIMS
    static std::string lastPath;
    static bool succeeded = false;
    
    if (!force && !lastPath.empty() && succeeded && path == lastPath)
        return true; // already loaded
    
    try
    {
        setup_bruker(path);
        // NOTE: this disables threading from TIMS-SDK, which otherwise seems to cause issues with the threading applied
        // in patRoon
        ThreadingManager::get_instance().set_opentims_threading();
        succeeded = true;
    }
    catch(const std::runtime_error &e)
    {
        if (lastPath != path) // only warn the first time
            Rcpp::warning("Failed to load Bruker TIMS library ('%s'): %s", path.c_str(), e.what());
        succeeded = false;
    }
    
    lastPath = path;
    
    return succeeded;
#else
    return false;
#endif
}
