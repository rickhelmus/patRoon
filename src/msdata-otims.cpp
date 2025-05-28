#include <Rcpp.h>

#include "msdata-otims.h"
#include "spectrum-raw.h"

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
    auto &tframe = handle->get_frame(metaMS.scans[scanSel.index]);
    if (tframe.num_peaks == 0)
        return SpectrumRaw();
    
    auto *bufs = reinterpret_cast<TIMSDecompBufferPair *>(tdata.get());
    
    tframe.decompress(bufs->second.get(), bufs->first.get());
    std::vector<uint32_t> IDs(tframe.num_peaks), intensities(tframe.num_peaks);
    std::vector<double> mzs(tframe.num_peaks), mobilities(tframe.num_peaks);
    tframe.save_to_buffs(nullptr, IDs.data(), nullptr, intensities.data(), mzs.data(), mobilities.data(), nullptr);
    tframe.close();

    SpectrumRaw ret;
    for (size_t i=0; i<tframe.num_peaks; ++i)
    {
        if (mobRange.isSet())
        {
            // NOTE: mobilities are sorted  from high to low
            if (mobilities[i] > mobRange.start)
                continue;
            if (mobilities[i] < mobRange.end)
                break;
        }
        
        if (!scanSel.MSMSFrameIndices.empty())
        {
            bool inRange = false;
            const auto curScanID = IDs[i];
            for (size_t j=0; j<scanSel.MSMSFrameIndices.size() && !inRange; ++j)
            {
                inRange = (curScanID >= meta.second.MSMSFrames[scanSel.index].subScans[j] &&
                    curScanID <= meta.second.MSMSFrames[scanSel.index].subScanEnds[j]);
            }
            if (!inRange)
            {
                // try again with next scan: increment until last element (i will be incremented again in main for loop)
                for (; i < (tframe.num_peaks-1) && IDs[i+1] == curScanID; ++i)
                    ;
                continue;
            }
        }
        ret.append(mzs[i], intensities[i], mobilities[i]);
    }
    
    return ret;
}

// [[Rcpp::export]]
bool initBrukerLibrary(const std::string &path)
{
    static std::string lastPath;
    static bool succeeded = false;
    
    if (!lastPath.empty() && succeeded && path == lastPath)
        return true; // already loaded
    
    try
    {
        setup_bruker(path);
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
}
