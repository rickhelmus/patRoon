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

SpectrumRaw MSReadBackendOTIMS::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const
{
    auto &tframe = handle->get_frame(scan);
    if (tframe.num_peaks == 0)
        return SpectrumRaw();
    
    auto *bufs = reinterpret_cast<TIMSDecompBufferPair *>(tdata.get());
    
    tframe.decompress(bufs->second.get(), bufs->first.get());
    std::vector<uint32_t> IDs(tframe.num_peaks), intensities(tframe.num_peaks);
    std::vector<double> mzs(tframe.num_peaks), mobilities(tframe.num_peaks);
    tframe.save_to_buffs(nullptr, IDs.data(), nullptr, intensities.data(), mzs.data(), mobilities.data(), nullptr);
    tframe.close();
    
    SpectrumRaw ret(tframe.num_peaks, true);
    for(size_t i=0; i<ret.size(); ++i)
        ret.setPeak(i, mzs[i], intensities[i], mobilities[i]);
    
    return ret;
}
