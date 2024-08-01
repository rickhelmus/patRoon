#include <Rcpp.h>

#include "msdata-sc.h"
#include "spectrum-raw.h"

#define PUGIXML_PATH "../pugixml/pugixml.hpp"
#include "StreamCraft/StreamCraft_mzml.cpp"
#include "StreamCraft/StreamCraft_utils.cpp"
#undef PUGIXML_PATH

namespace {

SpectrumRaw getSCSpectrum(sc::MZML *mzml, SpectrumRawTypes::Scan scan,
                          const SpectrumRawTypes::MobilityRange &mobRange)
{
    const auto s = mzml->get_spectrum(scan);
    
    SpectrumRaw ret(s.array_length);
    // UNDONE: handle IMS data/mobRange
    for(int i=0; i<s.array_length; ++i)
        ret.setPeak(i, s.binary_data[0][i], s.binary_data[1][i]); // UNDONE: OK to assume it's always this order?
    return ret;
}

}

// [[Rcpp::interfaces(r, cpp)]]

MSReadBackend::ThreadDataType MSReadBackendSC::doGetThreadData(void) const
{
    return std::make_shared<sc::MZML>(getCurrentFile());
}

SpectrumRaw MSReadBackendSC::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                            const SpectrumRawSelection &scanSel,
                                            const SpectrumRawTypes::MobilityRange &mobRange) const
{
    auto *mzml = reinterpret_cast<sc::MZML *>(tdata.get());
    const auto &meta = getSpecMetadata();
    
    if (MSLevel == SpectrumRawTypes::MSLevel::MS1)
        return getSCSpectrum(mzml, meta.first.scans[scanSel.index], mobRange);
    if (scanSel.MSMSFrameIndices.empty())
        return getSCSpectrum(mzml, meta.second.scans[scanSel.index], mobRange);
    
    // if we are here we need to get MS2 data from an IMS frame...
    
    SpectrumRaw ret;
    for (const auto i : scanSel.MSMSFrameIndices)
        ret.append(getSCSpectrum(mzml, meta.second.MSMSFrames[scanSel.index].subScans[i], mobRange));
    
    return ret;
}

void MSReadBackendSC::generateSpecMetadata(void)
{
    if (getCurrentFile().empty())
        return;

    sc::MZML mzml(getCurrentFile());
    const auto hd = mzml.get_spectra_headers();
    
    SpectrumRawMetadata meta;
    for (size_t i=0; i<hd.index.size(); ++i)
    {
        const bool isMS1 = hd.level[i] == 1;
        
        SpectrumRawMetadataMS *curMS1MD = (isMS1) ? &meta.first : &meta.second;
        curMS1MD->scans.push_back(hd.index[i]);
        curMS1MD->times.push_back(hd.rt[i]);
        curMS1MD->TICs.push_back(hd.tic[i]);
        curMS1MD->BPCs.push_back(hd.bpint[i]);
        
        if (!isMS1)
        {
            const SpectrumRawTypes::Mass l = hd.precursor_mz[i] - hd.window_mzlow[i], h = hd.precursor_mz[i] + hd.window_mzhigh[i];
            meta.second.isolationRanges.emplace_back(l, h);
        }
    }
    
    setSpecMetadata(std::move(meta));
}
