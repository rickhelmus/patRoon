#include <Rcpp.h>

#include "msdata-sc.h"
#include "spectrum-raw.h"

#define PUGIXML_PATH "../pugixml/pugixml.hpp"
#define STREAMCRAFT_HEADER_ONLY
#include "StreamCraft/StreamCraft_lib.hpp"
#undef PUGIXML_PATH

namespace {

SpectrumRaw getSCSpectrum(sc::MS_ANALYSIS *analysis, SpectrumRawTypes::Scan scan, const SpectrumRawTypes::MobilityRange &mobRange)
{
    const auto s = analysis->get_spectrum(scan);
    
    const auto mobArrayIt = std::find(s.binary_names.cbegin(), s.binary_names.cend(), "ion_mobility");
    const bool hasMob = mobArrayIt != s.binary_names.cend();
    
    SpectrumRaw ret(s.array_length, hasMob);
    
    // UNDONE: always safe to assume that m/z / intensity are in first two arrays?
    
    if (hasMob)
    {
        const size_t mobArrayInd = std::distance(s.binary_names.begin(), mobArrayIt);
        for (int i=0; i<s.array_length; ++i)
        {
            const auto mob = s.binary_data[mobArrayInd][i];
            if (!mobRange.isSet() || mobRange.within(mob))
                ret.setPeak(i, s.binary_data[0][i], s.binary_data[1][i], mob);
        }
    }
    else
    {
        for (int i=0; i<s.array_length; ++i)
            ret.setPeak(i, s.binary_data[0][i], s.binary_data[1][i]);
    }
    
    return ret;
}

}

// [[Rcpp::interfaces(r, cpp)]]

MSReadBackend::ThreadDataType MSReadBackendSC::doGetThreadData(void) const
{
    return std::make_shared<sc::MS_ANALYSIS>(getCurrentFile());
}

SpectrumRaw MSReadBackendSC::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                            const SpectrumRawSelection &scanSel,
                                            const SpectrumRawTypes::MobilityRange &mobRange) const
{
    auto *analysis = reinterpret_cast<sc::MS_ANALYSIS *>(tdata.get());
    const auto &meta = getSpecMetadata();
    
    if (MSLevel == SpectrumRawTypes::MSLevel::MS1)
        return getSCSpectrum(analysis, meta.first.scans[scanSel.index], mobRange);
    if (scanSel.MSMSFrameIndices.empty())
        return getSCSpectrum(analysis, meta.second.scans[scanSel.index], mobRange);
    
    // if we are here we need to get MS2 data from an IMS frame...
    
    SpectrumRaw ret;
    for (const auto i : scanSel.MSMSFrameIndices)
        ret.append(getSCSpectrum(analysis, meta.second.MSMSFrames[scanSel.index].subScans[i], mobRange));
    
    return ret;
}

void MSReadBackendSC::generateSpecMetadata(void)
{
    if (getCurrentFile().empty())
        return;

    sc::MS_ANALYSIS analysis(getCurrentFile());
    const auto hd = analysis.get_spectra_headers();
    
    const bool hasMob = analysis.has_ion_mobility();
    
    SpectrumRawMetadata meta;
    for (size_t i=0; i<hd.index.size(); ++i)
    {
        const bool isMS1 = hd.level[i] == 1;
        
        SpectrumRawMetadataMS *curMS1MD = (isMS1) ? &meta.first : &meta.second;
        curMS1MD->scans.push_back(hd.index[i]);
        curMS1MD->times.push_back(hd.rt[i]);
        curMS1MD->TICs.push_back(hd.tic[i]);
        curMS1MD->BPCs.push_back(hd.bpint[i]);
        // NOTE: a simple cast is sufficient as the enum values match the int values from SC
        curMS1MD->polarities.push_back(static_cast<SpectrumRawTypes::MSPolarity>(hd.polarity[i]));
        
        if (!isMS1)
        {
            if (!hasMob)
            {
                const SpectrumRawTypes::Mass l = hd.precursor_mz[i] - hd.window_mzlow[i], h = hd.precursor_mz[i] + hd.window_mzhigh[i];
                meta.second.isolationRanges.emplace_back(l, h);
            }
            else
            {
                // For IMS-MS/MS data, the different spectra inside a frame are stored in separate spectra
                // --> move these spectra to MSMSFrames structs, so our main table only contains separate frames
                frameMSMSInfo fi;
                const auto curTime = hd.rt[i];
                bool finished = false;
                while (!finished)
                {
                    const SpectrumRawTypes::Mass l = hd.precursor_mz[i] - hd.window_mzlow[i], h = hd.precursor_mz[i] + hd.window_mzhigh[i];
                    fi.isolationRanges.emplace_back(l, h);
                    fi.subScans.push_back(hd.index[i]);
                    ++i;
                    if (i >= hd.index.size())
                        break;
                    if (!compareTol(curTime, hd.rt[i]))
                        finished = true; // reached next frame
                }
                meta.second.MSMSFrames.push_back(std::move(fi));
                if (finished)
                    --i; // we reached one spec to far to find out we're finished
            }
        }
    }
    
    setSpecMetadata(std::move(meta));
}
