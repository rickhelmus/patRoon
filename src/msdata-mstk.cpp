#include <Rcpp.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "msdata-mstk.h"
#include "spectrum-raw.h"

namespace {

MSToolkit::MSReader getMSTKReader(void)
{
    MSToolkit::MSReader ret;
    ret.addFilter(MSToolkit::MS1);
    ret.addFilter(MSToolkit::MS2);
    return ret;
}

SpectrumRaw getMSTKSpec(MSToolkit::MSReader *msr, const std::string &file, SpectrumRawTypes::Scan scan,
                        const SpectrumRawTypes::MobilityRange &mobRange)
{
    MSToolkit::Spectrum s;
    
    if (!msr->readFile(file.c_str(), s, scan) || s.getScanNumber() == 0)
        Rcpp::stop("Abort: invalid spectrum scan index: %d", scan);
    
    const bool hasMob = s.hasIonMobilityArray();
    SpectrumRaw ret(s.size(), hasMob);
    
    if (hasMob)
    {
        for(size_t i=0; i<ret.size(); ++i)
        {
            if (!mobRange.isSet() || mobRange.within(s.atIM(i)))
                ret.setPeak(i, s.at(i).mz, s.at(i).intensity, s.atIM(i));
        }
    }
    else
    {
        for(size_t i=0; i<ret.size(); ++i)
            ret.setPeak(i, s.at(i).mz, s.at(i).intensity);
    }
    
    return ret;
}

}


int MSReadBackendMSTK::backends = 0;

// [[Rcpp::interfaces(r, cpp)]]

MSReadBackend::ThreadDataType MSReadBackendMSTK::doGetThreadData(void) const
{
    return std::make_shared<MSToolkit::MSReader>(getMSTKReader());
}

SpectrumRaw MSReadBackendMSTK::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                              const SpectrumRawSelection &scanSel,
                                              const SpectrumRawTypes::MobilityRange &mobRange) const
{
    auto *msr = reinterpret_cast<MSToolkit::MSReader *>(tdata.get());
    const auto &meta = getSpecMetadata();
    
    if (MSLevel == SpectrumRawTypes::MSLevel::MS1)
        return getMSTKSpec(msr, getCurrentFile(), meta.first.scans[scanSel.index], mobRange);
    if (scanSel.MSMSFrameIndices.empty())
        return getMSTKSpec(msr, getCurrentFile(), meta.second.scans[scanSel.index], mobRange);
    
    // if we are here we need to get MS2 data from an IMS frame...
 
    SpectrumRaw ret;
    for (const auto i : scanSel.MSMSFrameIndices)
        ret.append(getMSTKSpec(msr, getCurrentFile(), meta.second.MSMSFrames[scanSel.index].subScans[i], mobRange));
    
    return ret;
}

void MSReadBackendMSTK::generateSpecMetadata(void)
{
    if (getCurrentFile().empty())
        return;

    auto MSTKReader = getMSTKReader();
    // load first spectrum to get file size and see if we have IMS data
    MSToolkit::Spectrum firstSpec;
    MSTKReader.readFile(getCurrentFile().c_str(), firstSpec, 0, true);
    const size_t lastScan = MSTKReader.getLastScan();
    const bool haveIMS = firstSpec.hasIonMobilityArray();
    
    SpectrumRawMetadata meta;
    
    ThreadExceptionHandler exHandler;
    
    // UNDONE: make num_threads configurable
    #pragma omp parallel num_threads(12)
    {
        auto rd = getMSTKReader();
        MSToolkit::Spectrum spec;
        SpectrumRawMetadata threadMeta;
        
        #pragma omp for schedule(static) nowait
        for (size_t i=1; i<=lastScan; ++i)
        {
            exHandler.run([&]
            {
                rd.readFile(getCurrentFile().c_str(), spec, i, true);
                const bool isMS1 = spec.getMsLevel() == 1;
                
                SpectrumRawMetadataMS *curMS1MD = (isMS1) ? &threadMeta.first : &threadMeta.second;
                curMS1MD->scans.push_back(spec.getScanNumber());
                curMS1MD->times.push_back(spec.getRTime() * 60); // NOTE: MSTK takes minutes
                curMS1MD->TICs.push_back(spec.getTIC());
                curMS1MD->BPCs.push_back(spec.getBPI());
                
                if (!isMS1)
                    threadMeta.second.isolationRanges.emplace_back(spec.getSelWindowLower(), spec.getSelWindowUpper());
            });
        }
        
        #pragma omp for schedule(static) ordered
        for(int i=0; i<getOMPNumThreads(); ++i)
        {
            exHandler.run([&]
            {
                #pragma omp ordered
                {
                    meta.first.append(threadMeta.first);
                    meta.second.append(threadMeta.second);
                }
            });
        }
    }

    exHandler.reThrow();
    

    if (haveIMS)
    {
        // For IMS-MS/MS data, the different spectra inside are frame are stored in separate spectra
        // --> move these spectra to MSMSFrames structs, so our main table only contains separate frames
        
        SpectrumRawMetadataMSMS metaMSMSFrames;
        for (size_t i=0; i<meta.second.scans.size(); ++i)
        {
            metaMSMSFrames.scans.push_back(meta.second.scans[i]);
            metaMSMSFrames.times.push_back(meta.second.times[i]);
            metaMSMSFrames.TICs.push_back(meta.second.TICs[i]);
            metaMSMSFrames.BPCs.push_back(meta.second.BPCs[i]);
            
            frameMSMSInfo fi;
            const auto curTime = meta.second.times[i];
            bool finished = false;
            while (!finished)
            {
                fi.isolationRanges.push_back(meta.second.isolationRanges[i]);
                fi.subScans.push_back(meta.second.scans[i]);
                ++i;
                if (i > lastScan)
                    break;
                if (!compareTol(curTime, meta.second.times[i]))
                    finished = true; // reached next frame
            }
            metaMSMSFrames.MSMSFrames.push_back(std::move(fi));
            if (finished)
                --i; // we reached one spec to far to find out we're finished
        }
        meta.second = std::move(metaMSMSFrames);
    }

    setSpecMetadata(std::move(meta));
}
