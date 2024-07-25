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

}


int MSReadBackendMSTK::backends = 0;

// [[Rcpp::interfaces(r, cpp)]]

MSReadBackend::ThreadDataType MSReadBackendMSTK::doGetThreadData(void) const
{
    return std::make_shared<MSToolkit::MSReader>(getMSTKReader());
}

SpectrumRaw MSReadBackendMSTK::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const
{
    MSToolkit::Spectrum s;
    auto *msr = reinterpret_cast<MSToolkit::MSReader *>(tdata.get());
    
    if (!msr->readFile(getCurrentFile().c_str(), s, scan) || s.getScanNumber() == 0)
        Rcpp::stop("Abort: invalid spectrum scan index: %d", scan);

    const bool hasMob = s.hasIonMobilityArray();
    SpectrumRaw ret(s.size(), hasMob);
    
    if (hasMob)
    {
        for(size_t i=0; i<ret.size(); ++i)
            ret.setPeak(i, s.at(i).mz, s.at(i).intensity, s.atIM(i));
    }
    else
    {
        for(size_t i=0; i<ret.size(); ++i)
            ret.setPeak(i, s.at(i).mz, s.at(i).intensity);
    }
    
    return ret;
}

void MSReadBackendMSTK::generateSpecMetadata(void)
{
    if (getCurrentFile().empty())
        return;

    auto rd = getMSTKReader();
    // load first spectrum to get file size and see if we have IMS data
    MSToolkit::Spectrum spec;
    rd.readFile(getCurrentFile().c_str(), spec, 0, true);
    const size_t lastScan = rd.getLastScan();
    const bool haveIMS = spec.hasIonMobilityArray();
    
    SpectrumRawMetadata meta;
    
    // UNDONE: disabled OpenMP for now
    for (size_t i=1; i<=lastScan; ++i)
    {
        if (i > 1)
            rd.readFile(getCurrentFile().c_str(), spec, i, true); // we already read the first scan
        const bool isMS1 = spec.getMsLevel() == 1;
        
        SpectrumRawMetadataMS *curMS1MD = (isMS1) ? &meta.first : &meta.second;
        curMS1MD->scans.push_back(spec.getScanNumber());
        curMS1MD->times.push_back(spec.getRTime() * 60); // NOTE: MSTK takes minutes
        curMS1MD->TICs.push_back(spec.getTIC());
        curMS1MD->BPCs.push_back(spec.getBPI());
        
        assert(curMS1MD->scans.size() == curMS1MD->times.size());
        
        if (!isMS1)
        {
            if (!haveIMS)
                meta.second.isolationRanges.emplace_back(NumRange<SpectrumRawTypes::Mass>(spec.getSelWindowLower(), spec.getSelWindowUpper()));
            else
            {
                frameMSMSInfo fi;
                const SpectrumRawTypes::Time curTime = spec.getRTime();
                bool finished = false;
                while (!finished)
                {
                    fi.isolationRanges.emplace_back(NumRange<SpectrumRawTypes::Mass>(spec.getSelWindowLower(),
                                                                                     spec.getSelWindowUpper()));
                    fi.subScans.push_back(spec.getScanNumber());
                    ++i;
                    if (i > lastScan)
                        break;
                    rd.readFile(getCurrentFile().c_str(), spec, i, true);
                    if (!compareTol(curTime, spec.getRTime()))
                        finished = true; // reached next frame
                }
                meta.second.MSMSFrames.emplace_back(std::move(fi));
                if (finished)
                    --i; // we reached one spec to far to find out we're finished
            }
        }
    }
    
#if 0
    // load first spectrum to get file size and see if we have IMS data
    auto rd = getMSTKReader();
    MSToolkit::Spectrum firstSpec;
    rd.readFile(getCurrentFile().c_str(), firstSpec, 0, true);
    const size_t lastScan = rd.getLastScan();
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
                
                assert(curMS1MD->scans.size() == curMS1MD->times.size());
                
                if (!isMS1)
                    threadMeta.second.isolationRanges.emplace_back(NumRange<SpectrumRawTypes::Mass>(spec.getSelWindowLower(), spec.getSelWindowUpper()));
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
    
#endif
    emplaceSpecMeta(std::move(meta));
}
