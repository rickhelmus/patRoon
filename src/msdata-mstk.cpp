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

int getMSTKLastScan(const std::string &file)
{
    // load first spectrum to get file size
    auto rd = getMSTKReader();
    MSToolkit::Spectrum firstSpec;
    rd.readFile(file.c_str(), firstSpec, 0, true);
    return rd.getLastScan();
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

const SpectrumRawMetadata &MSReadBackendMSTK::doGetSpectrumRawMetadata(void) const
{
    if (!metadataLoaded && !getCurrentFile().empty())
    {
        metadataLoaded = true;
        
        const size_t lastScan = getMSTKLastScan(getCurrentFile());
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
        
        emplaceSpecMeta(std::move(meta));
    }
    
    return MSReadBackend::doGetSpectrumRawMetadata();
}
