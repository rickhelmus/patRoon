#include <Rcpp.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "msdata-mstk.h"
#include "spectrum-raw.h"

int MSReadBackendMSTK::backends = 0;

// [[Rcpp::interfaces(r, cpp)]]

MSReadBackend::ThreadDataType MSReadBackendMSTK::doGetThreadData(void) const
{
    auto ret = std::make_shared<MSToolkit::MSReader>();
    ret->addFilter(MSToolkit::MS1);
    ret->addFilter(MSToolkit::MS2);
    return ret;
}

SpectrumRaw MSReadBackendMSTK::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const
{
    MSToolkit::Spectrum s;
    auto *msr = reinterpret_cast<MSToolkit::MSReader *>(tdata.get());
    
    if (!msr->readFile(getCurrentFile().c_str(), s, scan) || s.getScanNumber() == 0)
        Rcpp::stop("Abort: invalid spectrum scan index: %d", scan);

    SpectrumRaw ret(s.size());
    for(int i=0; i<s.size(); ++i)
        ret.setPeak(i, s.at(i).mz, s.at(i).intensity);
        
    return ret;
}

const SpectrumRawMetadata &MSReadBackendMSTK::doGetSpectrumRawMetadata(void) const
{
    if (!metadataLoaded && !getCurrentFile().empty())
    {
        metadataLoaded = true;
        
        MSToolkit::MSReader reader;
        reader.addFilter(MSToolkit::MS1);
        reader.addFilter(MSToolkit::MS2);
        MSToolkit::Spectrum spec;
        
        SpectrumRawMetadata meta;
        
        while (reader.readFile(getCurrentFile().c_str(), spec, 0, true))
        {
            const bool isMS1 = spec.getMsLevel() == 1;
            
            SpectrumRawMetadataMS *curMS1MD = (isMS1) ? &meta.first : &meta.second;
            curMS1MD->scans.push_back(spec.getScanNumber());
            curMS1MD->times.push_back(spec.getRTime());
            curMS1MD->TICs.push_back(spec.getTIC());
            curMS1MD->BPCs.push_back(spec.getBPI());
            
            assert(curMS1MD->scans.size() == curMS1MD->times.size());
            
            if (!isMS1)
                meta.second.isolationRanges.emplace_back(NumRange<SpectrumRawTypes::Mass>(spec.getSelWindowLower(), spec.getSelWindowUpper()));
        }
        
        emplaceSpecMeta(std::move(meta));
    }
    
    return MSReadBackend::doGetSpectrumRawMetadata();
}


#if 0
using namespace MSToolkit;

// [[nope Rcpp::export]]
void MSToolkitTest(const std::string &path, int startIndex, int stopIndex)
{
    MSReader r;
    Spectrum s;
    int j;
    
    r.addFilter(MS1);
    r.addFilter(MS2);
    r.addFilter(MSX);
    r.addFilter(SRM);
    
    char nativeID[257];
    for (int i=startIndex; i<=stopIndex; ++i)
    {
        if (!r.readFile(path.c_str(), s, i) || s.getScanNumber() == 0)
        {
            Rcpp::Rcout << "Abort: invalid idx\n";
            return;
        }
            
        Rcpp::Rcout << "scan: " << s.getScanNumber() << "\n";
        
        char szNativeID[128];
        if (s.getNativeID(szNativeID, 128))
            Rcpp::Rcout << "success:  scan " << s.getScanNumber() << "nativeID: " << szNativeID << "\n";
        else
            Rcpp::Rcout << "failure:  scan " << s.getScanNumber() << "\n";
        
        Rcpp::Rcout << "size: " << s.size() << "\n";
        
        s.getNativeID(nativeID, 256);
        Rcpp::Rcout << nativeID << "\n";
        Rcpp::Rcout << "S\t" << s.getScanNumber() << "\t" << s.getScanNumber();
        for(j=0;j<s.sizeMZ();j++){
            Rcpp::Rcout << "\t" << s.getMZ(j);
        }
        Rcpp::Rcout << "\n";
        if(s.getRTime()>0) Rcpp::Rcout << "I\tRTime\t" << s.getRTime() << "\n";
        for(j=0;j<s.sizeZ();j++){
            Rcpp::Rcout << "Z\t" << s.atZ(j).z << "\t" << s.atZ(j).mh << "\n";
        }
        
        /*for(j=0;j<s.size();j++){
            Rcpp::Rcout << s.at(j).mz << "\t" << s.at(j).intensity << "\n";
        }*/
    }
}
#endif
