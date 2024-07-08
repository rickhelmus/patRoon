#ifndef PATROON_SC_H
#define PATROON_SC_H

#include <Rcpp.h>

#define PUGIXML_PATH "../pugixml/pugixml.hpp"
#include "StreamCraft/StreamCraft_mzml.hpp"
#undef PUGIXML_PATH

#include "msdata.h"

class MSReadBackendSC: public MSReadBackend
{
    mutable bool metadataLoaded = false; // mutable: to allow lazy loading

    void doOpen(const std::string &) override { }
    void doClose(void) override { metadataLoaded = false; }
    ThreadDataType doGetThreadData(void) const override;
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const override;
    
protected:
    const SpectrumRawMetadata &doGetSpectrumRawMetadata(void) const override;
    
public:
    MSReadBackendSC(void) { }
    ~MSReadBackendSC(void) { }
};

RCPP_EXPOSED_CLASS(MSReadBackendSC)

#endif
