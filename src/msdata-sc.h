#ifndef PATROON_SC_H
#define PATROON_SC_H

#include <Rcpp.h>

#define PUGIXML_PATH "../pugixml/pugixml.hpp"
#include "StreamCraft/StreamCraft_mzml.hpp"
#undef PUGIXML_PATH

#include "msdata.h"

class MSReadBackendSC: public MSReadBackend
{
    void doOpen(const std::string &) override { }
    void doClose(void) override { }
    ThreadDataType doGetThreadData(void) const override;
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const override;
    
public:
    MSReadBackendSC(void) { }
    ~MSReadBackendSC(void) { }
    
    void generateSpecMetadata(void);
};

RCPP_EXPOSED_CLASS(MSReadBackendSC)

#endif
