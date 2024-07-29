#ifndef PATROON_MEM_H
#define PATROON_MEM_H

#include <Rcpp.h>

#include <vector>

#include "msdata.h"
#include "spectrum-raw.h"

class MSReadBackendMem: public MSReadBackend
{
    std::vector<SpectrumRaw> spectraMS, spectraMS2;
    
    void doOpen(const std::string &) override { }
    void doClose(void) override { spectraMS.clear(); spectraMS2.clear(); }
    ThreadDataType doGetThreadData(void) const override { return nullptr; };
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                               const SpectrumRawSelection &scanSel) const override;
    
public:
    MSReadBackendMem(void) { }
    ~MSReadBackendMem(void) { }
    
    void setSpectra(const Rcpp::List &specsMS, const Rcpp::List &specsMS2);
};

RCPP_EXPOSED_CLASS(MSReadBackendMem)
    
#endif
