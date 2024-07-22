#ifndef PATROON_MEM_H
#define PATROON_MEM_H

#include <Rcpp.h>

#include <vector>

#include "msdata.h"
#include "spectrum-raw.h"

class MSReadBackendMem: public MSReadBackend
{
    std::vector<SpectrumRaw> spectra;
    
    void doOpen(const std::string &) override { }
    void doClose(void) override { spectra.clear(); }
    ThreadDataType doGetThreadData(void) const override { return nullptr; };
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const override;
    
public:
    MSReadBackendMem(void) { }
    ~MSReadBackendMem(void) { }
    
    void setSpectra(const Rcpp::List &specList);
};

RCPP_EXPOSED_CLASS(MSReadBackendMem)
    
#endif
