#ifndef PATROON_MSTK_H
#define PATROON_MSTK_H

#include <Rcpp.h>

#include <memory>
#include <string>

#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

#include "msdata.h"

class MSReadBackendMSTK: public MSReadBackend
{
    static int backends; // UNDONE: for debugging
    mutable bool metadataLoaded = false; // mutable: to allow lazy loading

    void doOpen(const std::string &) override { }
    void doClose(void) override { metadataLoaded = false; }
    ThreadDataType doGetThreadData(void) const override;
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const override;
    
protected:
    const SpectrumRawMetadata &doGetSpectrumRawMetadata(void) const override;
    
public:
    MSReadBackendMSTK(void) { ++backends; Rcpp::Rcout << "constr: backends:" << backends << "\n"; }
    MSReadBackendMSTK(const MSReadBackendMSTK &) = delete;
    ~MSReadBackendMSTK(void) { --backends; Rcpp::Rcout << "destr: backends:" << backends << "\n"; }
    
    int getBackends(void) const { return backends; }
};

RCPP_EXPOSED_CLASS(MSReadBackendMSTK)

#endif
