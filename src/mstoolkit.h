#ifndef PATROON_MSTOOLKIT_H
#define PATROON_MSTOOLKIT_H

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

    void doOpen(const std::string &) override { }
    void doClose(void) override { }
    ThreadDataType doGetThreadData(void) const override;
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, int index) const override;
    
public:
    MSReadBackendMSTK(void) { ++backends; Rcpp::Rcout << "backends:" << backends << "\n"; }
    MSReadBackendMSTK(const MSReadBackendMSTK &) = delete;
    ~MSReadBackendMSTK(void) { --backends; };
    
    int getBackends(void) const { return backends; }
};

RCPP_EXPOSED_CLASS(MSReadBackendMSTK)

#endif
