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
    std::string currentFile;
    static int backends; // UNDONE: for debugging
    
public:
    MSReadBackendMSTK(void) { ++backends; Rcpp::Rcout << "backends:" << backends << "\n"; }
    MSReadBackendMSTK(const MSReadBackendMSTK &o) = delete;
    ~MSReadBackendMSTK(void) { --backends; };
    
    int getBackends(void) const { return backends; }
    
    void open(const std::string &file) override { if (!currentFile.empty()) close(); currentFile = file; }
    void close(void) override { currentFile.clear(); }
    ThreadDataType getThreadData(void) const override;
    SpectrumRaw readSpectrum(const ThreadDataType &tdata, int index) const override;
};

RCPP_EXPOSED_CLASS(MSReadBackendMSTK)

#endif
