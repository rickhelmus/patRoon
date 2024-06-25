#ifndef PATROON_MSTOOLKIT_H
#define PATROON_MSTOOLKIT_H

#include <Rcpp.h>

#include <memory>
#include <string>

#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

class SpectrumRaw;

class MSToolkitBackend
{
public:
    using ThreadDataType = std::unique_ptr<MSToolkit::MSReader>;
    
private:
    std::string currentFile;
    static int backends; // UNDONE: for debugging
    
public:
    MSToolkitBackend(void) { ++backends; Rcpp::Rcout << "backends:" << backends << "\n"; }
    MSToolkitBackend(const MSToolkitBackend &o) = delete;
    ~MSToolkitBackend(void) { --backends; };
    
    int getBackends(void) const { return backends; }
    
    void open(const std::string &file) { if (!currentFile.empty()) close(); currentFile = file; }
    void close(void) { currentFile.clear(); }
    
    ThreadDataType getThreadData(void) const;
    
    SpectrumRaw readSpectrum(ThreadDataType &tdata, int index) const;
};

RCPP_EXPOSED_CLASS(MSToolkitBackend)

#endif
