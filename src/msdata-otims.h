#ifndef PATROON_OTIMS_H
#define PATROON_OTIMS_H

#include <Rcpp.h>

#include "msdata.h"
#include "opentims++/opentims.h"

#include <memory>

class MSReadBackendOTIMS: public MSReadBackend
{
    std::unique_ptr<TimsDataHandle> handle;
    
    void doOpen(const std::string &file) override;
    void doClose(void) override;
    ThreadDataType doGetThreadData(void) const override;
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const override;
};

RCPP_EXPOSED_CLASS(MSReadBackendOTIMS)

#endif
