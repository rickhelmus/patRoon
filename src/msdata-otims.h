#ifndef PATROON_OTIMS_H
#define PATROON_OTIMS_H

#ifdef WITH_OTIMS

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
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                               const SpectrumRawSelection &scanSel,
                               const SpectrumRawTypes::MobilityRange &mobRange) const override;
};

#else

#include "msdata.h"

class MSReadBackendOTIMS: public MSReadBackendUnavailable
{
public:
    MSReadBackendOTIMS(void) : MSReadBackendUnavailable("opentims") { }
    ~MSReadBackendOTIMS(void) { }
};

#endif // WITH_OTIMS

RCPP_EXPOSED_CLASS(MSReadBackendOTIMS)

#endif
