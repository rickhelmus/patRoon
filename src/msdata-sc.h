#ifndef PATROON_SC_H
#define PATROON_SC_H

#include <Rcpp.h>

#include "msdata.h"

class MSReadBackendSC: public MSReadBackend
{
    void doOpen(const std::string &) override { }
    void doClose(void) override { }
    ThreadDataType doGetThreadData(void) const override;
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                               const SpectrumRawSelection &scanSel,
                               const SpectrumRawTypes::MobilityRange &mobRange,
                               SpectrumRawTypes::Intensity minIntensityIMS) const override;
    
public:
    MSReadBackendSC(void) { }
    ~MSReadBackendSC(void) { }
    
    void generateSpecMetadata(void);
};

RCPP_EXPOSED_CLASS(MSReadBackendSC)

#endif
