#ifndef PATROON_MSTK_H
#define PATROON_MSTK_H

#include <Rcpp.h>

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
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                               const SpectrumRawSelection &scanSel,
                               const SpectrumRawTypes::MobilityRange &mobRange) const override;
    
public:
    MSReadBackendMSTK(void) { ++backends; Rcpp::Rcout << "constr: backends:" << backends << "\n"; }
    ~MSReadBackendMSTK(void) { --backends; Rcpp::Rcout << "destr: backends:" << backends << "\n"; }

    void generateSpecMetadata(void);
    int getBackends(void) const { return backends; }
};

RCPP_EXPOSED_CLASS(MSReadBackendMSTK)

#endif
