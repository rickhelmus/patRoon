#include <Rcpp.h>

#include "msdata-mem.h"
#include "spectrum-raw.h"

// [[Rcpp::interfaces(r, cpp)]]

SpectrumRaw MSReadBackendMem::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                             const SpectrumRawSelection &scanSel,
                                             const SpectrumRawTypes::MobilityRange &mobRange) const
{
    // UNDONE: handle mobRange
    
    if (MSLevel == SpectrumRawTypes::MSLevel::MS1)
        return spectraMS[scanSel.index];
    if (scanSel.MSMSFrameIndices.empty())
        return spectraMS2[scanSel.index];
    
    // if we are here we need to get MS2 data from an IMS frame...
    
    SpectrumRaw ret;
    for (const auto i : scanSel.MSMSFrameIndices)
        ret.append(spectraMS2[i]);
    return ret;
}

void MSReadBackendMem::setSpectra(const Rcpp::List &specsMS, const Rcpp::List &specsMS2)
{
    std::vector<SpectrumRaw> spsMS(specsMS.size()), spsMS2(specsMS2.size());
    
    for (int i=0; i<specsMS.size(); ++i)
    {
        const Rcpp::NumericMatrix &spm = Rcpp::as<Rcpp::NumericMatrix>(specsMS[i]);
        const Rcpp::NumericVector &mzs = spm(Rcpp::_, 0), ints = spm(Rcpp::_, 1);
        spsMS[i] = SpectrumRaw(Rcpp::as<std::vector<SpectrumRawTypes::Mass>>(mzs),
                               Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(ints));
    }
    spectraMS = std::move(spsMS);
    
    for (int i=0; i<specsMS2.size(); ++i)
    {
        const Rcpp::NumericMatrix &spm = Rcpp::as<Rcpp::NumericMatrix>(specsMS2[i]);
        const Rcpp::NumericVector &mzs = spm(Rcpp::_, 0), ints = spm(Rcpp::_, 1);
        spsMS2[i] = SpectrumRaw(Rcpp::as<std::vector<SpectrumRawTypes::Mass>>(mzs),
                                Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(ints));
    }
    spectraMS2 = std::move(spsMS2);
}
