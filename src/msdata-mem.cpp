#include <Rcpp.h>

#include "msdata-mem.h"
#include "spectrum-raw.h"

// [[Rcpp::interfaces(r, cpp)]]

SpectrumRaw MSReadBackendMem::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const
{
    const auto meta = getSpecMetadata();
    
    // first check if it's an MS spectrum
    auto it = std::lower_bound(getSpecMetadata().first.scans.begin(), getSpecMetadata().first.scans.end(), scan);
    if (it != getSpecMetadata().first.scans.end())
        return spectra[std::distance(getSpecMetadata().first.scans.begin(), it)];

    // MS/MS spectrum then?
    it = std::lower_bound(getSpecMetadata().second.scans.begin(), getSpecMetadata().second.scans.end(), scan);
    if (it != getSpecMetadata().second.scans.end())
        return spectra[std::distance(getSpecMetadata().second.scans.begin(), it)];
    
    Rcpp::stop("Abort: invalid spectrum scan index: %d", scan);
    
    return SpectrumRaw(); // unreached
}

void MSReadBackendMem::setSpectra(const Rcpp::List &specList)
{
    std::vector<SpectrumRaw> sps(specList.size());
    for (int i=0; i<specList.size(); ++i)
    {
        const Rcpp::NumericMatrix spm = Rcpp::as<Rcpp::NumericMatrix>(specList[i]);
        const Rcpp::NumericVector mzs = spm(Rcpp::_, 0), ints = spm(Rcpp::_, 1);
        sps[i] = SpectrumRaw(Rcpp::as<std::vector<SpectrumRawTypes::Mass>>(mzs),
                             Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(ints));
    }
    spectra = std::move(sps);
}
