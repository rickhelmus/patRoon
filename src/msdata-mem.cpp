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

void MSReadBackendMem::setSpecMetadata(const Rcpp::DataFrame &mdMS, const Rcpp::DataFrame &mdMSMS)
{
    SpectrumRawMetadata meta;
    
    // MS
    meta.first.scans = Rcpp::as<std::vector<SpectrumRawTypes::Scan>>(mdMS["scan"]);
    meta.first.times = Rcpp::as<std::vector<SpectrumRawTypes::Time>>(mdMS["time"]);
    meta.first.TICs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMS["TIC"]);
    meta.first.BPCs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMS["BPC"]);
    
    // MSMS
    meta.second.scans = Rcpp::as<std::vector<SpectrumRawTypes::Scan>>(mdMSMS["scan"]);
    meta.second.times = Rcpp::as<std::vector<SpectrumRawTypes::Time>>(mdMSMS["time"]);
    meta.second.TICs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMSMS["TIC"]);
    meta.second.BPCs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMSMS["BPC"]);
    const std::vector<SpectrumRawTypes::Mass> isoStarts = mdMSMS["isolationStart"];
    const std::vector<SpectrumRawTypes::Mass> isoEnds = mdMSMS["isolationEnd"];
    for (size_t i=0; i<isoStarts.size(); ++i)
        meta.second.isolationRanges.push_back(NumRange<SpectrumRawTypes::Mass>(isoStarts[i], isoEnds[i]));
    
    emplaceSpecMeta(std::move(meta));
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
