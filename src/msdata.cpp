#include <Rcpp.h>
#include <vector>

#include "mstoolkit.h"
#include "spectrum-raw.h"

namespace {

template<typename OutType, typename BackendType, typename FuncType>
OutType applyMSData(const BackendType &backend, const std::vector<int> &indices, FuncType func)
{
    OutType ret;
    
    auto tdata = backend.getThreadData();
    
    for (auto i : indices)
    {
        ret.push_back(func(backend.readSpectrum(tdata, i)));
    }
    
    return ret;
}
    
}


// [[Rcpp::export]]
Rcpp::DataFrame getMSSpectrum(const MSToolkitBackend &backend, int index)
{
    const std::vector<int> indices(1, index);
    
    const auto sfunc = [](const SpectrumRaw &spec)
    {
        return spec;
    };
    
    const std::vector<SpectrumRaw> spectra = applyMSData<std::vector<SpectrumRaw>>(backend, indices, sfunc);
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spectra[0].getMZs(),
                                   Rcpp::Named("intensity") = spectra[0].getIntensities());
}
