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

// [[Rcpp::export]]
Rcpp::List getEICList(const MSReadBackend &backend, const std::vector<double> &startMZs,
                      const std::vector<double> &endMZs)
{
    // UNDONE
    std::vector<int> indices(300);
    std::iota(indices.begin(), indices.end(), 1);
    
    // UNDONE
    struct EIC { std::vector<double> time, mz, intensity; };
    struct EICPoint
    {
        double time;
        std::vector<double> mzs, intensities;
        EICPoint(size_t s) : mzs(s), intensities(s) { }
    };
    
    const auto EICCount = startMZs.size();
    
    const auto sfunc = [&](const SpectrumRaw &spec)
    {
        EICPoint ret(EICCount);
        
        ret.time = spec.getTime();

        for (size_t i=0; i<EICCount; ++i)
        {
            ret.mzs[i] = ret.intensities[i] = 0.0;

            for (size_t j = 0; j<spec.size(); ++j)
            {
                const auto mz = spec.getMZs()[j];
                
                if (mz < startMZs[i])
                    continue;
                if (mz > endMZs[i])
                    break; // assume spec is mz sorted
                
                const auto inten = spec.getIntensities()[j];
                ret.intensities[i] += inten;
                ret.mzs[i] += mz * inten;
            }
            
            if (ret.intensities[i] > 0)
                ret.mzs[i] /= ret.intensities[i]; // weighted mean
        }
        
        return ret;
    };
    
    const auto EICPoints = applyMSData<std::vector<EICPoint>>(backend, indices, sfunc);
    const auto EPSize = EICPoints.size();
    
    std::vector<double> times(EPSize);
    for (size_t i=0; i<EPSize; ++i)
        times[i] = EICPoints[i].time;
    
    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        std::vector<double> mzs(EPSize), intensities(EPSize);
        for (size_t j=0; j<EPSize; ++j)
        {
            mzs[j] = EICPoints[j].mzs[i];
            intensities[j] = EICPoints[j].intensities[i];
        }

        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("time") = times,
                                         Rcpp::Named("intensity") = intensities,
                                         Rcpp::Named("mz") = mzs);
    }
    
    
    return ret;
}

