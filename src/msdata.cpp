#include <Rcpp.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "msdata.h"
#include "msdata-mem.h"
#include "msdata-mstk.h"
#include "msdata-sc.h"
#include "spectrum-raw.h"

namespace {

template<typename OutType, typename FuncType>
OutType applyMSData(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Scan> &scans, FuncType func)
{
    OutType ret(scans.size());
    ThreadExceptionHandler exHandler;

    // UNDONE: make num_threads configurable
    #pragma omp parallel num_threads(12)
    {
        auto tdata = backend.getThreadData();
        #pragma omp for
        for (size_t i=0; i<scans.size(); ++i)
        {
            exHandler.run([&]{ ret[i] = func(backend.readSpectrum(tdata, scans[i])); });
        }
    }
    
    exHandler.reThrow();
    
    return ret;
}
    
}


void MSReadBackend::open(const std::string &file)
{
    close();
    doOpen(file);
    currentFile = file;
}

void MSReadBackend::close(void)
{
    if (!currentFile.empty())
    {
        doClose();
        currentFile.clear();
        specMetadata.first.clear();
        specMetadata.second.clear();
    }
}

RCPP_MODULE(MSReadBackend)
{
    Rcpp::class_<MSReadBackend>("MSReadBackend")
        .method("open", &MSReadBackend::open)
        .method("close", &MSReadBackend::close)
    ;
    Rcpp::class_<MSReadBackendMem>("MSReadBackendMem")
        .derives<MSReadBackend>("MSReadBackend")
        .constructor()
        .method("setSpecMetadata", &MSReadBackendMem::setSpecMetadata)
        .method("setSpectra", &MSReadBackendMem::setSpectra)
    ;
    Rcpp::class_<MSReadBackendMSTK>("MSReadBackendMSTK")
        .derives<MSReadBackend>("MSReadBackend")
        .constructor()
        .method("getBackends", &MSReadBackendMSTK::getBackends)
    ;
    Rcpp::class_<MSReadBackendSC>("MSReadBackendSC")
        .derives<MSReadBackend>("MSReadBackend")
        .constructor()
    ;
}

// [[Rcpp::export]]
Rcpp::DataFrame getMSSpectrum(const MSReadBackend &backend, int index)
{
    const std::vector<SpectrumRawTypes::Scan> indices(1, index);
    
    const auto sfunc = [](const SpectrumRaw &spec)
    {
        return spec;
    };
    
    const std::vector<SpectrumRaw> spectra = applyMSData<std::vector<SpectrumRaw>>(backend, indices, sfunc);
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spectra[0].getMZs(),
                                   Rcpp::Named("intensity") = spectra[0].getIntensities());
}

// [[Rcpp::export]]
std::vector<SpectrumRawTypes::Scan> getScans(const MSReadBackend &backend, SpectrumRawTypes::Mass timeStart,
                                             SpectrumRawTypes::Mass timeEnd, int MSLevel,
                                             SpectrumRawTypes::Mass isoStart, SpectrumRawTypes::Mass isoEnd)
{
    return getSpecScanIndices(backend.getSpecMetadata(), NumRange<SpectrumRawTypes::Time>(timeStart, timeEnd),
                              (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2,
                              NumRange<SpectrumRawTypes::Mass>(isoStart, isoEnd));
}

// [[Rcpp::export]]
Rcpp::List getEICList(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Mass> &startMZs,
                      const std::vector<SpectrumRawTypes::Mass> &endMZs)
{
    // UNDONE: keep this local here?
    struct EICPoint
    {
        std::vector<SpectrumRawTypes::Mass> mzs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        EICPoint(void) = default;
        EICPoint(size_t s) : mzs(s), intensities(s) { }
    };
    
    const auto EICCount = startMZs.size();
    
    const auto sfunc = [&](const SpectrumRaw &spec)
    {
        EICPoint ret(EICCount);
        
        for (size_t i=0; i<EICCount; ++i)
        {
            ret.mzs[i] = ret.intensities[i] = 0.0;

            // NOTE: assume spec is mz sorted
            const auto mzIt = std::lower_bound(spec.getMZs().begin(), spec.getMZs().end(), startMZs[i]);
            if (mzIt == spec.getMZs().end())
                continue;
            
            for (size_t j=std::distance(spec.getMZs().begin(), mzIt); j<spec.size(); ++j)
            {
                const auto mz = spec.getMZs()[j];
                
                // NOTE: assume spec is mz sorted
                if (mz > endMZs[i])
                    break; 
                
                const auto inten = spec.getIntensities()[j];
                ret.intensities[i] += inten;
                ret.mzs[i] += mz * inten;
            }
            
            if (ret.intensities[i] > 0)
                ret.mzs[i] /= ret.intensities[i]; // weighted mean
        }
        
        return ret;
    };

    const auto &specMetaMS = backend.getSpecMetadata().first;
    const auto EICPoints = applyMSData<std::vector<EICPoint>>(backend, specMetaMS.scans, sfunc);
    const auto EPSize = EICPoints.size();
    
    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        std::vector<SpectrumRawTypes::Mass> mzs(EPSize);
        std::vector<SpectrumRawTypes::Intensity> intensities(EPSize);
        for (size_t j=0; j<EPSize; ++j)
        {
            mzs[j] = EICPoints[j].mzs[i];
            intensities[j] = EICPoints[j].intensities[i];
        }

        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("time") = specMetaMS.times,
                                         Rcpp::Named("intensity") = intensities,
                                         Rcpp::Named("mz") = mzs);
    }
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::DataFrame getMSMetadata(MSReadBackend &backend, int msLevel)
{
    const auto &meta = backend.getSpecMetadata();
    
    if (msLevel == 1)
        return(Rcpp::DataFrame::create(Rcpp::Named("scan") = meta.first.scans,
                                       Rcpp::Named("time") = meta.first.times,
                                       Rcpp::Named("TIC") = meta.first.TICs,
                                       Rcpp::Named("BPC") = meta.first.BPCs));
    
    // msLevel == 2
    
    const auto size = meta.second.isolationRanges.size();
    std::vector<SpectrumRawTypes::Mass> isolationMins(size), isolationMaxs(size);
    for (size_t i=0; i<size; ++i)
    {
        isolationMins[i] = meta.second.isolationRanges[i].start;
        isolationMaxs[i] = meta.second.isolationRanges[i].end;
    }
    
    return Rcpp::DataFrame::create(Rcpp::Named("time") = meta.second.times,
                                   Rcpp::Named("TIC") = meta.second.TICs,
                                   Rcpp::Named("BPC") = meta.second.BPCs,
                                   Rcpp::Named("isolationRangeMin") = isolationMins,
                                   Rcpp::Named("isolationRangeMax") = isolationMaxs);
}
