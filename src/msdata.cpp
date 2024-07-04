#include <Rcpp.h>
#include <vector>

#include "msdata.h"
#include "mstoolkit.h"
#include "spectrum-raw.h"

namespace {

template<typename OutType, typename FuncType>
OutType applyMSData(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Scan> &scans, FuncType func)
{
    OutType ret;
    
    auto tdata = backend.getThreadData();
    
    for (auto i : scans)
    {
        ret.push_back(func(backend.readSpectrum(tdata, i)));
    }
    
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
    Rcpp::class_<MSReadBackendMSTK>("MSReadBackendMSTK")
        .derives<MSReadBackend>("MSReadBackend")
        .constructor()
        .method("getBackends", &MSReadBackendMSTK::getBackends)
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
    // UNDONE
    std::vector<int> indices(300);
    std::iota(indices.begin(), indices.end(), 1);
    
    // UNDONE
    struct EICPoint
    {
        SpectrumRawTypes::Time time;
        std::vector<SpectrumRawTypes::Mass> mzs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
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
    
    std::vector<SpectrumRawTypes::Time> times(EPSize);
    for (size_t i=0; i<EPSize; ++i)
        times[i] = EICPoints[i].time;
    
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

        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("time") = times,
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
