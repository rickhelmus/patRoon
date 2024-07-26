#include <Rcpp.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "msdata.h"
#include "msdata-mem.h"
#include "msdata-mstk.h"
#include "msdata-otims.h"
#include "msdata-sc.h"
#include "spectrum-raw.h"

namespace {

template<typename OutType, typename FuncType, typename... Args>
OutType applyMSData(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Scan> &scans, FuncType func,
                    Args... args)
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
            exHandler.run([&]{ ret[i] = func(backend.readSpectrum(tdata, scans[i]), args...); });
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
        .method("setSpectra", &MSReadBackendMem::setSpectra)
    ;
    Rcpp::class_<MSReadBackendMSTK>("MSReadBackendMSTK")
        .derives<MSReadBackend>("MSReadBackend")
        .constructor()
        .method("generateSpecMetadata", &MSReadBackendMSTK::generateSpecMetadata)
        .method("getBackends", &MSReadBackendMSTK::getBackends)
    ;
    Rcpp::class_<MSReadBackendOTIMS>("MSReadBackendOTIMS")
        .derives<MSReadBackend>("MSReadBackend")
        .constructor()
    ;
    Rcpp::class_<MSReadBackendSC>("MSReadBackendSC")
        .derives<MSReadBackend>("MSReadBackend")
        .constructor()
        .method("generateSpecMetadata", &MSReadBackendSC::generateSpecMetadata)
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
    
    if (!spectra[0].getMobilities().empty())
        return Rcpp::DataFrame::create(Rcpp::Named("mz") = spectra[0].getMZs(),
                                       Rcpp::Named("intensity") = spectra[0].getIntensities(),
                                       Rcpp::Named("mobility") = spectra[0].getMobilities());
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spectra[0].getMZs(),
                                   Rcpp::Named("intensity") = spectra[0].getIntensities());
}

// [[Rcpp::export]]
std::vector<SpectrumRawTypes::Scan> getScans(const MSReadBackend &backend, SpectrumRawTypes::Mass timeStart,
                                             SpectrumRawTypes::Mass timeEnd, int MSLevel,
                                             SpectrumRawTypes::Mass isoStart, SpectrumRawTypes::Mass isoEnd)
{
    return getSpecScanIndices(backend.getSpecMetadata(), makeNumRange(timeStart, timeEnd),
                              (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2,
                              makeNumRange(isoStart, isoEnd));
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
Rcpp::DataFrame getMSMetadata(const MSReadBackend &backend, int msLevel)
{
    const auto &meta = backend.getSpecMetadata();
    
    if (msLevel == 1)
        return(Rcpp::DataFrame::create(Rcpp::Named("scan") = meta.first.scans,
                                       Rcpp::Named("time") = meta.first.times,
                                       Rcpp::Named("TIC") = meta.first.TICs,
                                       Rcpp::Named("BPC") = meta.first.BPCs));
    
    // msLevel == 2
    
    if (!meta.second.isolationRanges.empty()) // non IMS
    {
        const auto size = meta.second.isolationRanges.size();
        std::vector<SpectrumRawTypes::Mass> isolationMins(size), isolationMaxs(size);
        for (size_t i=0; i<size; ++i)
        {
            isolationMins[i] = meta.second.isolationRanges[i].start;
            isolationMaxs[i] = meta.second.isolationRanges[i].end;
        }
        
        return Rcpp::DataFrame::create(Rcpp::Named("scan") = meta.second.scans,
                                       Rcpp::Named("time") = meta.second.times,
                                       Rcpp::Named("TIC") = meta.second.TICs,
                                       Rcpp::Named("BPC") = meta.second.BPCs,
                                       Rcpp::Named("isolationRangeMin") = isolationMins,
                                       Rcpp::Named("isolationRangeMax") = isolationMaxs);
    }
    
    // MS2 / IMS --> convert to tidy format
    
    std::vector<SpectrumRawTypes::Scan> scans;
    std::vector<SpectrumRawTypes::Time> times;
    std::vector<SpectrumRawTypes::Intensity> TICs, BPCs;
    std::vector<SpectrumRawTypes::Mass> isolationMins, isolationMaxs;
    std::vector<SpectrumRawTypes::Scan> subScans, subScanEnds;
    
    for (size_t i=0; i<meta.second.scans.size(); ++i)
    {
        const frameMSMSInfo &fi = meta.second.MSMSFrames[i];
        for (size_t j=0; j<fi.isolationRanges.size(); ++j)
        {
            scans.push_back(meta.second.scans[i]);
            times.push_back(meta.second.times[i]);
            TICs.push_back(meta.second.TICs[i]);
            BPCs.push_back(meta.second.BPCs[i]);
            isolationMins.push_back(fi.isolationRanges[j].start);
            isolationMaxs.push_back(fi.isolationRanges[j].end);
            subScans.push_back(fi.subScans[j]);
            if (!fi.subScanEnds.empty())
                subScanEnds.push_back(fi.subScanEnds[j]);
        }
    }
    
    auto ret = Rcpp::DataFrame::create(Rcpp::Named("scan") = scans,
                                       Rcpp::Named("time") = times,
                                       Rcpp::Named("TIC") = TICs,
                                       Rcpp::Named("BPC") = BPCs,
                                       Rcpp::Named("isolationRangeMin") = isolationMins,
                                       Rcpp::Named("isolationRangeMax") = isolationMaxs,
                                       Rcpp::Named("subScan") = subScans);
    if (!subScanEnds.empty())
        ret["subScanEnd"] = subScanEnds;
    
    return ret;
}

// [[Rcpp::export]]
void setSpecMetadata(MSReadBackend &backend, const Rcpp::DataFrame &mdMS, const Rcpp::DataFrame &mdMSMS)
{
    SpectrumRawMetadata meta;

    // MS
    meta.first.scans = Rcpp::as<std::vector<SpectrumRawTypes::Scan>>(mdMS["scan"]);
    meta.first.times = Rcpp::as<std::vector<SpectrumRawTypes::Time>>(mdMS["time"]);
    meta.first.TICs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMS["TIC"]);
    meta.first.BPCs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMS["BPC"]);
    
    // MSMS
    std::vector<SpectrumRawTypes::Scan> R_scans = mdMSMS["scan"];
    std::vector<SpectrumRawTypes::Time> R_times = mdMSMS["time"];
    std::vector<SpectrumRawTypes::Intensity> R_TICs = mdMSMS["TIC"], R_BPCs = mdMSMS["BPC"];
    std::vector<SpectrumRawTypes::Mass> R_isoStarts = mdMSMS["isolationRangeMin"], R_isoEnds = mdMSMS["isolationRangeMax"];
    
    const std::vector<std::string> cn = mdMSMS.names();
    
    if (std::find(cn.begin(), cn.end(), "isolationStart") != cn.end()) // non-IMS
    {
        meta.second.scans = std::move(R_scans);
        meta.second.times = std::move(R_times);
        meta.second.TICs = std::move(R_TICs);
        meta.second.BPCs = std::move(R_BPCs);
        
        for (size_t i=0; i<R_isoStarts.size(); ++i)
            meta.second.isolationRanges.emplace_back(makeNumRange(R_isoStarts[i], R_isoEnds[i]));
    }
    else
    {
        std::vector<SpectrumRawTypes::Scan> R_subScans = mdMSMS["subScan"], R_subScanEnds;
        if (std::find(cn.begin(), cn.end(), "subScanEnd") != cn.end())
            R_subScanEnds = Rcpp::as<std::vector<SpectrumRawTypes::Scan>>(mdMSMS["subScanEnd"]);
        
        SpectrumRawTypes::Scan curScan;
        frameMSMSInfo curFI;
        for (size_t i=0; i<R_scans.size(); ++i)
        {
            const SpectrumRawTypes::Scan sc = R_scans[i];
            if (i == 0 || curScan != sc)
            {
                curScan = sc;
                meta.second.scans.push_back(sc);
                meta.second.times.push_back(R_times[i]);
                meta.second.TICs.push_back(R_TICs[i]);
                meta.second.BPCs.push_back(R_BPCs[i]);
                if (!curFI.empty())
                {
                    meta.second.MSMSFrames.emplace_back(std::move(curFI));
                    curFI.clear();
                }
            }
            curFI.isolationRanges.emplace_back(makeNumRange(R_isoStarts[i], R_isoEnds[i]));
            curFI.subScans.push_back(R_subScans[i]);
            if (!R_subScanEnds.empty())
                curFI.subScanEnds.push_back(R_subScanEnds[i]);
        }
        if (!curFI.empty())
            meta.second.MSMSFrames.emplace_back(std::move(curFI));
    }
    
    backend.emplaceSpecMeta(std::move(meta));
}

// [[Rcpp::export]]
Rcpp::List getMSPeakLists(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Time> &startTimes,
                          const std::vector<SpectrumRawTypes::Time> &endTimes,
                          const std::vector<SpectrumRawTypes::Mass> &precursorMZs, int MSLevel,
                          SpectrumRawTypes::Mass isoWindow, const std::string &method,
                          SpectrumRawTypes::Mass mzWindow, unsigned minAbundance, unsigned topMost,
                          SpectrumRawTypes::Intensity minIntensityPre, SpectrumRawTypes::Intensity minIntensityPost)
{
    const auto entries = startTimes.size();
    const auto clMethod = clustMethodFromStr(method);
    const auto specMeta = backend.getSpecMetadata();
    const auto specFilter = SpectrumRawFilter().setMinIntensity(minIntensityPre).setTopMost(topMost);
    
    std::vector<SpectrumRaw> averagedSpectra(entries);
    ThreadExceptionHandler exHandler;

    // NOTE: we could use applyMSData() here, but parallelizing the outer loop below is much faster
    #pragma omp parallel num_threads(12)
    {
        auto tdata = backend.getThreadData();
        
        #pragma omp for
        for (size_t i=0; i<entries; ++i)
        {
            exHandler.run([&]{
                const auto scans = getSpecScanIndices(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                                      (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2,
                                                      makeNumRange(precursorMZs[i] - isoWindow, precursorMZs[i] + isoWindow));
                std::vector<SpectrumRaw> spectra(scans.size());
                for (size_t j=0; j<scans.size(); ++j)
                    spectra[j] = filterSpectrumRaw(backend.readSpectrum(tdata, scans[j]), specFilter, precursorMZs[i]);
                averagedSpectra[i] = averageSpectraRaw(spectra, clMethod, mzWindow, true, minIntensityPost, minAbundance);
            });
        }
    }
    
    exHandler.reThrow();
    
    Rcpp::List ret(entries);
    for (size_t i=0; i<entries; ++i)
    {
        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = averagedSpectra[i].getMZs(),
                                         Rcpp::Named("intensity") = averagedSpectra[i].getIntensities());
    }
        
    return ret;
}
