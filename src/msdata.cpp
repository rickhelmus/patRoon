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
OutType applyMSData(const MSReadBackend &backend, SpectrumRawTypes::MSLevel MSLevel,
                    const std::vector<SpectrumRawSelection> &scanSels,
                    const SpectrumRawTypes::MobilityRange &mobRange, FuncType func, Args... args)
{
    OutType ret(scanSels.size());
    ThreadExceptionHandler exHandler;

    // UNDONE: make num_threads configurable
    #pragma omp parallel num_threads(12)
    {
        auto tdata = backend.getThreadData();
        #pragma omp for
        for (size_t i=0; i<scanSels.size(); ++i)
        {
            exHandler.run([&]{ ret[i] = func(backend.readSpectrum(tdata, MSLevel, scanSels[i], mobRange), scanSels[i],
                                             args...); });
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
    SpectrumRawSelection sel; sel.index = index;
    const std::vector<SpectrumRawSelection> sels(1, sel);
    
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &)
    {
        return spec;
    };
    
    const auto spectra = applyMSData<std::vector<SpectrumRaw>>(backend, SpectrumRawTypes::MSLevel::MS1, sels,
                                                               SpectrumRawTypes::MobilityRange(), sfunc);
    
    if (!spectra[0].getMobilities().empty())
        return Rcpp::DataFrame::create(Rcpp::Named("mz") = spectra[0].getMZs(),
                                       Rcpp::Named("intensity") = spectra[0].getIntensities(),
                                       Rcpp::Named("mobility") = spectra[0].getMobilities());
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spectra[0].getMZs(),
                                   Rcpp::Named("intensity") = spectra[0].getIntensities());
}

// [[Rcpp::export]]
Rcpp::DataFrame getScans(const MSReadBackend &backend, SpectrumRawTypes::Mass timeStart, SpectrumRawTypes::Mass timeEnd,
                         int MSLevel, SpectrumRawTypes::Mass isoStart, SpectrumRawTypes::Mass isoEnd)
{
    const auto sels = getSpecRawSelections(backend.getSpecMetadata(), makeNumRange(timeStart, timeEnd),
                                           (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2,
                                           makeNumRange(isoStart, isoEnd));
    std::vector<SpectrumRawTypes::Scan> inds, frInds;
    for (const auto &sel : sels)
    {
        if (sel.MSMSFrameIndices.empty())
            inds.push_back(sel.index);
        else
        {
            for (const auto fri : sel.MSMSFrameIndices)
            {
                inds.push_back(sel.index);
                frInds.push_back(fri);
            }
        }
    }
    
    auto ret = Rcpp::DataFrame::create(Rcpp::Named("index") = inds);
    if (!frInds.empty())
        ret["MSMSFrameInd"] = frInds;
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::List getEICList(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Mass> &startMZs,
                      const std::vector<SpectrumRawTypes::Mass> &endMZs,
                      const std::vector<SpectrumRawTypes::Time> &startTimes,
                      const std::vector<SpectrumRawTypes::Time> &endTimes,
                      const std::vector<SpectrumRawTypes::Mobility> &startMobs,
                      const std::vector<SpectrumRawTypes::Mobility> &endMobs)
{
    // UNDONE: keep this local here?
    struct EICPoint
    {
        std::vector<SpectrumRawTypes::Mass> mzs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        EICPoint(void) = default;
        EICPoint(size_t s, bool hasMob) : mzs(s), intensities(s), mobilities((hasMob) ? s : 0) { }
    };
    
    const auto EICCount = startMZs.size();
    const auto &specMetaMS = backend.getSpecMetadata().first;
    
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel)
    {
        const auto time = specMetaMS.times[ssel.index];
        const bool hasMob = spec.hasMobilities();
        EICPoint ret(EICCount, hasMob);
        
        for (size_t i=0; i<EICCount; ++i)
        {
            if (time < startTimes[i] || time > endTimes[i])
                continue;
            
            size_t startInd = 0;
            
            // NOTE: for non-IMS data assume spec is mz sorted
            
            if (!hasMob)
            {
                // use lower bound to speedup search for first mass
                const auto mzIt = std::lower_bound(spec.getMZs().begin(), spec.getMZs().end(), startMZs[i]);
                if (mzIt == spec.getMZs().end())
                    continue;
                startInd = std::distance(spec.getMZs().begin(), mzIt);
            }
            
            for (size_t j=startInd; j<spec.size(); ++j)
            {
                const auto mz = spec.getMZs()[j];
                const auto mob = (hasMob) ? spec.getMobilities()[j] : 0;
                
                if (hasMob)
                {
                    if (mz < startMZs[i] || mz > endMZs[i])
                        continue;
                    if (mob < startMobs[i] || mob > endMobs[i])
                        continue; // UNDONE: optimize eg for TIMS data where mobilities are sorted?
                }
                else if (mz > endMZs[i])
                    break; 
                
                const auto inten = spec.getIntensities()[j];
                ret.intensities[i] += inten;
                ret.mzs[i] += mz * inten;
                if (hasMob)
                    ret.mobilities[i] += mob * inten;
            }

            if (ret.intensities[i] > 0)
            {
                // weighted mean
                ret.mzs[i] /= ret.intensities[i];
                if (hasMob)
                    ret.mobilities[i] /= ret.intensities[i];
            }
        }
        
        return ret;
    };
    
    std::vector<SpectrumRawTypes::TimeRange> timeRanges(EICCount);
    for (size_t i=0; i<EICCount; ++i)
        timeRanges.emplace_back(startTimes[i], endTimes[i]);
    const auto sels = getSpecRawSelections(backend.getSpecMetadata(), timeRanges, SpectrumRawTypes::MSLevel::MS1,
                                           SpectrumRawTypes::IsolationRange());
    const auto EICPoints = applyMSData<std::vector<EICPoint>>(backend, SpectrumRawTypes::MSLevel::MS1, sels,
                                                              SpectrumRawTypes::MobilityRange(), sfunc);
    const auto EPSize = EICPoints.size();
    
    if (EPSize == 0)
        return Rcpp::List();
    
    // NOTE: assume all MS spectra have or have not IMS data
    const bool hasMob = !EICPoints[0].mobilities.empty();
    
    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        std::vector<SpectrumRawTypes::Time> times;
        std::vector<SpectrumRawTypes::Mass> mzs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        for (size_t j=0; j<EPSize; ++j)
        {
            const auto t = specMetaMS.times[sels[j].index];
            if (t < startTimes[i] || t > endTimes[i])
                continue;
            times.push_back(t);
            mzs.push_back(EICPoints[j].mzs[i]);
            intensities.push_back(EICPoints[j].intensities[i]);
            if (hasMob)
                mobilities.push_back(EICPoints[j].mobilities[i]);
        }
        
        auto df = Rcpp::DataFrame::create(Rcpp::Named("time") = times,
                                          Rcpp::Named("intensity") = intensities,
                                          Rcpp::Named("mz") = mzs);
        if (hasMob)
        {
            df["mobility"] = mobilities;
            // HACK: above assignment converts df to a list: https://stackoverflow.com/a/59369233
            df = Rcpp::DataFrame(df);
        }
        
        ret[i] = df;
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
    
    if (std::find(cn.begin(), cn.end(), "subScan") == cn.end()) // non-IMS
    {
        meta.second.scans = std::move(R_scans);
        meta.second.times = std::move(R_times);
        meta.second.TICs = std::move(R_TICs);
        meta.second.BPCs = std::move(R_BPCs);
        
        for (size_t i=0; i<R_isoStarts.size(); ++i)
            meta.second.isolationRanges.emplace_back(R_isoStarts[i], R_isoEnds[i]);
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
                    meta.second.MSMSFrames.push_back(std::move(curFI));
                    curFI.clear();
                }
            }
            curFI.isolationRanges.emplace_back(R_isoStarts[i], R_isoEnds[i]);
            curFI.subScans.push_back(R_subScans[i]);
            if (!R_subScanEnds.empty())
                curFI.subScanEnds.push_back(R_subScanEnds[i]);
        }
        if (!curFI.empty())
            meta.second.MSMSFrames.push_back(std::move(curFI));
    }
    
    backend.setSpecMetadata(std::move(meta));
}

// [[Rcpp::export]]
Rcpp::List getMSPeakLists(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Time> &startTimes,
                          const std::vector<SpectrumRawTypes::Time> &endTimes,
                          const std::vector<SpectrumRawTypes::Mass> &precursorMZs, bool withPrecursor, int MSLevel,
                          SpectrumRawTypes::Mass isoWindow, const std::string &method,
                          SpectrumRawTypes::Mass mzWindow, const std::vector<SpectrumRawTypes::Mobility> startMobs,
                          const std::vector<SpectrumRawTypes::Mobility> endMobs, unsigned minAbundance, unsigned topMost,
                          SpectrumRawTypes::Intensity minIntensityIMS, SpectrumRawTypes::Intensity minIntensityPre,
                          SpectrumRawTypes::Intensity minIntensityPost)
{
    // UNDONE: add mobility column to output?
    
    const auto entries = startTimes.size();
    const auto clMethod = clustMethodFromStr(method);
    const auto specMeta = backend.getSpecMetadata();
    const auto specFilter = SpectrumRawFilter().setMinIntensity(minIntensityPre).setTopMost(topMost);
    const auto specFilterIMS = SpectrumRawFilter()
        .setMinIntensity(minIntensityIMS)
        .setTopMost(topMost)
        .setWithPrecursor(withPrecursor);
    const auto MSLev = (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2;
    
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
                const auto sels = getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                                       MSLev,
                                                       makeNumRange(precursorMZs[i] - isoWindow, precursorMZs[i] + isoWindow));
                std::vector<SpectrumRaw> spectra(sels.size());
                const auto mobRange = makeNumRange(startMobs[i], endMobs[i]);
                for (size_t j=0; j<sels.size(); ++j)
                {
                    auto spec = backend.readSpectrum(tdata, MSLev, sels[j], mobRange);
                    if (!spec.getMobilities().empty())
                    {
                        spec = filterIMSFrame(spec, specFilterIMS, precursorMZs[i]);
                        spec = averageSpectraRaw(spec, frameSubSpecCount(spec), clMethod,
                                                 mzWindow, false, minIntensityPre, minAbundance);
                    }
                    else
                        spec = filterSpectrumRaw(spec, specFilter, precursorMZs[i]);
                    spectra[j] = std::move(spec);
                }
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
