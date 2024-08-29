#include <Rcpp.h>
#include <algorithm>
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
std::vector<std::vector<OutType>> applyMSData(const MSReadBackend &backend, SpectrumRawTypes::MSLevel MSLevel,
                                              const std::vector<std::vector<SpectrumRawSelection>> &scanSels,
                                              FuncType func, Args... args)
{
    /* This function will apply a callback on selected spectra. Multiple sets of spectra selections are supported, and
     * the function is optimized to avoid reading the same spectra more than once in case of overlap between sets.
     * 
     * The code makes the following assumptions:
     * - each set of SpectrumRawSelection objects is unique and sorted by index from low to high
     * - the index of each SpectrumRawSelection in a set is unique
     * - SpectrumRawSelection objects with the same index but different MSMSFrameIndices may occur across sets. Thus,
     *   across sets the full object instead of merely the index should be compared.
     * 
     * As this function is the main driver for eg EICs and peak lists, this function was heavily optimized and therefore
     * more complex.
    */
    
    const auto entries = scanSels.size();
    
    // to compare spectra selections by just the index, i.e. within a selection set
    const auto selCompInd = [](const SpectrumRawSelection &s, size_t i) { return s.index < i; };
    
    // collect all unique indices so we (1) can avoid >1 reads of a spectrum and (2) can easily use an OpenMP for loop below.
    std::vector<size_t> allScanSelInds;
    const size_t nspecs = (MSLevel == SpectrumRawTypes::MSLevel::MS1) ? backend.getSpecMetadata().first.scans.size() :
        backend.getSpecMetadata().second.scans.size();
    for (size_t i=0; i<nspecs; ++i)
    {
        // check if _any_ of the selection sets has this index
        const bool inSels = std::any_of(scanSels.begin(), scanSels.end(), [&](const auto &sels)
        {
            // NOTE: std::binary_search doesn't seem to support comparator with different types (ie the selCompInd lambda)
            // --> use std::lower_bound instead
            const auto it = std::lower_bound(sels.begin(), sels.end(), i, selCompInd);
            return it != sels.end() && it->index == i;
        });
        if (inSels)
            allScanSelInds.push_back(i);
    }
    
    // pre-allocate ret vector
    std::vector<std::vector<OutType>> ret;
    for (size_t i=0; i<entries; ++i)
        ret.emplace_back(scanSels[i].size());
    
    ThreadExceptionHandler exHandler;

    // UNDONE: make num_threads configurable
    #pragma omp parallel num_threads(12)
    {
        auto tdata = backend.getThreadData();
        #pragma omp for
        for (const auto selInd : allScanSelInds)
        {
            // UNDONE: handle mobility range? For now we stick with post-filtering, which makes the mob range argument
            // to readSpectrum currently unused...
            
            /* In case multiple sets contain spectrum selection with index selInd, there is a chance that the
             * SpectrumRawSelection objects still differ (see above). We therefore keep track of the selection of the
             * previous set and if it differs perform a new spectrum read.
             * 
             * UNDONE: this could be optmized further by also storing and re-using spectra before the previous selection set
            */ 
            
            SpectrumRaw curSpec;
            SpectrumRawSelection curSel;
            bool initSel = true;
            
            for (size_t j=0; j<entries; ++j)
            {
                const auto it = std::lower_bound(scanSels[j].begin(), scanSels[j].end(), selInd, selCompInd);
                if (it == scanSels[j].end() || it->index != selInd)
                    continue; // not found
                
                if (initSel || curSel != *it)
                {
                    initSel = false;
                    curSel = *it;
                    curSpec = backend.readSpectrum(tdata, MSLevel, curSel, SpectrumRawTypes::MobilityRange());
                }
                
                // UNDONE: optimization could be further pushed by not using a 2d vector here?
                const auto ind = std::distance(scanSels[j].begin(), it);
                exHandler.run([&] { ret[j][ind] = func(curSpec, curSel, j, args...); });
            }
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
    const std::vector<std::vector<SpectrumRawSelection>> sels(1, std::vector<SpectrumRawSelection>{SpectrumRawSelection(index)});
    
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        return spec;
    };
    
    const auto spectra = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc);
    const auto &spec = spectra[0][0];
    
    if (!spec.getMobilities().empty())
        return Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                       Rcpp::Named("intensity") = spec.getIntensities(),
                                       Rcpp::Named("mobility") = spec.getMobilities());
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                   Rcpp::Named("intensity") = spec.getIntensities());
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
    struct EICPoint
    {
        SpectrumRawTypes::Mass mz = 0.0;
        SpectrumRawTypes::Intensity intensity = 0;
        SpectrumRawTypes::Mobility mobility = 0.0;
    };
    
    const auto EICCount = startMZs.size();
    bool anySpecHasMob = false;
    
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t entry)
    {
        const bool hasMob = spec.hasMobilities();
        EICPoint ret;
        
        size_t startInd = 0;
        
        // NOTE: for non-IMS data assume spec is mz sorted
        
        if (!hasMob)
        {
            // use lower bound to speedup search for first mass
            const auto mzIt = std::lower_bound(spec.getMZs().begin(), spec.getMZs().end(), startMZs[entry]);
            if (mzIt == spec.getMZs().end())
                return ret;
            startInd = std::distance(spec.getMZs().begin(), mzIt);
        }
        else
            anySpecHasMob = true;
        
        for (size_t j=startInd; j<spec.size(); ++j)
        {
            const auto mz = spec.getMZs()[j];
            const auto mob = (hasMob) ? spec.getMobilities()[j] : 0;
            
            if (hasMob)
            {
                if (mz < startMZs[entry] || mz > endMZs[entry])
                    continue;
                if (mob < startMobs[entry] || mob > endMobs[entry])
                    continue; // UNDONE: optimize eg for TIMS data where mobilities are sorted?
            }
            else if (mz > endMZs[entry])
                break; 
            
            const auto inten = spec.getIntensities()[j];
            ret.intensity += inten;
            ret.mz += mz * inten;
            if (hasMob)
                ret.mobility += mob * inten;
        }

        if (ret.intensity > 0)
        {
            // weighted mean
            ret.mz /= ret.intensity;
            if (hasMob)
                ret.mobility /= ret.intensity;
        }
        
        return ret;
    };
    
    const auto &specMeta = backend.getSpecMetadata();
    std::vector<std::vector<SpectrumRawSelection>> scanSels;
    for (size_t i=0; i<EICCount; ++i)
    {
        scanSels.push_back(getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                                SpectrumRawTypes::MSLevel::MS1, SpectrumRawTypes::IsolationRange()));
    }
    
    const auto allEICPoints = applyMSData<EICPoint>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels, sfunc);
    

    if (allEICPoints.empty())
        return Rcpp::List();

    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        // NOTE: assume all MS spectra have or have not IMS data (UNDONE?)
        
        const auto &points = allEICPoints[i];
        std::vector<SpectrumRawTypes::Time> times(points.size());
        std::vector<SpectrumRawTypes::Mass> mzs(points.size());
        std::vector<SpectrumRawTypes::Intensity> intensities(points.size());
        std::vector<SpectrumRawTypes::Mobility> mobilities((anySpecHasMob) ? points.size() : 0);
        for (size_t j=0; j<points.size(); ++j)
        {
            times[j] = specMeta.first.times[scanSels[i][j].index];
            mzs[j] = points[j].mz;
            intensities[j] = points[j].intensity;
            if (anySpecHasMob)
                mobilities[j] = points[j].mobility;
        }
        
        auto df = Rcpp::DataFrame::create(Rcpp::Named("time") = times,
                                          Rcpp::Named("intensity") = intensities,
                                          Rcpp::Named("mz") = mzs);
        if (anySpecHasMob)
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
    const auto MSLev = (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2;
    const auto specMeta = backend.getSpecMetadata();
    const auto specFilter = SpectrumRawFilter().setMinIntensity(minIntensityPre).setTopMost(topMost);
    const auto specFilterIMS = SpectrumRawFilter()
        .setMinIntensity(minIntensityIMS)
        .setTopMost(topMost)
        .setWithPrecursor(withPrecursor);
    
    const auto &sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t e)
    {
        const bool hasMob = spec.hasMobilities();
        
        if (hasMob)
        {
            const auto specf = filterIMSFrame(spec, specFilterIMS, precursorMZs[e],
                                              makeNumRange(startMobs[e], endMobs[e]));
            return averageSpectraRaw(specf, frameSubSpecCount(specf), clMethod, mzWindow, false, minIntensityPre,
                                     minAbundance);
        }
        return filterSpectrumRaw(spec, specFilter, precursorMZs[e]);
    };
    
    std::vector<std::vector<SpectrumRawSelection>> scanSels;
    for (size_t i=0; i<entries; ++i)
    {
        scanSels.push_back(getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]), MSLev,
                                                makeNumRange(precursorMZs[i] - isoWindow, precursorMZs[i] + isoWindow)));
    }
    
    const auto allSpectra = applyMSData<SpectrumRaw>(backend, MSLev, scanSels, sfunc);
    
    std::vector<SpectrumRaw> averagedSpectra(entries);
    #pragma omp parallel for num_threads(12)
    for (size_t i=0; i<entries; ++i)
        averagedSpectra[i] = averageSpectraRaw(allSpectra[i], clMethod, mzWindow, true, minIntensityPost, minAbundance);

    Rcpp::List ret(entries);
    for (size_t i=0; i<entries; ++i)
    {
        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = averagedSpectra[i].getMZs(),
                                         Rcpp::Named("intensity") = averagedSpectra[i].getIntensities());
    }

    return ret;
}
