#include <Rcpp.h>

#include <algorithm>
#include <vector>
#include <chrono>

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

class SimpleTimer
{
    using clock = std::chrono::steady_clock;
    
    clock::time_point startTime;
    bool running, enabled;
    
public:
    SimpleTimer(bool e = true) : running(false), enabled(e) {}
    
    void start(const std::string &msg = "")
    {
        if (!enabled)
            return;
        
        if (!msg.empty())
            std::cerr << msg << ": ";
        startTime = clock::now();
        running = true;
    }
    
    double stop(void)
    {
        if (!running || !enabled)
            return 0.0;
        const auto end = clock::now();
        running = false;
        const double ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - startTime).count();
        
        if (ms >= 1000.0)
            std::cerr << std::fixed << std::setprecision(3) << (ms / 1000.0) << " s\n";
        else
            std::cerr << std::fixed << std::setprecision(3) << ms << " ms\n";
        
        return ms;
    }
};

template<typename OutType, typename FuncType, typename... Args>
std::vector<std::vector<OutType>> applyMSData(const MSReadBackend &backend, SpectrumRawTypes::MSLevel MSLevel,
                                              const std::vector<std::vector<SpectrumRawSelection>> &scanSels,
                                              FuncType func, SpectrumRawTypes::Intensity minIntensityIMS,
                                              SpectrumRawTypes::MSSortType IMSSortType,
                                              std::function<SpectrumRaw(const SpectrumRaw &)> prepFunc = {},
                                              Args... args)
{
    /* This function will apply a callback on selected spectra. Multiple sets of spectra selections are supported, and
     * the function is optimized to avoid reading the same spectra more than once in case of overlap between sets.
     * 
     * The code makes the following assumptions:
     * - each set of SpectrumRawSelection objects is unique and sorted by index from low to high
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

    #pragma omp parallel
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
             * UNDONE: this could be optimized further by also storing and re-using spectra before the previous selection set
            */ 
            
            exHandler.run([&] {
                SpectrumRaw curSpec;
                SpectrumRawSelection curSel;
                bool initSel = true;
                
                for (size_t j=0; j<entries; ++j)
                {
                    auto it = std::lower_bound(scanSels[j].begin(), scanSels[j].end(), selInd, selCompInd);
                    while (it != scanSels[j].end() && it->index == selInd)
                    {
                       if (initSel || curSel != *it)
                        {
                            initSel = false;
                            curSel = *it;
                            curSpec = backend.readSpectrum(tdata, MSLevel, curSel, SpectrumRawTypes::MobilityRange(),
                                                           minIntensityIMS);
                            if (curSpec.hasMobilities()) // NOTE: for non-IMS data we assume data is always mz sorted
                                curSpec.sort(IMSSortType);
                            if (prepFunc)
                                curSpec = prepFunc(curSpec);
                        }
                        
                        // UNDONE: optimization could be further pushed by not using a 2d vector here?
                        const auto ind = std::distance(scanSels[j].begin(), it);
                        ret[j][ind] = func(curSpec, curSel, j, args...);
                        ++it;
                    }
                }
            });
        }
    }
    
    exHandler.reThrow();
    
    return ret;
}

std::pair<std::vector<SpectrumRawTypes::Time>, std::vector<SpectrumRawTypes::Intensity>>
    padEIC(const std::vector<SpectrumRawTypes::Time> &allTimes,
           SpectrumRawTypes::Time startTime, SpectrumRawTypes::Time endTime,
           const std::vector<SpectrumRawTypes::Time> &times,
           const std::vector<SpectrumRawTypes::Intensity> &intensities)
{
    std::vector<SpectrumRawTypes::Time> outTimes;
    std::vector<SpectrumRawTypes::Intensity> outIntensities;
    
    if (allTimes.empty())
        return std::make_pair(outTimes, outIntensities);
    if (times.empty())
    {
        outTimes = { allTimes.front(), allTimes.back() };
        outIntensities.resize(2, 0.0);
        return std::make_pair(outTimes, outIntensities);
    }
    
    // snap start/end ranges to allTimes
    auto startTimeIt = std::lower_bound(allTimes.begin(), allTimes.end(), startTime);
    if (startTimeIt == allTimes.end())
        startTimeIt = std::prev(allTimes.end());
    startTime = *startTimeIt;
    auto endTimeIt = std::lower_bound(startTimeIt, allTimes.end(), endTime);
    if (endTimeIt == allTimes.end() || (endTimeIt != allTimes.begin() && *endTimeIt > endTime))
        --endTimeIt;
    endTime = *endTimeIt;
    
    auto it = times.begin();
    auto allIt = allTimes.begin();
    
    while (it != times.end() && allIt != allTimes.end())
    {
        // corresponding iterator in allTimes
        const auto allItCorresp = std::lower_bound(allIt, allTimes.end(), *it);
        
        // add zero point if previous was missing in EIC and within the RT range
        if (allItCorresp != allTimes.begin())
        {
            const auto aicPrv = std::prev(allItCorresp);
            if ((it == times.begin() || (!compareTol(*(std::prev(it)), *aicPrv) && !compareTol(*aicPrv, outTimes.back()))) &&
                numberGTE(*aicPrv, startTime))
                
            {
                outTimes.push_back(*aicPrv);
                outIntensities.push_back(0.0);
            }
        }
        
        // add EIC point
        outTimes.push_back(*it);
        outIntensities.push_back(intensities[std::distance(times.begin(), it)]);
        
        allIt = std::next(allItCorresp);
        if (allIt == allTimes.end())
            break;
        
        // add zero point if next is missing in EIC
        ++it;
        if ((it == times.end() || !compareTol(*it, *allIt)) && numberLTE(*allIt, endTime))
        {
            outTimes.push_back(*allIt);
            outIntensities.push_back(0.0);
            ++allIt;
        }
    }
    
    // add final point if needed
    if (allIt != allTimes.end() && numberLTE(allTimes.back(), endTime))
    {
        outTimes.push_back(allTimes.back());
        outIntensities.push_back(0.0);
    }
    
    // add start/stop points if needed
    if (outTimes.empty() || outTimes.front() > startTime)
    {
        outTimes.insert(outTimes.begin(), startTime);
        outIntensities.insert(outIntensities.begin(), 0.0);
    }
    if (outTimes.empty() || outTimes.back() < endTime)
    {
        outTimes.push_back(endTime);
        outIntensities.push_back(0.0);
    }
    
    return std::make_pair(outTimes, outIntensities);
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
        .method("setNeedIMS", &MSReadBackend::setNeedIMS)
        .method("open", &MSReadBackend::open)
        .method("close", &MSReadBackend::close)
        .method("getCurrentFile", &MSReadBackend::getCurrentFile)
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

MSReadBackendUnavailable::MSReadBackendUnavailable(const char *n)
{
    Rcpp::stop("Backend '%s' is not available!", n);
}

// [[Rcpp::export]]
bool backendAvailable(const std::string &b)
{
    if (b == "mstoolkit")
    {
#ifndef WITH_MSTK
        return false;
#endif
    }
    else if (b == "opentims")
    {
#ifndef WITH_OTIMS
        return false;
#endif
    }
    
    return true;
}

// [[Rcpp::export]]
int walkSpectra(const MSReadBackend &backend)
{
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t) { return spec.size(); };
    const auto &meta = backend.getSpecMetadata();
    
    std::vector<std::vector<SpectrumRawSelection>> sels(1);
    for (size_t i=0; i<meta.first.scans.size(); ++i)
        sels[0].emplace_back(i);
    
    const auto ret = applyMSData<size_t>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc, 0,
                                         SpectrumRawTypes::MSSortType::NONE);
    return std::accumulate(ret[0].begin(), ret[0].end(), 0);
    
#if 0
    size_t ret = 0;
    #pragma omp parallel
    {
        auto tdata = backend.getThreadData();
        #pragma omp for
        for (size_t i=0; i<50 /*meta.first.scans.size()*/; ++i)
        {
            //#pragma omp critical (StupidNameHere1)
            {
            const auto s = backend.readSpectrum(tdata, SpectrumRawTypes::MSLevel::MS1, SpectrumRawSelection(i), SpectrumRawTypes::MobilityRange());
            #pragma omp atomic
            ret += s.size();
            }
        }
    }
    return ret;
#endif
}

// [[Rcpp::export]]
int countAllPeaks(const MSReadBackend &backend)
{
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t) { return spec.size(); };
    const auto &meta = backend.getSpecMetadata();
    
    std::vector<std::vector<SpectrumRawSelection>> sels(1);
    for (size_t i=0; i<meta.first.scans.size(); ++i)
        sels[0].emplace_back(i);
    const auto peaksMS1 = applyMSData<size_t>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc, 0,
                                              SpectrumRawTypes::MSSortType::NONE);
    const auto totPeaksMS1 = std::accumulate(peaksMS1[0].begin(), peaksMS1[0].end(), 0);

    sels[0].clear();
    for (size_t i=0; i<meta.second.scans.size(); ++i)
    {
        if (meta.second.MSMSFrames.empty())
            sels[0].emplace_back(i);
        else
        {
            for (size_t j=0; j<meta.second.MSMSFrames[i].precursorMZs.size(); ++j)
            {
                SpectrumRawSelection ssel(i);
                ssel.MSMSFrameIndices.push_back(j);
                sels[0].push_back(std::move(ssel));
            }
        }
    }
    const auto peaksMS2 = applyMSData<size_t>(backend, SpectrumRawTypes::MSLevel::MS2, sels, sfunc, 0,
                                              SpectrumRawTypes::MSSortType::NONE);
    const auto totPeaksMS2 = std::accumulate(peaksMS2[0].begin(), peaksMS2[0].end(), 0);
    
    return totPeaksMS1 + totPeaksMS2;
}

// [[Rcpp::export]]
Rcpp::DataFrame getMSSpectrum(const MSReadBackend &backend, int index, int MSLevel, int frameIndex = -1,
                              double minIntensity = 0)
{
    SpectrumRawSelection sel(index);
    if (MSLevel == 2 && frameIndex != -1)
        sel.MSMSFrameIndices.push_back(frameIndex);
    std::vector<std::vector<SpectrumRawSelection>> sels(1, std::vector<SpectrumRawSelection>{sel});
    
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        return spec;
    };
    
    const auto MSLev = (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2;
    const auto spectra = applyMSData<SpectrumRaw>(backend, MSLev, sels, sfunc, minIntensity,
                                                  SpectrumRawTypes::MSSortType::NONE);
    const auto &spec = spectra[0][0];
    
    if (!spec.getMobilities().empty())
        return Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                       Rcpp::Named("intensity") = spec.getIntensities(),
                                       Rcpp::Named("mobility") = spec.getMobilities());
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                   Rcpp::Named("intensity") = spec.getIntensities());
}

// [[Rcpp::export]]
Rcpp::DataFrame getCollapsedFrame(const MSReadBackend &backend, int index, SpectrumRawTypes::Mass mzWindow,
                                  SpectrumRawTypes::Intensity minIntensityIMS,
                                  SpectrumRawTypes::Intensity minIntensityPre,
                                  SpectrumRawTypes::PeakAbundance minAbundanceRel,
                                  SpectrumRawTypes::PeakAbundance minAbundanceAbs, const std::string &method,
                                  double mzMin = 0.0, double mzMax = 0.0, double minInt = 0.0, unsigned topMost = 0,
                                  double prec = 0.0, double mobMin = 0.0, double mobMax = 0.0)
{
    const auto clMethod = clustMethodFromStr(method);
    
    SpectrumRawSelection sel(index);
    std::vector<std::vector<SpectrumRawSelection>> sels(1, std::vector<SpectrumRawSelection>{sel});
    
    SpectrumRawFilter filter;
    if (mzMin > 0.0 || mzMax > 0.0)
        filter.setMZRange(mzMin, mzMax);
    if (minInt > 0.0)
        filter.setMinIntensity(minInt);
    if (topMost > 0)
        filter.setTopMost(topMost);
    
    
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        return filterIMSFrame(spec, filter, prec, SpectrumRawTypes::MobilityRange(mobMin, mobMax));
    };
    
    const auto spectra = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc,
                                                  minIntensityIMS, SpectrumRawTypes::MSSortType::MOBILITY_MZ);
    const auto &spec = averageSpectraRaw(spectra[0][0], frameSubSpecIDs(spectra[0][0]), clMethod, mzWindow, false,
                                         minIntensityPre, minAbundanceRel, minAbundanceAbs);
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                   Rcpp::Named("intensity") = spec.getIntensities(),
                                   Rcpp::Named("abundance_rel") = spec.getAbundancesRel(),
                                   Rcpp::Named("abundance_abs") = spec.getAbundancesAbs());
}

// [[Rcpp::export]]
Rcpp::DataFrame getFilteredFrame(const MSReadBackend &backend, int index, double mzMin = 0.0, double mzMax = 0.0,
                                 double minInt = 0.0, unsigned topMost = 0, double prec = 0.0, double mobMin = 0.0,
                                 double mobMax = 0.0)
{
    SpectrumRawSelection sel(index);
    std::vector<std::vector<SpectrumRawSelection>> sels(1, std::vector<SpectrumRawSelection>{sel});
    
    SpectrumRawFilter filter;
    if (mzMin > 0.0 || mzMax > 0.0)
        filter.setMZRange(mzMin, mzMax);
    if (minInt > 0.0)
        filter.setMinIntensity(minInt);
    if (topMost > 0)
        filter.setTopMost(topMost);
    
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        return filterIMSFrame(spec, filter, prec, SpectrumRawTypes::MobilityRange(mobMin, mobMax));
    };
    
    const auto spec = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc,
                                                  0, SpectrumRawTypes::MSSortType::MOBILITY_MZ)[0][0];

    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                   Rcpp::Named("intensity") = spec.getIntensities(),
                                   Rcpp::Named("mobility") = spec.getMobilities());
}

// [[Rcpp::export]]
Rcpp::DataFrame getCentroidedFrame(const MSReadBackend &backend, int index, SpectrumRawTypes::Mass mzWindow,
                                   SpectrumRawTypes::Mobility mobWindow, SpectrumRawTypes::Intensity minIntensity,
                                   const std::string &method)
{
    const auto clMethod = clustMethodFromStr(method);
    
    SpectrumRawSelection sel(index);
    std::vector<std::vector<SpectrumRawSelection>> sels(1, std::vector<SpectrumRawSelection>{sel});
    
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        return spec;
    };
    
    const auto spectra = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc, 0,
                                                  SpectrumRawTypes::MSSortType::MOBILITY_MZ);
    const auto &spec = centroidIMSFrame(spectra[0][0], clMethod, mzWindow, mobWindow, minIntensity);
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                   Rcpp::Named("intensity") = spec.getIntensities(),
                                   Rcpp::Named("mobility") = spec.getMobilities());
}

// [[Rcpp::export]]
Rcpp::DataFrame getScans(const MSReadBackend &backend, SpectrumRawTypes::Time timeStart, SpectrumRawTypes::Time timeEnd,
                         int MSLevel, SpectrumRawTypes::Mass prec, SpectrumRawTypes::Mass fixedIsoWidth = 0.0)
{
    const auto sels = getSpecRawSelections(backend.getSpecMetadata(), makeNumRange(timeStart, timeEnd),
                                           (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2,
                                           prec, fixedIsoWidth);
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

#if 0
// older version that loads EICs directly in applyMSData(); pros: doesn't need to load full spectra, cons: may more use
// more memory for many EICs

Rcpp::List getEICList(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Mass> &startMZs,
                      const std::vector<SpectrumRawTypes::Mass> &endMZs,
                      const std::vector<SpectrumRawTypes::Time> &startTimes,
                      const std::vector<SpectrumRawTypes::Time> &endTimes,
                      const std::vector<SpectrumRawTypes::Mobility> &startMobs,
                      const std::vector<SpectrumRawTypes::Mobility> &endMobs,
                      SpectrumRawTypes::Mass mzExpIMSWindow, SpectrumRawTypes::Intensity minIntensityIMS, bool compress,
                      bool showProgress = false, bool withBP = false, SpectrumRawTypes::Intensity minEICIntensity = 0,
                      SpectrumRawTypes::Time minEICAdjTime = 0,
                      SpectrumRawTypes::Intensity minEICAdjIntensity = 0)
{
    struct EICPoint
    {
        SpectrumRawTypes::Time time = 0.0;
        SpectrumRawTypes::Mass mz = 0.0, mzMin = 0.0, mzMax = 0.0;
        SpectrumRawTypes::Intensity intensity = 0;
        SpectrumRawTypes::Mobility mobility = 0.0, mobMin = 0.0, mobMax = 0.0;
        SpectrumRawTypes::Mass mzBP = 0.0;
        SpectrumRawTypes::Mobility mobilityBP = 0.0;
        SpectrumRawTypes::Intensity intensityBP = 0.0;
    };
    
    const auto EICCount = startMZs.size();
    const auto &specMeta = backend.getSpecMetadata();
    bool anySpecHasMob = false;
    
    const auto minMZ = *(std::min_element(startMZs.begin(), startMZs.end()));
    const auto maxMZ = *(std::max_element(endMZs.begin(), endMZs.end()));
    const auto minMob = (startMobs.empty()) ? 0 : *(std::min_element(startMobs.begin(), startMobs.end()));
    const auto maxMob = (endMobs.empty()) ? 0 : *(std::max_element(endMobs.begin(), endMobs.end()));
    
    const auto specPrepFunc = [&](const SpectrumRaw &spec)
    {
        if (spec.hasMobilities())
        {
            const auto sInds = getSortedInds(spec.getMZs());
            SpectrumRaw sortedSpec;
            for (size_t k=0; k<spec.size(); ++k)
            {
                if (spec.getMZs()[sInds[k]] < minMZ || spec.getMZs()[sInds[k]] > maxMZ)
                    continue;
                if (spec.getMobilities()[sInds[k]] < minMob ||
                    (maxMob != 0.0 && spec.getMobilities()[sInds[k]] > maxMob))
                    continue;
                    
                sortedSpec.append(spec.getMZs()[sInds[k]], spec.getIntensities()[sInds[k]],
                                  spec.getMobilities()[sInds[k]]);
            }
            return sortedSpec;
        }
        return spec;
    };
    
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t entry)
    {
        const bool hasMob = spec.hasMobilities();
        EICPoint ret;
        
        ret.time = specMeta.first.times[ssel.index];

        // NOTE: assume spec is mz sorted
        // use lower bound to speedup search for first mass
        const auto mzStart = startMZs[entry] - ((hasMob) ? mzExpIMSWindow : 0.0);
        const auto mzEnd = endMZs[entry] + ((hasMob) ? mzExpIMSWindow : 0.0);
        const auto mzIt = std::lower_bound(spec.getMZs().begin(), spec.getMZs().end(), mzStart);
        if (mzIt == spec.getMZs().end())
            return ret;
        const auto startInd = std::distance(spec.getMZs().begin(), mzIt);
        
        if (hasMob)
            anySpecHasMob = true;
        
        for (size_t j=startInd; j<spec.size(); ++j)
        {
            const auto mz = spec.getMZs()[j];
            if (mz > mzEnd)
                break; 
            
            const auto mob = (hasMob) ? spec.getMobilities()[j] : 0;
            if (hasMob)
            {
                if (mob < startMobs[entry] || (endMobs[entry] != 0.0 && mob > endMobs[entry]))
                    continue;
            }
            
            const auto inten = spec.getIntensities()[j];
            if (inten > 0)
            {
                ret.intensity += inten;
                ret.mz += mz * inten;
                if (ret.mzMin == 0.0 || mz < ret.mzMin)
                    ret.mzMin = mz;
                if (mz > ret.mzMax)
                    ret.mzMax = mz;
                if (hasMob)
                {
                    ret.mobility += mob * inten;
                    if (ret.mobMin == 0.0 || mob < ret.mobMin)
                        ret.mobMin = mob;
                    if (mob > ret.mobMax)
                        ret.mobMax = mob;
                }
                if (withBP && inten > ret.intensityBP)
                {
                    ret.intensityBP = inten;
                    ret.mzBP = mz;
                    ret.mobilityBP = mob;
                }
                //Rcpp::Rcout << "EIC: " << ret.time << "/" << j << "\t" << std::fixed << mz << "\t" << inten << "\t" << ret.mzMin << "/" << ret.mzMax << "\t" << mob << "\n";
            }
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
    
    std::vector<std::vector<SpectrumRawSelection>> scanSels;
    for (size_t i=0; i<EICCount; ++i)
    {
        scanSels.push_back(getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                                SpectrumRawTypes::MSLevel::MS1, 0));
    }
    
    auto allEICPoints = applyMSData<EICPoint>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels, sfunc, minIntensityIMS,
                                              showProgress, specPrepFunc);

    if (allEICPoints.empty())
        return Rcpp::List();

    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        // NOTE: assume all MS spectra have or have not IMS data (UNDONE?)
        
        auto &points = allEICPoints[i];
        
        // start with some pre-filtering
        bool valid = true;
        
        if (minEICIntensity > 0)
        {
            const auto maxIt = std::max_element(points.begin(), points.end(), [](const auto &a, const auto &b)
            {
                return a.intensity < b.intensity;
            });
            valid = numberGTE(maxIt->intensity, minEICIntensity);
        }
        if (valid && minEICAdjTime > 0 && minEICAdjIntensity > 0)
        {
            SpectrumRawTypes::Time startTime = 0;
            valid = false; // until proven otherwise
            for (size_t j=0; j<points.size(); ++j)
            {
                if (numberGTE(points[j].intensity, minEICAdjIntensity))
                {
                    if (startTime == 0)
                        startTime = points[j].time;
                    else if (numberGTE(points[j].time - startTime, minEICAdjTime))
                    {
                        valid = true;
                        break;
                    }
                }
                else
                    startTime = 0;
            }
        }
        
        if (!valid)
        {
            ret[i] = Rcpp::List();
            continue;
        }
        
        if (compress && points.size() >= 3)
        {
            for (auto it=std::next(points.begin()); it!=std::prev(points.end()); )
            {
                if (it->intensity == 0)
                {
                    if (std::prev(it)->intensity == 0 && std::next(it)->intensity == 0)
                    {
                        it = points.erase(it);
                        continue;
                    }
                }
                ++it;
            }
        }

        std::vector<SpectrumRawTypes::Time> times(points.size());
        std::vector<SpectrumRawTypes::Mass> mzs(points.size()), mzMins(points.size()), mzMaxs(points.size());
        std::vector<SpectrumRawTypes::Intensity> intensities(points.size());
        std::vector<SpectrumRawTypes::Mobility> mobilities((anySpecHasMob) ? points.size() : 0);
        std::vector<SpectrumRawTypes::Mobility> mobMins((anySpecHasMob) ? points.size() : 0);
        std::vector<SpectrumRawTypes::Mobility> mobMaxs((anySpecHasMob) ? points.size() : 0);
        std::vector<SpectrumRawTypes::Intensity> intensitiesBP((withBP) ? points.size() : 0);
        std::vector<SpectrumRawTypes::Mass> mzsBP((withBP) ? points.size() : 0);
        std::vector<SpectrumRawTypes::Mass> mobilitiesBP((anySpecHasMob && withBP) ? points.size() : 0);
        for (size_t j=0; j<points.size(); ++j)
        {
            times[j] = points[j].time;
            mzs[j] = points[j].mz; mzMins[j] = points[j].mzMin; mzMaxs[j] = points[j].mzMax;
            intensities[j] = points[j].intensity;
            if (anySpecHasMob)
            {
                mobilities[j] = points[j].mobility;
                mobMins[j] = points[j].mobMin;
                mobMaxs[j] = points[j].mobMax;
            }
            if (withBP)
            {
                intensitiesBP[j] = points[j].intensityBP;
                mzsBP[j] = points[j].mzBP;
                if (anySpecHasMob)
                    mobilitiesBP[j] = points[j].mobilityBP;
            }
        }
        
        auto df = Rcpp::List::create(Rcpp::Named("time") = times,
                                     Rcpp::Named("intensity") = intensities,
                                     Rcpp::Named("mz") = mzs,
                                     Rcpp::Named("mzmin") = mzMins,
                                     Rcpp::Named("mzmax") = mzMaxs);
        
        if (anySpecHasMob)
        {
            df["mobility"] = mobilities;
            df["mobmin"] = mobMins;
            df["mobmax"] = mobMaxs;
        }
        
        if (withBP)
        {
            df["intensityBP"] = intensitiesBP;
            df["mzBP"] = mzsBP;
            if (anySpecHasMob)
                df["mobilityBP"] = mobilitiesBP;
        }

        ret[i] = df;
    }
    
    return ret;
}
#endif

// [[Rcpp::export]]
Rcpp::List getEICList(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Mass> &startMZs,
                      const std::vector<SpectrumRawTypes::Mass> &endMZs,
                      const std::vector<SpectrumRawTypes::Time> &startTimes,
                      const std::vector<SpectrumRawTypes::Time> &endTimes,
                      const std::vector<SpectrumRawTypes::Mobility> &startMobs,
                      const std::vector<SpectrumRawTypes::Mobility> &endMobs,
                      SpectrumRawTypes::Time gapFactor, SpectrumRawTypes::Mass mzExpIMSWindow,
                      SpectrumRawTypes::Intensity minIntensityIMS, const std::string &mode = "simple",
                      bool pad = false, SpectrumRawTypes::Intensity minEICIntensity = 0,
                      SpectrumRawTypes::Time minEICAdjTime = 0, unsigned minEICAdjPoints = 0,
                      SpectrumRawTypes::Intensity minEICAdjIntensity = 0, unsigned topMost = 0)
{
    // NOTE: startTimes/endTimes may be length one vectors, in which case they are used for all EICs
    
    struct EICPoint
    {
        SpectrumRawTypes::Time time = 0.0;
        SpectrumRawTypes::Mass mz = 0.0, mzMin = 0.0, mzMax = 0.0;
        SpectrumRawTypes::Intensity intensity = 0;
        SpectrumRawTypes::Mobility mobility = 0.0, mobMin = 0.0, mobMax = 0.0;
        SpectrumRawTypes::Mass mzBP = 0.0;
        SpectrumRawTypes::Mobility mobilityBP = 0.0;
        SpectrumRawTypes::Intensity intensityBP = 0.0;
    };
    
    struct EIC
    {
        std::vector<SpectrumRawTypes::Time> times;
        std::vector<SpectrumRawTypes::Mass> mzs, mzMins, mzMaxs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        std::vector<SpectrumRawTypes::Mobility> mobilities, mobMins, mobMaxs;
        std::vector<SpectrumRawTypes::Mass> mzsBP;
        std::vector<SpectrumRawTypes::Mobility> mobilitiesBP;
        std::vector<SpectrumRawTypes::Intensity> intensitiesBP;
        void addPoint(const EICPoint &p)
        {
            times.push_back(p.time);
            mzs.push_back(p.mz);
            mzMins.push_back(p.mzMin);
            mzMaxs.push_back(p.mzMax);
            intensities.push_back(p.intensity);
            mobilities.push_back(p.mobility);
            mobMins.push_back(p.mobMin);
            mobMaxs.push_back(p.mobMax);
            mzsBP.push_back(p.mzBP);
            mobilitiesBP.push_back(p.mobilityBP);
            intensitiesBP.push_back(p.intensityBP);
        }
        void addPoint(SpectrumRawTypes::Time t, SpectrumRawTypes::Intensity i) { times.push_back(t); intensities.push_back(i); }
        void clear(void)
        {
            times.clear();
            mzs.clear();
            mzMins.clear();
            mzMaxs.clear();
            intensities.clear();
            mobilities.clear();
            mobMins.clear();
            mobMaxs.clear();
            mzsBP.clear();
            mobilitiesBP.clear();
            intensitiesBP.clear();
        }
        size_t size(void) const { return times.size(); }
        bool empty(void) const { return times.empty(); }
    };
    
    struct AllPeaks
    {
        std::vector<SpectrumRawTypes::Scan> indices;
        std::vector<SpectrumRawTypes::Mass> mzs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        AllPeaks() = default;
        AllPeaks(size_t size, bool withMob) : indices(size), mzs(size), intensities(size), mobilities((withMob) ? size : 0) { }
        void clear(void)
        {
            indices.clear();
            mzs.clear();
            intensities.clear();
            mobilities.clear();
        }
    };
    
    const auto EICCount = startMZs.size();
    if (EICCount == 0)
        return Rcpp::List();
    
    enum class EICMode { SIMPLE, FULL, FULL_MZ, TEST };
    const auto eicMode = (mode == "simple") ? EICMode::SIMPLE : (mode == "full") ? EICMode::FULL : (mode == "full_mz") ? EICMode::FULL_MZ : EICMode::TEST;
    
    const auto &specMeta = backend.getSpecMetadata();
    
    const auto minMZ = *(std::min_element(startMZs.begin(), startMZs.end()));
    const auto maxMZ = *(std::max_element(endMZs.begin(), endMZs.end()));
    const auto minMob = (startMobs.empty()) ? 0 : *(std::min_element(startMobs.begin(), startMobs.end()));
    const auto maxMob = (endMobs.empty()) ? 0 : *(std::max_element(endMobs.begin(), endMobs.end()));
    
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t)
    {
        if (minMZ == 0.0 && maxMZ == 0.0 && minMob == 0.0 && maxMob == 0.0)
            return spec;
        
        const bool hasMob = spec.hasMobilities();
        
        // NOTE: specs are mz sorted
        const auto startIt = std::lower_bound(spec.getMZs().begin(), spec.getMZs().end(), minMZ);
        
        SpectrumRaw specf;
        for (size_t i=std::distance(spec.getMZs().begin(), startIt); i<spec.size(); ++i)
        {
            if (spec.getMZs()[i] > maxMZ)
                break;
            
            if (hasMob)
            {
                if (spec.getMobilities()[i] < minMob || (maxMob != 0.0 && spec.getMobilities()[i] > maxMob))
                    continue;
                specf.append(spec.getMZs()[i], spec.getIntensities()[i], spec.getMobilities()[i]);
            }
            else
                specf.append(spec.getMZs()[i], spec.getIntensities()[i]);
        }
        return specf;
    };
    
    SimpleTimer timer(true);
    
    timer.start("Collecting all scans");
    std::set<SpectrumRawTypes::Scan> allScans;
    for (size_t i=0; i<startTimes.size(); ++i)
    {
        const auto sels = getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                               SpectrumRawTypes::MSLevel::MS1, 0);
        for (const auto &sel : sels)
            allScans.insert(sel.index);
    }
    std::vector<std::vector<SpectrumRawSelection>> scanSels(1);
    for (const auto &scan : allScans)
        scanSels[0].emplace_back(scan);
    timer.stop();
    
    timer.start("Loading all MS spectra");
    auto allSpectra = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels, sfunc,
                                               minIntensityIMS, SpectrumRawTypes::MSSortType::MZ);
    timer.stop();
    if (allSpectra.empty())
        return Rcpp::List();
    
    timer.start("Preparing EICs (AllPeaks)");
    const auto totalPeaks = std::accumulate(allSpectra[0].begin(), allSpectra[0].end(), 0,
        [](size_t sum, const SpectrumRaw &spec) { return sum + spec.size(); });
    
    const bool anySpecHasMob = allSpectra[0][0].hasMobilities(); // NOTE: assume all MS spectra have or have not IMS data (UNDONE?)
    AllPeaks allPeaks(totalPeaks, anySpecHasMob);
    for (size_t i=0, peaksi=0; i<allSpectra[0].size(); ++i)
    {
        if (allSpectra[0][i].getMZs().empty())
            continue; // no peaks in this spectrum
      
        std::fill(allPeaks.indices.begin() + peaksi, allPeaks.indices.begin() + peaksi + allSpectra[0][i].size(),
                  scanSels[0][i].index);
        std::copy(allSpectra[0][i].getMZs().begin(), allSpectra[0][i].getMZs().end(),
                  allPeaks.mzs.begin() + peaksi);
        std::copy(allSpectra[0][i].getIntensities().begin(), allSpectra[0][i].getIntensities().end(),
                  allPeaks.intensities.begin() + peaksi);
        if (anySpecHasMob)
        {
            std::copy(allSpectra[0][i].getMobilities().begin(), allSpectra[0][i].getMobilities().end(),
                      allPeaks.mobilities.begin() + peaksi);
        }
        
        peaksi += allSpectra[0][i].size();
        allSpectra[0][i].clear();
    }
    timer.stop();

    timer.start("Sorting all peaks");
    const auto sortedInds = getSortedInds(allPeaks.mzs);

    AllPeaks allPeaksSorted(allPeaks.indices.size(), anySpecHasMob);
    for (size_t i=0; i<sortedInds.size(); ++i)
    {
        const auto j = sortedInds[i];
        allPeaksSorted.indices[i] = allPeaks.indices[j];
        allPeaksSorted.mzs[i] = allPeaks.mzs[j];
        allPeaksSorted.intensities[i] = allPeaks.intensities[j];
        if (anySpecHasMob)
            allPeaksSorted.mobilities[i] = allPeaks.mobilities[j];
    }
    allPeaks.clear();
    timer.stop();

    timer.start("Processing EICs");
    std::vector<EIC> allEICs(EICCount);
    std::vector<SpectrumRawTypes::Intensity> allEICMaxIntensities(EICCount);
    #pragma omp parallel for
    for (size_t i=0; i<EICCount; ++i)
    {
        EIC eic;
        SpectrumRawTypes::Intensity maxInten = 0;
        SpectrumRawTypes::Time startTimeAboveThr = 0;
        unsigned adjPointsAboveThr = 0;
        bool enoughTimeAboveThr = false, enoughPointsAboveThr = false;
     
        const auto mzStart = startMZs[i] - ((anySpecHasMob) ? mzExpIMSWindow : 0.0);
        const auto mzEnd = endMZs[i] + ((anySpecHasMob) ? mzExpIMSWindow : 0.0);
        const auto timeStart = (startTimes.size() == 1) ? startTimes[0] : startTimes[i];
        const auto timeEnd = (endTimes.size() == 1) ? endTimes[0] : endTimes[i];
     
        const auto itStart = std::lower_bound(allPeaksSorted.mzs.cbegin(), allPeaksSorted.mzs.cend(), mzStart);
        if (itStart == allPeaksSorted.mzs.cend())
            continue;
        const auto itEnd = std::prev(std::upper_bound(itStart, allPeaksSorted.mzs.cend(), mzEnd));
        const auto startInd = std::distance(allPeaksSorted.mzs.cbegin(), itStart);
        auto endInd = std::distance(allPeaksSorted.mzs.cbegin(), itEnd);
        if (startInd > endInd)
            endInd = startInd; // only one peak at the end
        const auto sortedInds = getSortedInds(allPeaksSorted.indices.cbegin() + startInd,
                                              allPeaksSorted.indices.cbegin() + endInd);

        EICPoint curPoint;
        size_t curScanInd = allPeaksSorted.indices.size(), prvScanInd = allPeaksSorted.indices.size();
        bool init = true;
        for (size_t j=0; ; ++j)
        {
            const bool ended = j == sortedInds.size();
            if (ended && init)
                break; // nothing was done
            
            const auto allPeaksSortedInd = (ended) ? 0 : (startInd + sortedInds[j]);
            const auto scanInd = allPeaksSorted.indices[allPeaksSortedInd];
            // Rcpp::Rcout << "EIC: " << j << "/" << sortedInds.size() << "/" << specMeta.first.scans[curScanInd] << "/" << curScanInd << "/" << init << "/" << ended << "/" << time << "/" << startInd << "/" << endInd << "\n";
            
            if (ended || (!init && scanInd != curScanInd))
            {
                const auto curTime = specMeta.first.times[curScanInd];
                if (curPoint.intensity > 0)
                {
                    curPoint.time = curTime;
                    if (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ)
                    {
                        curPoint.mz /= curPoint.intensity;
                        if (anySpecHasMob && eicMode != EICMode::FULL_MZ)
                            curPoint.mobility /= curPoint.intensity;
                        eic.addPoint(curPoint);
                    }
                    else
                        eic.addPoint(curPoint.time, curPoint.intensity);
                    if (curPoint.intensity > maxInten)
                        maxInten = curPoint.intensity;
                }
                
                if ((!enoughTimeAboveThr && minEICAdjTime > 0.0) || (!enoughPointsAboveThr && minEICAdjPoints > 0))
                {
                    if (prvScanInd != (curScanInd - 1))
                    {
                        // Rcpp::Rcout << "Stop! EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
                        startTimeAboveThr = 0.0;
                        adjPointsAboveThr = 0;   
                    }
                    
                    if (prvScanInd != allPeaksSorted.indices.size() &&
                        numberGTE(curPoint.intensity, minEICAdjIntensity))
                    {
                        if (minEICAdjTime > 0.0)
                        {
                            if (startTimeAboveThr == 0.0)
                                startTimeAboveThr = curTime;
                            else if (numberGTE(curTime - startTimeAboveThr, minEICAdjTime))
                                enoughTimeAboveThr = true;
                        }
                        if (minEICAdjPoints > 0)
                        {
                            ++adjPointsAboveThr;
                            // Rcpp::Rcout << "above: EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
                            if (adjPointsAboveThr >= minEICAdjPoints)
                                enoughPointsAboveThr = true;
                        }
                    }
                }
                
                if (ended)
                    break;
                
                curPoint = EICPoint();
                prvScanInd = curScanInd;
                curScanInd = scanInd;
            }

            const auto time = specMeta.first.times[scanInd];
            if (!timeGTE(time, timeStart) || (timeEnd != 0.0 && !timeLTE(time, timeEnd)))
                continue;
            
            //Rcpp::Rcout << "EIC: " << allPeaksSortedInd << "/" << j << "/" << specMeta.first.times[allPeaksSorted.indices[allPeaksSortedInd]] << "/" << mz << "/" << mzStart << "/" << mzEnd << "\n";
            const auto mob = (anySpecHasMob) ? allPeaksSorted.mobilities[allPeaksSortedInd] : 0;
            if (anySpecHasMob)
            {
                if (mob < startMobs[i] || (endMobs[i] != 0.0 && mob > endMobs[i]))
                    continue;
            }
            
            if (init)
            {
                curScanInd = scanInd;
                init = false;
            }
            
            const auto mz = allPeaksSorted.mzs[allPeaksSortedInd];
            const auto inten = allPeaksSorted.intensities[allPeaksSortedInd];
            if (inten > 0)
            {
                curPoint.intensity += inten;
                if (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ)
                {
                    curPoint.mz += mz * inten;
                    if (curPoint.mzMin == 0.0 || mz < curPoint.mzMin)
                        curPoint.mzMin = mz;
                    if (mz > curPoint.mzMax)
                        curPoint.mzMax = mz;
                    if (anySpecHasMob && eicMode != EICMode::FULL_MZ)
                    {
                        curPoint.mobility += mob * inten;
                        if (curPoint.mobMin == 0.0 || mob < curPoint.mobMin)
                            curPoint.mobMin = mob;
                        if (mob > curPoint.mobMax)
                            curPoint.mobMax = mob;
                    }
                    if (inten > curPoint.intensityBP)
                    {
                        curPoint.intensityBP = inten;
                        curPoint.mzBP = mz;
                        curPoint.mobilityBP = mob;
                    }
                    //Rcpp::Rcout << "EIC: " << curPoint.time << "/" << j << "\t" << std::fixed << mz << "\t" << inten << "\t" << curPoint.mzMin << "/" << curPoint.mzMax << "\t" << mob << "\n";
                }
            }
        }
        
        if (minEICIntensity != 0.0 && maxInten < minEICIntensity)
            continue;
        if (!enoughTimeAboveThr && minEICAdjTime > 0.0)
            continue;
        if (!enoughPointsAboveThr && minEICAdjPoints > 0)
            continue;
        
        if (eicMode == EICMode::SIMPLE && pad && !eic.empty())
        {
            const auto padded = padEIC(specMeta.first.times, timeStart, timeEnd, eic.times, eic.intensities);
            eic.times = std::move(padded.first);
            eic.intensities = std::move(padded.second);
        }
        
        allEICs[i] = std::move(eic);
        allEICMaxIntensities[i] = maxInten;
    }
    
    allPeaksSorted.clear(); // no need for this anymore
    timer.stop();
    
    SpectrumRawTypes::Intensity minMaxIntens = 0.0;
    if (topMost > 0 && allEICMaxIntensities.size() > topMost)
    {
        timer.start("Determining EIC topMost...");
        auto sortedAEMI = allEICMaxIntensities;
        std::nth_element(sortedAEMI.begin(), sortedAEMI.begin() + (topMost-1), sortedAEMI.end(),
                         std::greater<SpectrumRawTypes::Intensity>());
        minMaxIntens = sortedAEMI[topMost-1];
        timer.stop();
    }
    
    timer.start("Converting EICs to R data");
    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        auto &eic = allEICs[i];
        
        if (allEICMaxIntensities[i] < minMaxIntens)
        {
            // just clear the EIC so an empty will be added below. Note that we would clear the EIC at the end anyway.
            eic.clear();
        }
        
        // NOTE: constructing the ouput list directly with create() is much faster than assigning items afterwards.
        // Hence, we do quite some if statements here...
    
        if (eicMode == EICMode::SIMPLE)
        {
            auto mat = Rcpp::NumericMatrix(eic.times.size(), 2);
            for (size_t i=0; i<eic.times.size(); ++i)
            {
                mat(i, 0) = eic.times[i];
                mat(i, 1) = eic.intensities[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity");
            ret[i] = mat;
        }
        else if (eicMode == EICMode::FULL && anySpecHasMob)
        {
            auto mat = Rcpp::NumericMatrix(eic.times.size(), 11);
            for (size_t i=0; i<eic.times.size(); ++i)
            {
                mat(i, 0) = eic.times[i];
                mat(i, 1) = eic.intensities[i];
                mat(i, 2) = eic.intensitiesBP[i];
                mat(i, 3) = eic.mzs[i];
                mat(i, 4) = eic.mzsBP[i];
                mat(i, 5) = eic.mzMins[i];
                mat(i, 6) = eic.mzMaxs[i];
                mat(i, 7) = eic.mobilities[i];
                mat(i, 8) = eic.mobMins[i];
                mat(i, 9) = eic.mobMaxs[i];
                mat(i, 10) = eic.mobilitiesBP[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity", "intensityBP", "mz", "mzBP",
                                                                "mzmin", "mzmax", "mobility", "mobmin", "mobmax",
                                                                "mobilityBP");
            ret[i] = mat;
        }
        else if (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ)
        {
            auto mat = Rcpp::NumericMatrix(eic.times.size(), 7);
            for (size_t i=0; i<eic.times.size(); ++i)
            {
                mat(i, 0) = eic.times[i];
                mat(i, 1) = eic.intensities[i];
                mat(i, 2) = eic.intensitiesBP[i];
                mat(i, 3) = eic.mzs[i];
                mat(i, 4) = eic.mzsBP[i];
                mat(i, 5) = eic.mzMins[i];
                mat(i, 6) = eic.mzMaxs[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity", "intensityBP", "mz", "mzBP",
                                                                "mzmin", "mzmax");
            ret[i] = mat;
        }
        else // if (eicMode == EICMode::TEST)
        {
            ret[i] = !eic.empty();
        }
        
        eic.clear(); // free memory as EICs may consume a lot
    }
    
    if (eicMode != EICMode::TEST)
    {
        if (specMeta.first.times.size() > 1 && gapFactor > 0.0)
        {
            // NOTE: Bruker TIMS data (and maybe others?) seem to omit zero intensity scans, leading to time gaps. Since
            // this leads to incorrect EICs, we pad here. To detect gaps, we take the median difference between scans
            // and assume a gap is at least X higher than that.For padding we just have to add additional time points,
            // EIC filling and padding functions will assume these are zero intensity points.
            std::vector<SpectrumRawTypes::Time> diffs(specMeta.first.times.size() - 1, 0);
            for (size_t i=1; i<specMeta.first.times.size(); ++i)
                diffs[i - 1] = specMeta.first.times[i] - specMeta.first.times[i - 1];
            const double medianDiff = median(diffs);
            std::vector<SpectrumRawTypes::Time> allXValues;
            for (size_t i=0; i<specMeta.first.times.size(); ++i)
            {
                allXValues.push_back(specMeta.first.times[i]);
                if (i == (specMeta.first.times.size() - 1))
                    break;
                
                const auto diff = specMeta.first.times[i + 1] - specMeta.first.times[i];
                if (diff > medianDiff * gapFactor) // add dummy point after current
                {
                    allXValues.push_back(specMeta.first.times[i] + medianDiff);
                    if (diff > medianDiff * 2.0) // add dummy point before next
                        allXValues.push_back(specMeta.first.times[i + 1] - medianDiff);
                }
            }
            ret.attr("allXValues") = allXValues;
        }
        else
            ret.attr("allXValues") = specMeta.first.times;
    }
    timer.stop();
    
    return ret;
}

// [[Rcpp::export]]
std::vector<SpectrumRawTypes::Intensity> doFillEIXIntensities(const std::vector<SpectrumRawTypes::Time> &allXValues,
                                                              const std::vector<SpectrumRawTypes::Time> &xvalues,
                                                              const std::vector<SpectrumRawTypes::Intensity> &intensities)
{
    return fillEIXIntensities(allXValues, xvalues, intensities);
}

// [[Rcpp::export]]
Rcpp::DataFrame getMSMetadata(const MSReadBackend &backend, int msLevel)
{
    const auto &meta = backend.getSpecMetadata();
    
    const auto polsToInts = [](const auto &v)
    {
        std::vector<int> ret(v.size());
        for (size_t i=0; i<v.size(); ++i)
            ret[i] = static_cast<int>(v[i]);
        return ret;
    };
    
    if (msLevel == 1)
        return Rcpp::DataFrame::create(Rcpp::Named("scan") = meta.first.scans,
                                       Rcpp::Named("time") = meta.first.times,
                                       Rcpp::Named("TIC") = meta.first.TICs,
                                       Rcpp::Named("BPC") = meta.first.BPCs,
                                       Rcpp::Named("polarity") = polsToInts(meta.first.polarities));
    
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
                                       Rcpp::Named("polarity") = polsToInts(meta.second.polarities),
                                       Rcpp::Named("isolationRangeMin") = isolationMins,
                                       Rcpp::Named("isolationRangeMax") = isolationMaxs,
                                       Rcpp::Named("precursorMZ") = meta.second.precursorMZs);
    }
    
    // MS2 / IMS --> convert to tidy format
    
    std::vector<SpectrumRawTypes::Scan> scans;
    std::vector<SpectrumRawTypes::Time> times;
    std::vector<SpectrumRawTypes::Intensity> TICs, BPCs;
    std::vector<int> polaritiesInt;
    std::vector<SpectrumRawTypes::Mass> isolationMins, isolationMaxs;
    std::vector<SpectrumRawTypes::Mass> precursorMZs;
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
            polaritiesInt.push_back(static_cast<int>(meta.second.polarities[i]));
            isolationMins.push_back(fi.isolationRanges[j].start);
            isolationMaxs.push_back(fi.isolationRanges[j].end);
            precursorMZs.push_back(fi.precursorMZs[j]);
            subScans.push_back(fi.subScans[j]);
            if (!fi.subScanEnds.empty())
                subScanEnds.push_back(fi.subScanEnds[j]);
        }
    }
    
    auto ret = Rcpp::DataFrame::create(Rcpp::Named("scan") = scans,
                                       Rcpp::Named("time") = times,
                                       Rcpp::Named("TIC") = TICs,
                                       Rcpp::Named("BPC") = BPCs,
                                       Rcpp::Named("polarity") = polaritiesInt,
                                       Rcpp::Named("isolationRangeMin") = isolationMins,
                                       Rcpp::Named("isolationRangeMax") = isolationMaxs,
                                       Rcpp::Named("precursorMZ") = precursorMZs,
                                       Rcpp::Named("subScan") = subScans);
    if (!subScanEnds.empty())
        ret["subScanEnd"] = subScanEnds;
    
    return ret;
}

// [[Rcpp::export]]
void setSpecMetadata(MSReadBackend &backend, const Rcpp::DataFrame &mdMS, const Rcpp::DataFrame &mdMSMS)
{
    const auto polsFromInts = [](const auto &v)
    {
        std::vector<SpectrumRawTypes::MSPolarity> ret(v.size());
        for (size_t i=0; i<v.size(); ++i)
            ret[i] = static_cast<SpectrumRawTypes::MSPolarity>(v[i]);
        return ret;
    };
    
    SpectrumRawMetadata meta;

    // MS
    meta.first.scans = Rcpp::as<std::vector<SpectrumRawTypes::Scan>>(mdMS["scan"]);
    meta.first.times = Rcpp::as<std::vector<SpectrumRawTypes::Time>>(mdMS["time"]);
    meta.first.TICs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMS["TIC"]);
    meta.first.BPCs = Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(mdMS["BPC"]);
    meta.first.polarities = polsFromInts(Rcpp::as<std::vector<int>>(mdMS["polarity"]));
    
    // MSMS
    std::vector<SpectrumRawTypes::Scan> R_scans = mdMSMS["scan"];
    std::vector<SpectrumRawTypes::Time> R_times = mdMSMS["time"];
    std::vector<SpectrumRawTypes::Intensity> R_TICs = mdMSMS["TIC"], R_BPCs = mdMSMS["BPC"];
    auto R_polarities = polsFromInts(Rcpp::as<std::vector<int>>(mdMSMS["polarity"]));
    std::vector<SpectrumRawTypes::Mass> R_isoStarts = mdMSMS["isolationRangeMin"], R_isoEnds = mdMSMS["isolationRangeMax"];
    std::vector<SpectrumRawTypes::Mass> R_precursorMZs = mdMSMS["precursorMZ"];
    
    const std::vector<std::string> cn = mdMSMS.names();
    
    if (std::find(cn.begin(), cn.end(), "subScan") == cn.end()) // non-IMS
    {
        meta.second.scans = std::move(R_scans);
        meta.second.times = std::move(R_times);
        meta.second.TICs = std::move(R_TICs);
        meta.second.BPCs = std::move(R_BPCs);
        meta.second.polarities = std::move(R_polarities);
        
        for (size_t i=0; i<R_isoStarts.size(); ++i)
            meta.second.isolationRanges.emplace_back(R_isoStarts[i], R_isoEnds[i]);
        
        meta.second.precursorMZs = std::move(R_precursorMZs);
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
                meta.second.polarities.push_back(R_polarities[i]);
                if (!curFI.empty())
                {
                    meta.second.MSMSFrames.push_back(std::move(curFI));
                    curFI.clear();
                }
            }
            curFI.isolationRanges.emplace_back(R_isoStarts[i], R_isoEnds[i]);
            curFI.precursorMZs.push_back(R_precursorMZs[i]);
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
                          const std::vector<SpectrumRawTypes::Mass> &precursorMZs,
                          SpectrumRawTypes::Mass fixedIsolationWidth, bool withPrecursor,
                          bool retainPrecursor, int MSLevel, const std::string &method, SpectrumRawTypes::Mass mzWindow,
                          const std::vector<SpectrumRawTypes::Mobility> startMobs,
                          const std::vector<SpectrumRawTypes::Mobility> endMobs,
                          SpectrumRawTypes::PeakAbundance minAbundanceRel,
                          SpectrumRawTypes::PeakAbundance minAbundanceAbs,
                          SpectrumRawTypes::PeakAbundance minAbundanceIMSRel,
                          SpectrumRawTypes::PeakAbundance minAbundanceIMSAbs, unsigned topMost,
                          SpectrumRawTypes::Intensity minIntensityIMS, SpectrumRawTypes::Intensity minIntensityPre,
                          SpectrumRawTypes::Intensity minIntensityPost, SpectrumRawTypes::Intensity minBPIntensity)
{
    const auto entries = startTimes.size();
    const auto clMethod = clustMethodFromStr(method);
    const auto MSLev = (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2;
    const auto specMeta = backend.getSpecMetadata();
    const auto baseSpecFilter = SpectrumRawFilter()
        .setTopMost(topMost)
        .setWithPrecursor(withPrecursor)
        .setRetainPrecursor(retainPrecursor);
    const auto specFilter = SpectrumRawFilter(baseSpecFilter).setMinIntensity(minIntensityPre);
    const auto specFilterIMS = SpectrumRawFilter(baseSpecFilter);
    
    // NOTE: for IMS data, averageSpectraRaw() is called which returns a SpectrumRawAveraged. Since we don't care about
    // the additional metadata from this class, we purposely slice it by explicitly specifying the lambda's return type.
    const auto &sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t e) -> SpectrumRaw
    {
        const bool hasMob = spec.hasMobilities();
        
        if (hasMob)
        {
            const auto specf = filterIMSFrame(spec, specFilterIMS, precursorMZs[e],
                                              makeNumRange(startMobs[e], endMobs[e]));
            return averageSpectraRaw(specf, frameSubSpecIDs(specf), clMethod, mzWindow, false, minIntensityPre,
                                     minAbundanceIMSRel, minAbundanceIMSAbs);
        }
        return filterSpectrumRaw(spec, specFilter, precursorMZs[e]);
    };
    
    std::vector<std::vector<SpectrumRawSelection>> scanSels;
    for (size_t i=0; i<entries; ++i)
    {
        scanSels.push_back(getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]), MSLev,
                                                precursorMZs[i], fixedIsolationWidth, minBPIntensity));
        /*Rcpp::Rcout << "ss " << i << "/" << precursorMZs[i] << ": ";
        for (const auto &ss : scanSels.back())
            Rcpp::Rcout << ss.index << "/" << specMeta.second.scans[ss.index] << " ";
        Rcpp::Rcout << "\n";*/
    }
    
    const auto allSpectra = applyMSData<SpectrumRaw>(backend, MSLev, scanSels, sfunc, minIntensityIMS,
                                                     SpectrumRawTypes::MSSortType::MOBILITY_MZ);
    
    std::vector<SpectrumRawAveraged> averagedSpectra(entries);
    #pragma omp parallel for
    for (size_t i=0; i<entries; ++i)
        averagedSpectra[i] = averageSpectraRaw(allSpectra[i], clMethod, mzWindow, true, minIntensityPost,
                                               minAbundanceRel, minAbundanceAbs);

    Rcpp::List ret(entries);
    for (size_t i=0; i<entries; ++i)
    {
        ret[i] = Rcpp::List::create(Rcpp::Named("mz") = averagedSpectra[i].getMZs(),
                                    Rcpp::Named("intensity") = averagedSpectra[i].getIntensities(),
                                    Rcpp::Named("abundance_rel") = averagedSpectra[i].getAbundancesRel(),
                                    Rcpp::Named("abundance_abs") = averagedSpectra[i].getAbundancesAbs());
    }

    return ret;
}

// [[Rcpp::export]]
Rcpp::List getEIMList(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Mass> &startMZs,
                      const std::vector<SpectrumRawTypes::Mass> &endMZs,
                      const std::vector<SpectrumRawTypes::Time> &startTimes,
                      const std::vector<SpectrumRawTypes::Time> &endTimes,
                      const std::vector<SpectrumRawTypes::Mobility> &startMobs,
                      const std::vector<SpectrumRawTypes::Mobility> &endMobs,
                      SpectrumRawTypes::Intensity minIntensity, SpectrumRawTypes::Mass mzExpIMSWindow, bool compress)
{
    const auto entries = startTimes.size();
    const auto specMeta = backend.getSpecMetadata();
    
    // tolerance used for summing up equal mobility peaks. Normally these values are equivalent between frames
    const SpectrumRawTypes::Mobility avgTol = 0.00001;
    
    struct EIM
    {
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        // UNDONE: also collect m/z?
        EIM(void) = default;
        EIM(size_t s) : mobilities(s), intensities(s) { }
        void append(SpectrumRawTypes::Mobility m, SpectrumRawTypes::Intensity i) { mobilities.push_back(m); intensities.push_back(i); }
    };
    
    const auto &sfunc = [mzExpIMSWindow, &startMZs, &endMZs, &startMobs, &endMobs](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t e)
    {
        if (!spec.empty() && !spec.hasMobilities())
            Rcpp::stop("Cannot load mobilogram: no mobility data found!");
        
        EIM ret;
        auto it = (startMobs[e] == 0.0) ?  spec.getMobilities().begin() :
            std::lower_bound(spec.getMobilities().begin(), spec.getMobilities().end(), startMobs[e]);
        while (it != spec.getMobilities().end() && (endMobs[e] == 0.0 || numberLTE(*it, endMobs[e])))
        {
            const auto nextMobIt = std::upper_bound(it, spec.getMobilities().end(), *it);
            const auto startInd = std::distance(spec.getMobilities().begin(), it);
            const auto endInd = std::distance(spec.getMobilities().begin(), nextMobIt);
            const auto endMzIt = std::next(spec.getMZs().begin(), endInd);
            auto mzIt = std::lower_bound(std::next(spec.getMZs().begin(), startInd), endMzIt,
                                         startMZs[e] - mzExpIMSWindow);
            SpectrumRawTypes::Intensity curIntensity = 0;
            for (; mzIt != endMzIt && numberLTE(*mzIt, endMZs[e] + mzExpIMSWindow); ++mzIt)
            {
                curIntensity += spec.getIntensities()[std::distance(spec.getMZs().begin(), mzIt)];
            }
            
            ret.mobilities.push_back(*it);
            ret.intensities.push_back(curIntensity);
            
            it = nextMobIt; // continue with next mobility
        }
        return ret;
    };
    
    std::vector<std::vector<SpectrumRawSelection>> scanSels;
    for (size_t i=0; i<entries; ++i)
    {
        scanSels.push_back(getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                                SpectrumRawTypes::MSLevel::MS1, 0));
    }
    
    const auto allEIMs = applyMSData<EIM>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels, sfunc, minIntensity,
                                          SpectrumRawTypes::MSSortType::MOBILITY_MZ);

    std::vector<EIM> averageEIMs(entries);

    #pragma omp parallel for
    for (size_t i=0; i<entries; ++i)
    {
        EIM flatEIM;
        for (const EIM &eim : allEIMs[i])
        {
            std::copy(eim.mobilities.begin(), eim.mobilities.end(), std::back_inserter(flatEIM.mobilities));
            std::copy(eim.intensities.begin(), eim.intensities.end(), std::back_inserter(flatEIM.intensities));
        }
        
        if (flatEIM.mobilities.size() == 0)
            continue;
        
        const auto clusts = clusterNums(flatEIM.mobilities, clusterMethod::DISTANCE, avgTol);
        const int maxClust = *(std::max_element(clusts.begin(), clusts.end()));
        EIM avgEIM(maxClust + 1);
        std::vector<size_t> clSizes(maxClust + 1);
        
        // sum data, averaging is done below
        // NOTE: for now don't consider weighting averages, since there are likely to be clusters with zero intensities
        for (size_t j=0; j<clusts.size(); ++j)
        {
            const size_t cl = clusts[j];
            avgEIM.mobilities[cl] += flatEIM.mobilities[j];
            avgEIM.intensities[cl] += flatEIM.intensities[j];
            ++clSizes[cl];
        }
        
        const auto EIMSize = avgEIM.mobilities.size();
        
        // average data
        for (size_t j=0; j<EIMSize; ++j)
        {
            avgEIM.mobilities[j] /= static_cast<SpectrumRawTypes::Mobility>(clSizes[j]); // mean of all values
            avgEIM.intensities[j] /= allEIMs[i].size(); // mean of values (including frames without this cluster)
        }
        
        // sort and compress data
        const auto sInds = getSortedInds(avgEIM.mobilities);
        EIM sortedEIM;
        const bool doComp = compress && EIMSize >= 3; // we always keep first/last point, so need >=3 points
        for (size_t i=0; i<EIMSize; ++i)
        {
            // if we want to (1) compress, (2) current intensity == 0, (3) is not the first sorted point, (4) is not the
            // last point to check and (5) and not the last sorted point.
            if (doComp && avgEIM.intensities[sInds[i]] == 0 && sInds[i] > 0 && i < (EIMSize-1) &&
                sInds[i] != (EIMSize-1))
            {
                const auto prevInt = avgEIM.intensities[sInds[i-1]], nextInt = avgEIM.intensities[sInds[i+1]];
                if (prevInt == 0 && nextInt == 0)
                    continue; // skip points with zero intensities that are neighbored by others.
            }
            sortedEIM.append(avgEIM.mobilities[sInds[i]], avgEIM.intensities[sInds[i]]);
        }
        
        averageEIMs[i] = std::move(sortedEIM);
    }
        
    Rcpp::List ret(entries);
    for (size_t i=0; i<entries; ++i)
    {
        auto mat = Rcpp::NumericMatrix(averageEIMs[i].mobilities.size(), 2);
        for (size_t j=0; j<averageEIMs[i].mobilities.size(); ++j)
        {
            mat(j, 0) = averageEIMs[i].mobilities[j];
            mat(j, 1) = averageEIMs[i].intensities[j];
        }
        Rcpp::colnames(mat) = Rcpp::CharacterVector::create("mobility", "intensity");
        ret[i] = mat;
    }
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix compressEIM(const std::vector<SpectrumRawTypes::Mobility> &mobilities,
                                const std::vector<SpectrumRawTypes::Intensity> &intensities)
{
    // compress EIM data by removing zero intensity points that are neighbored by other zero intensity points.
    // NOTE: we always keep the first and last point, so we need at least 3 points to compress.
    if (mobilities.size() < 3)
    {
        auto mat = Rcpp::NumericMatrix(mobilities.size(), 2);
        for (size_t i=0; i<mobilities.size(); ++i)
        {
            mat(i, 0) = mobilities[i];
            mat(i, 1) = intensities[i];
        }
        Rcpp::colnames(mat) = Rcpp::CharacterVector::create("mobility", "intensity");
        return mat;
    }
    
    std::vector<SpectrumRawTypes::Mobility> outMobs;
    std::vector<SpectrumRawTypes::Intensity> outInts;
    
    for (size_t i=0; i<mobilities.size(); ++i)
    {
        if (intensities[i] == 0 && i > 0 && i < (mobilities.size() - 1) &&
            intensities[i - 1] == 0 && intensities[i + 1] == 0)
            continue; // skip zero intensity point that is neighbored by other zero intensity points
        
        outMobs.push_back(mobilities[i]);
        outInts.push_back(intensities[i]);
    }
    
    auto mat = Rcpp::NumericMatrix(outMobs.size(), 2);
    for (size_t i=0; i<outMobs.size(); ++i)
    {
        mat(i, 0) = outMobs[i];
        mat(i, 1) = outInts[i];
    }
    Rcpp::colnames(mat) = Rcpp::CharacterVector::create("mobility", "intensity");
    return mat;
}

// [[Rcpp::export]]
Rcpp::NumericVector getPeakIntensities(const MSReadBackend &backend,
                                       const std::vector<SpectrumRawTypes::Mass> &startMZs,
                                       const std::vector<SpectrumRawTypes::Mass> &endMZs,
                                       const std::vector<SpectrumRawTypes::Time> &times)
{
    const auto &specMeta = backend.getSpecMetadata().first;
    const auto entries = startMZs.size();
    std::vector<std::vector<SpectrumRawSelection>> scanSels(entries);
    
    if (entries == 0)
        return Rcpp::NumericVector();

    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t e)
    {
        const auto startIt = std::lower_bound(spec.getMZs().begin(), spec.getMZs().end(), startMZs[e]);
        if (startIt == spec.getMZs().end())
            return 0;
        const auto endIt = std::upper_bound(startIt, spec.getMZs().end(), endMZs[e]);
        const auto startInd = std::distance(spec.getMZs().begin(), startIt);
        const auto endInd = std::distance(spec.getMZs().begin(), endIt);
        return std::accumulate(spec.getIntensities().begin() + startInd, spec.getIntensities().begin() + endInd, 0);
    };
        
    for (size_t i=0; i<entries; ++i)
    {
        // use lower bound to quickly find the index that matches the given scan time. However, since given RTs are
        // possibly not exact and will not exactly match the actual scan times, it may be that the lower_bound() returns
        // the first scan that is higher, while the previous may actually be closer. Hence, also check the previous scan
        // and take the closest.
        const auto it = std::lower_bound(specMeta.times.begin(), specMeta.times.end(), times[i]);
        auto ind = std::distance(specMeta.times.begin(), it);
        if (it != specMeta.times.begin()) // only check previous if there is one
        {
            const auto prevIt = std::prev(it);
            if (it == specMeta.times.end() || (std::fabs(times[i] - *it) > std::fabs(times[i] - *prevIt)))
                --ind; // use the previous if the RT was not matched or is less close
        }
        scanSels[i].emplace_back(ind);
    }
    
    const auto ints = applyMSData<SpectrumRawTypes::Intensity>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels,
                                                               sfunc, 0, SpectrumRawTypes::MSSortType::MZ);
    
    auto ret = Rcpp::NumericVector(entries);
    for (size_t i=0; i<entries; ++i)
        ret[i] = (!ints[i].empty()) ? ints[i][0] : 0;
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::List collapseIMSFrames(const MSReadBackend &backend, SpectrumRawTypes::Mass mzStart, SpectrumRawTypes::Mass mzEnd,
                             SpectrumRawTypes::Mobility mobilityStart, SpectrumRawTypes::Mobility mobilityEnd,
                             const std::string &method, SpectrumRawTypes::Mass mzWindow,
                             SpectrumRawTypes::PeakAbundance minAbundanceRel,
                             SpectrumRawTypes::PeakAbundance minAbundanceAbs, unsigned topMost,
                             SpectrumRawTypes::Intensity minIntensityIMS, SpectrumRawTypes::Intensity minIntensityPre,
                             bool includeMSMS)
{
    const auto clMethod = clustMethodFromStr(method);
    const auto filterP = SpectrumRawFilter()
        .setMZRange(mzStart, mzEnd)
        .setTopMost(topMost);
    const auto mobRange = makeNumRange(mobilityStart, mobilityEnd);
    
    const auto &sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        if (!spec.hasMobilities() && !spec.empty()) // NOTE: we cannot distinguish between non-IMS and empty spectra at the moment
            Rcpp::stop("Tried to collapse non-IMS data!");
        
        const auto specf = filterIMSFrame(spec, filterP, 0.0, mobRange);
        return averageSpectraRaw(specf, frameSubSpecIDs(specf), clMethod, mzWindow, false, minIntensityPre,
                                 minAbundanceRel, minAbundanceAbs);
    };
    
    const auto &specMetaMS = backend.getSpecMetadata().first;
    std::vector<std::vector<SpectrumRawSelection>> scanSelsMS(1);
    for (size_t i=0; i<specMetaMS.scans.size(); ++i)
        scanSelsMS[0].emplace_back(i);
    
    const auto spectraMS = applyMSData<SpectrumRawAveraged>(backend, SpectrumRawTypes::MSLevel::MS1, scanSelsMS, sfunc,
                                                            minIntensityIMS,
                                                            SpectrumRawTypes::MSSortType::MOBILITY_MZ)[0];
    
    const auto &getSpecRList = [](const auto &spectra)
    {
        // NOTE: we return matrices so these can be directly consumed by mzR
        Rcpp::List ret(spectra.size());
        const auto coln = Rcpp::CharacterVector::create("mz", "intensity");
        for (size_t i=0; i<spectra.size(); ++i)
        {
            Rcpp::NumericMatrix m(spectra[i].size(), 2);
            Rcpp::NumericVector mzs = Rcpp::wrap(spectra[i].getMZs()), ints = Rcpp::wrap(spectra[i].getIntensities());
            m(Rcpp::_, 0) = mzs; m(Rcpp::_, 1) = ints;
            Rcpp::colnames(m) = coln;
            ret[i] = m;
        }
        return ret;
    };
    
    if (!includeMSMS)
        return Rcpp::List::create(Rcpp::Named("MS1") = getSpecRList(spectraMS));

        
    const auto &specMetaMS2 = backend.getSpecMetadata().second;
    std::vector<std::vector<SpectrumRawSelection>> scanSelsMS2(1);
    std::vector<SpectrumRawTypes::Scan> framesMS2;
    std::vector<SpectrumRawTypes::Mass> scanPrecursorMZs, isolationStarts, isolationEnds;
    for (size_t i=0; i<specMetaMS2.scans.size(); ++i)
    {
        // special case, eg DIA
        if (specMetaMS2.MSMSFrames[i].isolationRanges.size() < 2)
        {
            scanSelsMS2[0].emplace_back(i);
            framesMS2.push_back(specMetaMS2.scans[i]);
            if (specMetaMS2.MSMSFrames[i].isolationRanges.empty())
            {
                scanPrecursorMZs.push_back(0.0); isolationStarts.push_back(0.0); isolationEnds.push_back(0.0);
            }
            else
            {
                const auto ir = specMetaMS2.MSMSFrames[i].isolationRanges[0];
                isolationStarts.push_back(ir.start); isolationEnds.push_back(ir.end);
                scanPrecursorMZs.push_back(specMetaMS2.MSMSFrames[i].precursorMZs[0]);
            }
            continue;
        }
        
        // for PASEF data we combine spectra with close precursor m/z
        const auto precMZClusts = clusterNums(specMetaMS2.MSMSFrames[i].precursorMZs, clMethod, mzWindow);
        const int maxClust = *(std::max_element(precMZClusts.begin(), precMZClusts.end()));
        
        for (int cl=0; cl<=maxClust; ++cl)
        {
            SpectrumRawSelection ssel(i);
            SpectrumRawTypes::Mass prec = 0.0;
            SpectrumRawTypes::IsolationRange ir(0.0, 0.0);
            for (size_t j=0; j<specMetaMS2.MSMSFrames[i].precursorMZs.size(); ++j)
            {
                if (cl == precMZClusts[j])
                {
                    ssel.MSMSFrameIndices.push_back(j);
                    prec += specMetaMS2.MSMSFrames[i].precursorMZs[j];
                    ir.start += specMetaMS2.MSMSFrames[i].isolationRanges[j].start;
                    ir.end += specMetaMS2.MSMSFrames[i].isolationRanges[j].end;
                }
            }
            scanSelsMS2[0].push_back(ssel);
            framesMS2.push_back(specMetaMS2.scans[i]);
            const auto size = static_cast<SpectrumRawTypes::Mass>(ssel.MSMSFrameIndices.size());
            scanPrecursorMZs.push_back(prec / size);
            isolationStarts.push_back(ir.start / size); isolationEnds.push_back(ir.end / size);
        }
    }
    
    const auto spectraMS2 = applyMSData<SpectrumRawAveraged>(backend, SpectrumRawTypes::MSLevel::MS2, scanSelsMS2, sfunc,
                                                             minIntensityIMS,
                                                             SpectrumRawTypes::MSSortType::MOBILITY_MZ)[0];
    
    return Rcpp::List::create(Rcpp::Named("MS1") = getSpecRList(spectraMS),
                              Rcpp::Named("MS2") = getSpecRList(spectraMS2),
                              Rcpp::Named("framesMS2") = framesMS2,
                              Rcpp::Named("precursorMZs") = scanPrecursorMZs,
                              Rcpp::Named("isolationStarts") = isolationStarts,
                              Rcpp::Named("isolationEnds") = isolationEnds);
}

// [[Rcpp::export]]
Rcpp::List getIMSIsolationInfo(const MSReadBackend &backend)
{
    const auto &specMeta = backend.getSpecMetadata().second;
    
    using MobRangeAndTIC = std::pair<SpectrumRawTypes::MobilityRange, SpectrumRawTypes::Intensity>;
    
    const auto &sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t)
    {
        if (!spec.empty() && !spec.hasMobilities())
            throw(std::runtime_error("Need IMS data!"));
        
        SpectrumRawTypes::Mobility minMob = -1.0, maxMob = 0.0;
        SpectrumRawTypes::Intensity totInt = 0.0;
        for (size_t i=0; i<spec.size(); ++i)
        {
            const auto m = spec.getMobilities()[i];
            if (minMob == -1.0)
                minMob = maxMob = m;
            else
            {
                minMob = std::min(minMob, m);
                maxMob = std::max(maxMob, m);
            }
            totInt += spec.getIntensities()[i];
        }
        
        return MobRangeAndTIC(SpectrumRawTypes::MobilityRange(minMob, maxMob), totInt);
    };
    
    std::vector<std::vector<SpectrumRawSelection>> scanSels(1);
    std::vector<SpectrumRawTypes::Mass> isolationMins, isolationMaxs, precursorMZs;
    std::vector<SpectrumRawTypes::Time> times;
    for (size_t i=0; i<specMeta.scans.size(); ++i)
    {
        for (size_t j=0; j<specMeta.MSMSFrames[i].precursorMZs.size(); ++j)
        {
            SpectrumRawSelection ssel(i);
            ssel.MSMSFrameIndices.push_back(j);
            scanSels[0].push_back(std::move(ssel));
            times.push_back(specMeta.times[i]);
            isolationMins.push_back(specMeta.MSMSFrames[i].isolationRanges[j].start);
            isolationMaxs.push_back(specMeta.MSMSFrames[i].isolationRanges[j].end);
            precursorMZs.push_back(specMeta.MSMSFrames[i].precursorMZs[j]);
        }
    }
    
    if (isolationMins.empty())
        return Rcpp::List();
    
    const auto mats = applyMSData<MobRangeAndTIC>(backend, SpectrumRawTypes::MSLevel::MS2, scanSels, sfunc, 0,
                                                  SpectrumRawTypes::MSSortType::NONE)[0];
    std::vector<SpectrumRawTypes::Mobility> mobStarts(mats.size());
    std::vector<SpectrumRawTypes::Mobility> mobEnds(mats.size());
    std::vector<SpectrumRawTypes::Intensity> TICs(mats.size());
    for (size_t i=0; i<mats.size(); ++i)
    {
        mobStarts[i] = mats[i].first.start;
        mobEnds[i] = mats[i].first.end;
        TICs[i] = mats[i].second;
    }

    return Rcpp::List::create(Rcpp::Named("time") = times,
                              Rcpp::Named("isolationRangeMin") = isolationMins,
                              Rcpp::Named("isolationRangeMax") = isolationMaxs,
                              Rcpp::Named("precursorMZ") = precursorMZs,
                              Rcpp::Named("mobStart") = mobStarts,
                              Rcpp::Named("mobEnd") = mobEnds,
                              Rcpp::Named("TIC") = TICs);
}

// [[Rcpp::export]]
void testMS1Writer(const MSReadBackend &backend, const std::string &out, SpectrumRawTypes::Mass mzStart,
                   SpectrumRawTypes::Mass mzEnd, SpectrumRawTypes::Mobility mobilityStart,
                   SpectrumRawTypes::Mobility mobilityEnd, const std::string &method, SpectrumRawTypes::Mass mzWindow,
                   SpectrumRawTypes::PeakAbundance minAbundance, unsigned topMost,
                   SpectrumRawTypes::Intensity minIntensityIMS, SpectrumRawTypes::Intensity minIntensityPre)
{
#if 0 // MSTK writer was removed
    
#ifdef WITH_MSTK
    const auto clMethod = clustMethodFromStr(method);
    const auto filterP = SpectrumRawFilter()
        .setMinIntensity(minIntensityIMS)
        .setMZRange(mzStart, mzEnd)
        .setTopMost(topMost);
    const auto mobRange = makeNumRange(mobilityStart, mobilityEnd);
    
    const auto &sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t) -> SpectrumRaw
    {
        if (!spec.hasMobilities())
            Rcpp::stop("Tried to collapse non-IMS data!");
        
        const auto specf = filterIMSFrame(spec, filterP, 0.0, mobRange);
        return averageSpectraRaw(specf, frameSubSpecIDs(specf), clMethod, mzWindow, false, minIntensityPre,
                                 minAbundance);
    };
    
    const auto &specMeta = backend.getSpecMetadata().first;
    std::vector<std::vector<SpectrumRawSelection>> scanSels(1);
    for (size_t i=0; i<specMeta.scans.size(); ++i)
        scanSels[0].emplace_back(i);
    
    const auto spectra = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels, sfunc, 0)[0];
    
    writeMS1SpectraMSTK(out, spectra, specMeta);
#endif
    
#endif
}

// [[Rcpp::export]]
Rcpp::List getChromPoints(const MSReadBackend &backend, SpectrumRawTypes::Time rtStart, SpectrumRawTypes::Time rtEnd,
                          SpectrumRawTypes::Mass mzStart, SpectrumRawTypes::Mass mzEnd, bool withMob,
                          SpectrumRawTypes::Mobility mobilityStart, SpectrumRawTypes::Mobility mobilityEnd)
{
    const auto specf = SpectrumRawFilter().setMZRange(mzStart, mzEnd);
    const auto mobRange = makeNumRange(mobilityStart, mobilityEnd);
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        return (spec.hasMobilities()) ? filterIMSFrame(spec, specf, 0.0, mobRange) : filterSpectrumRaw(spec, specf, 0.0);
    };
    
    const auto &meta = backend.getSpecMetadata();
    
    std::vector<std::vector<SpectrumRawSelection>> sels(1);
    auto startIt = std::lower_bound(meta.first.times.begin(), meta.first.times.end(), rtStart);
    std::vector<SpectrumRawTypes::Time> times;
    for (; startIt != meta.first.times.end() && numberLTE(*startIt, rtEnd); ++startIt)
    {
        const auto ind = std::distance(meta.first.times.begin(), startIt);
        sels[0].emplace_back(ind);
        times.push_back(*startIt);
    }

    const auto specs = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc, 25,
                                                SpectrumRawTypes::MSSortType::MOBILITY_MZ)[0];
    SpectrumRaw flatSpec;
    std::vector<SpectrumRawTypes::Time> flatTimes;
    for (size_t i=0; i<specs.size(); ++i)
    {
        flatSpec.append(specs[i]);
        flatTimes.insert(flatTimes.end(), specs[i].size(), times[i]);
    }
    
    if (!withMob)
    {
        return Rcpp::List::create(Rcpp::Named("time") = flatTimes,
                                  Rcpp::Named("mz") = flatSpec.getMZs(),
                                  Rcpp::Named("intensity") = flatSpec.getIntensities());
    }
    
    return Rcpp::List::create(Rcpp::Named("time") = flatTimes,
                              Rcpp::Named("mz") = flatSpec.getMZs(),
                              Rcpp::Named("mobility") = flatSpec.getMobilities(),
                              Rcpp::Named("intensity") = flatSpec.getIntensities());
}
