/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <Rcpp.h>

#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "msdata.h"
#include "msdata.hpp"
#include "msdata-mem.h"
#include "msdata-mstk.h"
#include "msdata-otims.h"
#include "msdata-sc.h"
#include "spectrum-raw.h"


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

std::vector<SpectrumRawTypes::Mobility> MSReadBackend::generateMobilities()
{
    if (!mobilities.empty())
        return mobilities; // may already have been loaded or pre-loaded by backend
    
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        return spec.getMobilities();
    };
    std::vector<std::vector<SpectrumRawSelection>> sels(1);
    const auto &meta = getSpecMetadata();
    for (size_t i=0; i<meta.first.scans.size(); ++i)
        sels[0].emplace_back(i);
    const auto ret = applyMSData<std::vector<SpectrumRawTypes::Mobility>>(*this, SpectrumRawTypes::MSLevel::MS1, sels,
                                                                          sfunc, 0, SpectrumRawTypes::MSSortType::NONE);
    
    // sort and make unique
    std::set<SpectrumRawTypes::Mobility> unMobs;
    for (const auto &mvec : ret[0])
        unMobs.insert(mvec.begin(), mvec.end());
    mobilities.assign(unMobs.begin(), unMobs.end());
    
    return mobilities;
}

RCPP_MODULE(MSReadBackend)
{
    Rcpp::class_<MSReadBackend>("MSReadBackend")
        .method("setNeedIMS", &MSReadBackend::setNeedIMS)
        .method("getNeedIMS", &MSReadBackend::getNeedIMS)
        .method("open", &MSReadBackend::open)
        .method("close", &MSReadBackend::close)
        .method("getCurrentFile", &MSReadBackend::getCurrentFile)
        .method("generateMobilities", &MSReadBackend::generateMobilities)
        // NOTE: select right overload
        .method("setMobilities", static_cast<void (MSReadBackend::*)(const std::vector<SpectrumRawTypes::Mobility> &)>(&MSReadBackend::setMobilities))
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
Rcpp::List getStats(const MSReadBackend &backend)
{
    const auto sfunc = [](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t)
    {
        size_t specs = 0, peaks = 0;
        SpectrumRawTypes::Mobility lastMob = -1;
        for (size_t i=0; i<spec.getMobilities().size(); ++i)
        {
            if (lastMob != spec.getMobilities()[i])
            {
                lastMob = spec.getMobilities()[i];
                ++specs;
            }
            ++peaks;
        }
        return std::make_pair(specs, peaks);
    };
    const auto &meta = backend.getSpecMetadata();
    
    std::vector<std::vector<SpectrumRawSelection>> sels(1);
    for (size_t i=0; i<meta.first.scans.size(); ++i)
        sels[0].emplace_back(i);
    const auto statsMS1 = applyMSData<std::pair<size_t, size_t>>(backend, SpectrumRawTypes::MSLevel::MS1, sels, sfunc, 0,
                                                                 SpectrumRawTypes::MSSortType::NONE)[0];
    const auto totMS1 = std::accumulate(statsMS1.begin(), statsMS1.end(), std::make_pair(0, 0),
                                        [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b)
                                      {
                                          return std::make_pair(a.first + b.first, a.second + b.second);
                                      });

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
    const auto statsMS2 = applyMSData<std::pair<size_t, size_t>>(backend, SpectrumRawTypes::MSLevel::MS2, sels, sfunc, 0,
                                              SpectrumRawTypes::MSSortType::NONE)[0];
    const auto totMS2 = std::accumulate(statsMS2.begin(), statsMS2.end(), std::make_pair(0, 0),
                                        [](const std::pair<size_t, size_t> &a, const std::pair<size_t, size_t> &b)
                                        {
                                            return std::make_pair(a.first + b.first, a.second + b.second);
                                        });
    
    
    return Rcpp::List::create(Rcpp::Named("specsMS1") = totMS1.first,
                              Rcpp::Named("pointsMS1") = totMS1.second,
                              Rcpp::Named("specsMS2") = totMS2.first,
                              Rcpp::Named("pointsMS2") = totMS2.second);
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
Rcpp::List getMSSpectra(const MSReadBackend &backend, SpectrumRawTypes::Time timeStart, SpectrumRawTypes::Time timeEnd,
                        SpectrumRawTypes::Mobility mobStart, SpectrumRawTypes::Mobility mobEnd,
                        int MSLevel, SpectrumRawTypes::Mass prec = 0.0, SpectrumRawTypes::Intensity minIntensityIMS = 25,
                        SpectrumRawTypes::Intensity minIntensityPre = 0, int topMost = 50, bool withPrecursor = false,
                        int minAbundanceIMSAbs = 2, bool retainPrecursor = true)
{
    
    const auto baseSpecFilter = SpectrumRawFilter()
        .setTopMost(topMost)
        .setWithPrecursor(withPrecursor)
        .setRetainPrecursor(retainPrecursor);
    const auto specFilter = SpectrumRawFilter(baseSpecFilter).setMinIntensity(minIntensityPre);
    const auto specFilterIMS = SpectrumRawFilter(baseSpecFilter);
    
    // NOTE: for IMS data, averageSpectraRaw() is called which returns a SpectrumRawAveraged. Since we don't care about
    // the additional metadata from this class, we purposely slice it by explicitly specifying the lambda's return type.
    const auto &sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t) -> SpectrumRaw
    {
        const bool hasMob = spec.hasMobilities();
        
        if (hasMob)
        {
            const auto specf = filterIMSFrame(spec, specFilterIMS, prec, makeNumRange(mobStart, mobEnd));
            // return specf;
            return averageSpectraRaw(specf, frameSubSpecIDs(specf), clusterMethod::DISTANCE, 0.005, false,
                                     minIntensityPre, 0, minAbundanceIMSAbs);
        }
        return filterSpectrumRaw(spec, specFilter, prec);
    };
    
    const auto MSLev = (MSLevel == 1) ? SpectrumRawTypes::MSLevel::MS1 : SpectrumRawTypes::MSLevel::MS2;
    const auto sels = getSpecRawSelections(backend.getSpecMetadata(), makeNumRange(timeStart, timeEnd), MSLev,
                                           prec, 0.0, 0);
    const std::vector<std::vector<SpectrumRawSelection>> selsList(1, std::vector<SpectrumRawSelection>{sels});
    
    const auto spectra = applyMSData<SpectrumRaw>(backend, MSLev, selsList, sfunc, minIntensityIMS,
                                                  SpectrumRawTypes::MSSortType::MOBILITY_MZ)[0];
    
    Rcpp::List ret(spectra.size());
    for (size_t i=0; i<spectra.size(); ++i)
    {
        const auto &spec = spectra[i];
        if (!spec.getMobilities().empty())
            ret[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                             Rcpp::Named("intensity") = spec.getIntensities(),
                                             Rcpp::Named("mobility") = spec.getMobilities());
        else
            ret[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = spec.getMZs(),
                                             Rcpp::Named("intensity") = spec.getIntensities());
    }
    
    return ret;
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
    
    struct EIM
    {
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        // UNDONE: also collect m/z?
        EIM(void) = default;
        EIM(size_t s) : mobilities(s), intensities(s) { }
        void append(SpectrumRawTypes::Mobility m, SpectrumRawTypes::Intensity i) { mobilities.push_back(m); intensities.push_back(i); }
        size_t size(void) const { return mobilities.size(); }
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
            
            if (curIntensity > 0)
            {
                ret.mobilities.push_back(*it);
                ret.intensities.push_back(curIntensity);
            }
            
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

    std::vector<EIM> summedEIMs(entries);
    const auto &allMobs = backend.getMobilities();

    #pragma omp parallel for
    for (size_t i=0; i<entries; ++i)
    {
        std::map<SpectrumRawTypes::Mobility, SpectrumRawTypes::Intensity> sumEIMMap;
        for (const EIM &eim : allEIMs[i])
        {
            for (size_t j=0; j<eim.size(); ++j)
                sumEIMMap[eim.mobilities[j]] += eim.intensities[j];
        }
        
        // merge in complete mobilities set to get gap free data
        auto it = (startMobs[i] == 0.0) ? allMobs.begin() : std::lower_bound(allMobs.begin(), allMobs.end(), startMobs[i]);
        for (; it != allMobs.end() && (endMobs[i] == 0.0 || numberLTE(*it, endMobs[i])); ++it)
            sumEIMMap.insert(std::make_pair(*it, 0)); // insert with zero intensity if not present
        
        // convert from map
        EIM sumEIM;
        for (const auto &p : sumEIMMap)
            sumEIM.append(p.first, p.second);
        
        // compress data
        const auto EIMSize = sumEIM.size();
        if (compress && EIMSize >= 3) // we always keep first/last point, so need >=3 points
        {
            EIM comprEIM;
            for (size_t j=0; j<EIMSize; ++j)
            {
                // if current intensity == 0, is not the first point, and is not the last point to check.
                if (sumEIM.intensities[j] == 0 && j != 0 && j < (EIMSize-1))
                {
                    const auto prevInt = comprEIM.intensities.back(), nextInt = sumEIM.intensities[j+1];
                    if (prevInt == 0 && nextInt == 0)
                        continue; // skip points with zero intensities that are neighbored by others.
                }
                comprEIM.append(sumEIM.mobilities[j], sumEIM.intensities[j]);
            }
            summedEIMs[i] = std::move(comprEIM);
        }
        else
            summedEIMs[i] = std::move(sumEIM);
    }
        
    Rcpp::List ret(entries);
    for (size_t i=0; i<entries; ++i)
    {
        auto mat = Rcpp::NumericMatrix(summedEIMs[i].mobilities.size(), 2);
        for (size_t j=0; j<summedEIMs[i].mobilities.size(); ++j)
        {
            mat(j, 0) = summedEIMs[i].mobilities[j];
            mat(j, 1) = summedEIMs[i].intensities[j];
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
Rcpp::List getChromPoints(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Time> &startTimes,
                          const std::vector<SpectrumRawTypes::Time> &endTimes,
                          const std::vector<SpectrumRawTypes::Mass> &startMZs,
                          const std::vector<SpectrumRawTypes::Mass> &endMZs,
                          const std::vector<SpectrumRawTypes::Mobility> &startMobs,
                          const std::vector<SpectrumRawTypes::Mobility> &endMobs,
                          bool withMob)
{
    const auto entries = startTimes.size();
    if (entries == 0)
        return Rcpp::List();

    const auto specMeta = backend.getSpecMetadata();

    std::vector<std::vector<SpectrumRawSelection>> scanSels(entries);
    std::vector<std::vector<SpectrumRawTypes::Time>> allTimes(entries);
    for (size_t i=0; i<entries; ++i)
    {
        const auto sels = getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                               SpectrumRawTypes::MSLevel::MS1, 0);
        std::vector<SpectrumRawTypes::Time> times(sels.size());
        for (size_t j=0; j<sels.size(); ++j)
            times[j] = specMeta.first.times[sels[j].index];
        allTimes[i] = std::move(times);
        scanSels[i] = std::move(sels);
    }

    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &, size_t entry)
    {
        const auto specf = SpectrumRawFilter().setMZRange(startMZs[entry], endMZs[entry]);
        const auto mobRange = makeNumRange(startMobs[entry], endMobs[entry]);
        return (spec.hasMobilities()) ? filterIMSFrame(spec, specf, 0.0, mobRange) : filterSpectrumRaw(spec, specf, 0.0);
    };

    const auto allSpectra = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels, sfunc, 25,
                                                     SpectrumRawTypes::MSSortType::MOBILITY_MZ);
    
    Rcpp::List ret(entries);
    for (size_t i = 0; i < entries; ++i)
    {
        const auto &specs = allSpectra[i];
        const auto &times = allTimes[i];
        SpectrumRaw flatSpec;
        std::vector<SpectrumRawTypes::Time> flatTimes;
        for (size_t j=0; j<specs.size(); ++j)
        {
            flatSpec.append(specs[j]);
            flatTimes.insert(flatTimes.end(), specs[j].size(), times[j]);
        }

        const auto nPoints = flatTimes.size();
        if (!withMob)
        {
            auto mat = Rcpp::NumericMatrix(nPoints, 3);
            for (size_t j = 0; j<nPoints; ++j)
            {
                mat(j, 0) = flatTimes[j];
                mat(j, 1) = flatSpec.getMZs()[j];
                mat(j, 2) = flatSpec.getIntensities()[j];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "mz", "intensity");
            ret[i] = mat;
        }
        else
        {
            auto mat = Rcpp::NumericMatrix(nPoints, 4);
            for (size_t j=0; j<nPoints; ++j)
            {
                mat(j, 0) = flatTimes[j];
                mat(j, 1) = flatSpec.getMZs()[j];
                mat(j, 2) = flatSpec.getMobilities()[j];
                mat(j, 3) = flatSpec.getIntensities()[j];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "mz", "mobility", "intensity");
            ret[i] = mat;
        }
    }

    return ret;
}
