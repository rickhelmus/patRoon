#include <Rcpp.h>

#include "spectrum-raw.h"
#include "utils.h"

#include <algorithm>

namespace{

template<typename Spec> Spec flattenSpectra(const std::vector<Spec> &spectra)
{
    Spec ret;
    for (const auto &sp : spectra)
        ret.append(sp);
    return ret;
}

template<typename Spec> std::vector<size_t> flattenedSpecIDs(const std::vector<Spec> &spectra)
{
    std::vector<size_t> ret;
    size_t curID = 0;
    for (const auto &sp : spectra)
    {
        ret.insert(ret.end(), sp.size(), curID);
        ++curID;
    }
    return ret;
}

}

// UNDONE: remove this overload?
std::vector<SpectrumRawSelection> getSpecRawSelections(const SpectrumRawMetadata &specMeta,
                                                       const NumRange<SpectrumRawTypes::Time> &timeRange,
                                                       SpectrumRawTypes::MSLevel MSLevel,
                                                       SpectrumRawTypes::Mass precursor)
{
    const SpectrumRawMetadataMS &metaMS = (MSLevel == SpectrumRawTypes::MSLevel::MS1) ? specMeta.first : specMeta.second;
    std::vector<SpectrumRawSelection> ret;
    const bool isMSMS = MSLevel == SpectrumRawTypes::MSLevel::MS2;
    const bool isIMSMSMS = isMSMS && specMeta.second.isolationRanges.empty();
    
    const auto startIt = std::lower_bound(metaMS.times.begin(), metaMS.times.end(), timeRange.start);
    if (startIt != metaMS.times.end())
    {
        for (size_t i=std::distance(metaMS.times.begin(), startIt);
             i<metaMS.times.size() && metaMS.times[i]<=timeRange.end; ++i)
        {
            SpectrumRawSelection sel;
            if (isIMSMSMS)
            {
                for (size_t j=0; j<specMeta.second.MSMSFrames[i].isolationRanges.size(); ++j)
                {
                    if (specMeta.second.MSMSFrames[i].isolationRanges[j].within(precursor))
                        sel.MSMSFrameIndices.push_back(j);
                }
                if (sel.MSMSFrameIndices.empty())
                    continue; // no MS/MS data for this one
            }
            else if (isMSMS && !specMeta.second.isolationRanges[i].within(precursor))
                continue;
            sel.index = i;
            ret.push_back(std::move(sel));
        }
    }
    
    return ret;
}

SpectrumRaw filterSpectrumRaw(const SpectrumRaw &spectrum, const SpectrumRawFilter &filter,
                              SpectrumRawTypes::Mass precursor)
{
    if (spectrum.empty())
        return spectrum;
    
    const SpectrumRawTypes::Mass precTol = 0.01; // UNDONE: configurable tolerance?
    size_t closestPrecMZInd = spectrum.size();
    SpectrumRawTypes::Mass minPrecMZDiff = -1.0;
    if (precursor != 0.0)
    {
        for (auto it = std::lower_bound(spectrum.getMZs().begin(), spectrum.getMZs().end(), precursor - precTol);
             it!=spectrum.getMZs().end(); ++it)
        {
            const SpectrumRawTypes::Mass d = std::abs(precursor - *it);
            if (d <= precTol && (minPrecMZDiff == -1.0 || d < minPrecMZDiff))
            {
                closestPrecMZInd = std::distance(spectrum.getMZs().begin(), it);
                minPrecMZDiff = d;
            }
            else
                break; // data is sorted: all next mzs will be greater and therefore with higher deviation
        }
    }
    
    if (filter.withPrecursor && minPrecMZDiff == -1.0)
        return SpectrumRaw(); // no precursor present
    
    size_t startInd = 0, endInd = spectrum.size();
    if (filter.mzRange.start > 0)
    {
        const auto it = std::lower_bound(spectrum.getMZs().begin(), spectrum.getMZs().end(), filter.mzRange.start);
        startInd = std::distance(spectrum.getMZs().begin(), it);
    }
    if (filter.mzRange.end > 0)
    {
        const auto it = std::upper_bound(spectrum.getMZs().begin() + startInd, spectrum.getMZs().end(), filter.mzRange.end);
        endInd = std::distance(spectrum.getMZs().begin(), it);
    }
    
    SpectrumRaw ret;
    SpectrumRawTypes::Intensity minInt = filter.minIntensity;
    const bool doTopMost = filter.topMost != 0 && (endInd - startInd) > filter.topMost;
    const bool hasMob = spectrum.hasMobilities();
    
    if (doTopMost)
    {
        std::vector<SpectrumRawTypes::Intensity> topMostInts(filter.topMost);
        std::partial_sort_copy(spectrum.getIntensities().begin() + startInd, spectrum.getIntensities().begin() + endInd,
                               topMostInts.begin(), topMostInts.end(),
                               [](SpectrumRawTypes::Intensity i1, SpectrumRawTypes::Intensity i2) { return i1 > i2; });
        minInt = std::max(minInt, topMostInts.back());
    }
    
    bool addedPrec = false;
    for (size_t i=startInd; i<endInd; i++)
    {
        if (spectrum.getIntensities()[i] >= minInt)
        {
            if (hasMob)
                ret.append(spectrum.getMZs()[i], spectrum.getIntensities()[i], spectrum.getMobilities()[i]);
            else
                ret.append(spectrum.getMZs()[i], spectrum.getIntensities()[i]);
            if (i == closestPrecMZInd)
                addedPrec = true;
        }
    }
    
    if (!addedPrec && closestPrecMZInd != spectrum.size())
    {
        const auto it = std::upper_bound(ret.getMZs().begin(), ret.getMZs().end(), spectrum.getMZs()[closestPrecMZInd]);
        if (hasMob)
        {
            ret.insert(std::distance(ret.getMZs().begin(), it), spectrum.getMZs()[closestPrecMZInd],
                       spectrum.getIntensities()[closestPrecMZInd], spectrum.getMobilities()[closestPrecMZInd]);
        }
        else
        {
            ret.insert(std::distance(ret.getMZs().begin(), it), spectrum.getMZs()[closestPrecMZInd],
                       spectrum.getIntensities()[closestPrecMZInd]);
        }
    }
    
    return ret;
}

SpectrumRaw filterIMSFrame(const SpectrumRaw &frame, const SpectrumRawFilter &filter,
                           SpectrumRawTypes::Mass precursor, const SpectrumRawTypes::MobilityRange &mobRange)
{
    if (frame.empty())
        return frame;
    
    SpectrumRaw ret;
    
    // filter each sub-spectrum and combine spectra back to output frame
    SpectrumRawTypes::Mobility curMob;
    SpectrumRaw curSpec;
    for (size_t i=0; ; ++i)
    {
        if (i == 0 || i >= frame.size() || curMob != frame.getMobilities()[i])
        {
            if (i > 0)
            {
                ret.append(filterSpectrumRaw(curSpec, filter, precursor));
                curSpec.clear();
            }
            
            if (i >= frame.size())
                break; // done
            
            curMob = frame.getMobilities()[i];
        }
        
        if (mobRange.isSet() && !mobRange.within(curMob))
            continue;
        
        curSpec.append(frame.getMZs()[i], frame.getIntensities()[i], frame.getMobilities()[i]);
    }
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::DataFrame testSpecFilter(const std::vector<SpectrumRawTypes::Mass> &mzs,
                               const std::vector<SpectrumRawTypes::Intensity> &ints, double mzMin, double mzMax,
                               double minInt, unsigned topMost, double prec)
{
    const SpectrumRaw spec(mzs, ints);
    const SpectrumRawFilter fil = SpectrumRawFilter().setMZRange(mzMin, mzMax).setMinIntensity(minInt).setTopMost(topMost);
    const SpectrumRaw specF = filterSpectrumRaw(spec, fil, prec);
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = specF.getMZs(),
                                   Rcpp::Named("intensity") = specF.getIntensities());
}

SpectrumRawAveraged averageSpectraRaw(const SpectrumRawAveraged &flattenedSpecs, const std::vector<size_t> &specIDs,
                                      clusterMethod method, SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance)
{
    if (flattenedSpecs.empty()) // all spectra are empty
        return SpectrumRawAveraged();
    
    const bool alreadyAveraged = !flattenedSpecs.getAbundances().empty();
    const std::vector<int> clusts = clusterNums(flattenedSpecs.getMZs(), method, window);
    const int maxClust = *(std::max_element(clusts.begin(), clusts.end()));
    SpectrumRawAveraged binnedSpectrum(maxClust + 1, false, alreadyAveraged);
    std::vector<std::set<size_t>> binIDs(maxClust + 1);
    
    // Rcpp::Rcout << "flattenedSpecs: " << flattenedSpecs.size() << "\n";
    // Rcpp::Rcout << "maxClust: " << maxClust << "\n";
    
    // sum data for each cluster
    for (size_t i=0; i<clusts.size(); ++i)
    {
        const size_t cl = clusts[i];
        const SpectrumRawTypes::Intensity intenfl = flattenedSpecs.getIntensities()[i];
        const SpectrumRawTypes::Mass mz = binnedSpectrum.getMZs()[cl] + (flattenedSpecs.getMZs()[i] *
                                                                static_cast<SpectrumRawTypes::Mass>(intenfl));
        const SpectrumRawTypes::Intensity inten = binnedSpectrum.getIntensities()[cl] + intenfl;
        
        binnedSpectrum.setPeak(cl, mz, inten);
        
        if (alreadyAveraged)
        {
            const SpectrumRawTypes::PeakAbundance ab = binnedSpectrum.getAvgPrevAbundances()[cl] + flattenedSpecs.getAbundances()[i];
            binnedSpectrum.setAvgPrevAbundance(cl, ab);
        }
        
        binIDs[cl].insert(specIDs[i]);
    }
    
    // NOTE: assume that IDs are sorted and a regular range without gaps
    const auto numSpecs = specIDs.back() + 1;
    
    // average data
    for (size_t i=0; i<binnedSpectrum.size(); ++i)
    {
        const SpectrumRawTypes::Mass mz = binnedSpectrum.getMZs()[i] / static_cast<SpectrumRawTypes::Mass>(binnedSpectrum.getIntensities()[i]);
        SpectrumRawTypes::Intensity inten = binnedSpectrum.getIntensities()[i];
        if (averageIntensities)
            inten /= static_cast<SpectrumRawTypes::Intensity>(numSpecs);
        binnedSpectrum.setPeak(i, mz, inten);
        
        if (alreadyAveraged)
        {
            binnedSpectrum.setAvgPrevAbundance(i, binnedSpectrum.getAvgPrevAbundances()[i] /
                static_cast<SpectrumRawTypes::PeakAbundance>(numSpecs));
        }
    }
    
    // Rcpp::Rcout << "binnedSpectrum: " << binnedSpectrum.size() << "\n";
    
    // sort spectrum && pre-treat
    const auto sortedInds = getSortedInds(binnedSpectrum.getMZs());
    SpectrumRawAveraged sortedSpectrum;
    for (size_t i=0; i<sortedInds.size(); ++i)
    {
        const auto j = sortedInds[i];
        if (minIntensity > 0 && binnedSpectrum.getIntensities()[j] < minIntensity)
            continue;
        
        const auto abundance = static_cast<SpectrumRawTypes::PeakAbundance>(binIDs[j].size()) /
            static_cast<SpectrumRawTypes::PeakAbundance>(numSpecs);
        if (abundance < minAbundance)
            continue;
        
        if (alreadyAveraged)
            sortedSpectrum.append(binnedSpectrum.getMZs()[j], binnedSpectrum.getIntensities()[j], abundance,
                                  binnedSpectrum.getAvgPrevAbundances()[j]);
        else
            sortedSpectrum.append(binnedSpectrum.getMZs()[j], binnedSpectrum.getIntensities()[j], abundance);
    }
    
    // Rcpp::Rcout << "sortedSpectrum: " << sortedSpectrum.size() << "\n";
    
    return sortedSpectrum;
}

SpectrumRawAveraged averageSpectraRaw(const SpectrumRaw &flattenedSpecs, const std::vector<size_t> &specIDs,
                                      clusterMethod method, SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance)
{
    return averageSpectraRaw(SpectrumRawAveraged(flattenedSpecs.getMZs(), flattenedSpecs.getIntensities()), specIDs,
                             method, window, averageIntensities, minIntensity, minAbundance);
}

SpectrumRawAveraged averageSpectraRaw(const std::vector<SpectrumRaw> &spectra, clusterMethod method,
                                      SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance)
{
    return averageSpectraRaw(flattenSpectra(spectra), flattenedSpecIDs(spectra), method, window, averageIntensities,
                             minIntensity, minAbundance);
}

SpectrumRawAveraged averageSpectraRaw(const std::vector<SpectrumRawAveraged> &spectra, clusterMethod method,
                                      SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance)
{
    return averageSpectraRaw(flattenSpectra(spectra), flattenedSpecIDs(spectra), method, window, averageIntensities,
                             minIntensity, minAbundance);
}

// [[Rcpp::export]]
Rcpp::DataFrame doAverageSpectra(Rcpp::List specs, const std::string &method, SpectrumRawTypes::Mass window,
                                 SpectrumRawTypes::Intensity minIntensity, SpectrumRawTypes::PeakAbundance minAbundance)
{
    const auto clMethod = clustMethodFromStr(method);
    
    std::vector<SpectrumRawAveraged> spectra(specs.size());
    bool alreadyAveraged = false, checkAveraged = false;
    for (size_t i=0; i<spectra.size(); ++i)
    {
        const Rcpp::DataFrame s = Rcpp::as<Rcpp::DataFrame>(specs[i]);
        const std::vector<double> mzs = s["mz"];
        const std::vector<double> intensities = s["intensity"];
        
        if (!checkAveraged)
        {
            const std::vector<std::string> cn = s.names();
            alreadyAveraged = std::find(cn.cbegin(), cn.cend(), "abundance") != cn.cend();
            checkAveraged = true;
        }

        const std::vector<double> abundances = (alreadyAveraged) ? s["abundance"] : std::vector<double>();

        for (size_t j=0; j<mzs.size(); ++j)
        {
            if (alreadyAveraged)
                spectra[i].append(mzs[j], intensities[j], abundances[j]);
            else
                spectra[i].append(mzs[j], intensities[j]);
        }
    }
    
    const auto avgsp = averageSpectraRaw(spectra, clMethod, window, true, minIntensity, minAbundance);
    
    if (alreadyAveraged)
    {
        return Rcpp::DataFrame::create(Rcpp::Named("mz") = avgsp.getMZs(),
                                       Rcpp::Named("intensity") = avgsp.getIntensities(),
                                       Rcpp::Named("abundance") = avgsp.getAbundances(),
                                       Rcpp::Named("abundance_prev") = avgsp.getAvgPrevAbundances());
    }
    
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = avgsp.getMZs(),
                                   Rcpp::Named("intensity") = avgsp.getIntensities(),
                                   Rcpp::Named("abundance") = avgsp.getAbundances());
}

// [[Rcpp::export]]
Rcpp::List doAverageSpectraList(Rcpp::List specsList, const std::string &method,
                                SpectrumRawTypes::Mass window, SpectrumRawTypes::Intensity minIntensity,
                                SpectrumRawTypes::PeakAbundance minAbundance)
{
    const size_t entries = specsList.size();
    const auto clMethod = clustMethodFromStr(method);

    std::vector<std::vector<SpectrumRawAveraged>> allSpectra;
    bool alreadyAveraged = false, checkedAveraged = false;
    for (size_t i=0; i<entries; ++i)
    {
        const Rcpp::List slist = specsList[i];
        const size_t ssize = slist.size();
        allSpectra.emplace_back(ssize);
        auto &spectra = allSpectra.back();
        
        for (size_t j=0; j<ssize; ++j)
        {
            const Rcpp::DataFrame s = Rcpp::as<Rcpp::DataFrame>(slist[j]);
            const std::vector<double> mzs = s["mz"];
            const std::vector<double> intensities = s["intensity"];
            
            if (!checkedAveraged)
            {
                const std::vector<std::string> cn = s.names();
                alreadyAveraged = std::find(cn.cbegin(), cn.cend(), "abundance") != cn.cend();
                checkedAveraged = true;
            }
            
            const std::vector<double> abundances = (alreadyAveraged) ? s["abundance"] : std::vector<double>();
            
            for (size_t k=0; k<mzs.size(); ++k)
            {
                if (alreadyAveraged)
                    spectra[j].append(mzs[k], intensities[k], abundances[k]);
                else
                    spectra[j].append(mzs[k], intensities[k]);
            }
        }
    }
    
    std::vector<SpectrumRawAveraged> averagedSpecList(entries);
    // UNDONE: num threads
    #pragma omp parallel for num_threads(12)
    for (size_t i=0; i<entries; ++i)
        averagedSpecList[i] = averageSpectraRaw(allSpectra[i], clMethod, window, true, minIntensity, minAbundance);

    Rcpp::List ret(entries);
    
    for (size_t i=0; i<entries; ++i)
    {
        if (alreadyAveraged)
        {
            ret[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = averagedSpecList[i].getMZs(),
                                             Rcpp::Named("intensity") = averagedSpecList[i].getIntensities(),
                                             Rcpp::Named("abundance") = averagedSpecList[i].getAbundances(),
                                             Rcpp::Named("abundance_prev") = averagedSpecList[i].getAvgPrevAbundances());
        }
        else
        {
            ret[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = averagedSpecList[i].getMZs(),
                                             Rcpp::Named("intensity") = averagedSpecList[i].getIntensities(),
                                             Rcpp::Named("abundance") = averagedSpecList[i].getAbundances());
        }
    }
    
    return ret;
}

std::vector<size_t> frameSubSpecIDs(const SpectrumRaw &frame)
{
    std::vector<size_t> ret(frame.size());
    size_t curID = 0;
    SpectrumRawTypes::Mobility curMob;
    for (size_t i=0; i<frame.size(); ++i)
    {
        const auto m = frame.getMobilities()[i];
        if (i == 0)
            curMob = m;
        else if (curMob != m)
        {
            ++curID;
            curMob = m;
        }
        ret[i] = curID;
    }
    return ret;
}
