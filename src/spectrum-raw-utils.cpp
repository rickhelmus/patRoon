#include <Rcpp.h>

#include "spectrum-raw.h"

#include <algorithm>

std::vector<SpectrumRawTypes::Scan> getSpecScanIndices(const SpectrumRawMetadata &specMeta,
                                                       const NumRange<SpectrumRawTypes::Time> &timeRange,
                                                       SpectrumRawTypes::MSLevel MSLevel,
                                                       const NumRange<SpectrumRawTypes::Mass> &isoRange)
{
    const SpectrumRawMetadataMS &metaMS = (MSLevel == SpectrumRawTypes::MSLevel::MS1) ? specMeta.first : specMeta.second;
    std::vector<SpectrumRawTypes::Scan> ret;
    
    const auto startIt = std::lower_bound(metaMS.times.begin(), metaMS.times.end(), timeRange.start);
    if (startIt != metaMS.times.end())
    {
        for (size_t i=std::distance(metaMS.times.begin(), startIt);
             i<metaMS.times.size() && metaMS.times[i]<=timeRange.end; ++i)
        {
            if (MSLevel == SpectrumRawTypes::MSLevel::MS2 && !isoRange.inside(specMeta.second.isolationRanges[i]))
                continue;
            ret.push_back(metaMS.scans[i]);
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
            if (minPrecMZDiff == -1.0 || d < minPrecMZDiff)
            {
                closestPrecMZInd = std::distance(spectrum.getMZs().begin(), it);
                minPrecMZDiff = d;
            }
            else
                break; // data is sorted: all next mzs will be greater and therefore with higher deviation
        }
    }
    
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
    
    if (doTopMost)
    {
        std::vector<SpectrumRawTypes::Intensity> topMostInts(filter.topMost);
        std::partial_sort_copy(spectrum.getIntensities().begin() + startInd, spectrum.getIntensities().begin() + endInd,
                               topMostInts.begin(), topMostInts.end(),
                               [](SpectrumRawTypes::Intensity i1, SpectrumRawTypes::Intensity i2) { return i1 > i2; });
        minInt = std::max(minInt, topMostInts.back());
        Rcpp::Rcout << "tm: " << topMostInts.back() << "/" << topMostInts.size() << "/" << topMostInts.front() << "\n";
    }
    
    bool addedPrec = false;
    for (size_t i=startInd; i<endInd; i++)
    {
        if (spectrum.getIntensities()[i] >= minInt)
        {
            ret.append(spectrum.getMZs()[i], spectrum.getIntensities()[i]);
            if (i == closestPrecMZInd)
                addedPrec = true;
        }
    }
    
    if (!addedPrec && closestPrecMZInd != spectrum.size())
    {
        const auto it = std::upper_bound(ret.getMZs().begin(), ret.getMZs().end(), spectrum.getMZs()[closestPrecMZInd]);
        ret.insert(std::distance(ret.getMZs().begin(), it), spectrum.getMZs()[closestPrecMZInd],
                   spectrum.getIntensities()[closestPrecMZInd]);
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
