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
