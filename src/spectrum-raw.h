#ifndef PATROON_SPECTRUM_RAW_H
#define PATROON_SPECTRUM_RAW_H

#include <stddef.h>
#include <vector>

#include "utils.h"

namespace SpectrumRawTypes {

using Scan = unsigned;
using Time = float;
using Mass = float;
using Intensity = float; // UNDONE: or unsigned?
using IsolationRange = NumRange<Mass>;

enum class MSLevel { MS1, MS2 };

}

struct SpectrumRawMetadataMS
{
    std::vector<SpectrumRawTypes::Scan> scans;
    std::vector<SpectrumRawTypes::Time> times;
    std::vector<SpectrumRawTypes::Intensity> TICs, BPCs;
    void clear(void) { scans.clear(); times.clear(); TICs.clear(); BPCs.clear(); }
};
struct SpectrumRawMetadataMSMS: public SpectrumRawMetadataMS
{
    std::vector<SpectrumRawTypes::IsolationRange> isolationRanges;
    void clear(void) { SpectrumRawMetadataMS::clear(); isolationRanges.clear(); } // NOTE: not virtual
};
using SpectrumRawMetadata = std::pair<SpectrumRawMetadataMS, SpectrumRawMetadataMSMS>;

class SpectrumRaw
{
    SpectrumRawTypes::Time time;
    std::vector<SpectrumRawTypes::Mass> mzs;
    std::vector<SpectrumRawTypes::Intensity> intensities;
    
public:
    SpectrumRaw(void) = default;
    SpectrumRaw(SpectrumRawTypes::Time t, const std::vector<SpectrumRawTypes::Mass> &m,
                const std::vector<SpectrumRawTypes::Intensity> &i) : time(t), mzs(m), intensities(i) { }
    SpectrumRaw(SpectrumRawTypes::Time t, size_t size) : time(t), mzs(size), intensities(size) { }
    
    SpectrumRawTypes::Time getTime(void) const { return time; }
    const auto &getMZs(void) const { return mzs; }
    const auto &getIntensities(void) const { return intensities; }
    
    void setTime(SpectrumRawTypes::Time t) { time = t; }
    void append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten);
    void append(const SpectrumRaw &sp);
    void setPeak(size_t i, SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten) { mzs[i] = mz; intensities[i] = inten; }
    
    auto size(void) const { return mzs.size(); }
    bool empty(void) const { return mzs.empty(); }
    void resize(size_t s) { mzs.resize(s); intensities.resize(s); }
    void clear(void) { time = 0.0f; mzs.clear(); intensities.clear(); }
};

std::vector<SpectrumRawTypes::Scan> getSpecScanIndices(const SpectrumRawMetadata &specMeta,
                                                       const NumRange<SpectrumRawTypes::Time> &timeRange,
                                                       SpectrumRawTypes::MSLevel MSLevel,
                                                       const NumRange<SpectrumRawTypes::Mass> &isoRange);
#endif
