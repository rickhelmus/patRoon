#ifndef PATROON_SPECTRUM_RAW_H
#define PATROON_SPECTRUM_RAW_H

#include <stddef.h>
#include <algorithm>
#include <vector>

#include "utils.h"

namespace SpectrumRawTypes {

using Scan = unsigned;
using Time = float;
using Mass = float;
using Mobility = float;
using Intensity = float; // UNDONE: or unsigned?
using TimeRange = NumRange<Time>;
using IsolationRange = NumRange<Mass>;
using MobilityRange = NumRange<Mobility>;
using PeakAbundance = float;

enum class MSLevel { MS1, MS2 };
enum class MSPolarity { POSITIVE = 1, NEGATIVE = -1, UNKNOWN = 0 };

}

struct frameMSMSInfo
{
    std::vector<SpectrumRawTypes::IsolationRange> isolationRanges;
    std::vector<SpectrumRawTypes::Scan> subScans, subScanEnds;
    bool empty(void) const { return isolationRanges.empty(); }
    void clear(void) { isolationRanges.clear(); subScans.clear(); subScanEnds.clear(); }
};

struct SpectrumRawMetadataMS
{
    std::vector<SpectrumRawTypes::Scan> scans;
    std::vector<SpectrumRawTypes::Time> times;
    std::vector<SpectrumRawTypes::Intensity> TICs, BPCs;
    std::vector<SpectrumRawTypes::MSPolarity> polarities;
    SpectrumRawMetadataMS(void) = default;
    SpectrumRawMetadataMS(size_t s) : scans(s), times(s), TICs(s), BPCs(s), polarities(s) { }
    void clear(void) { scans.clear(); times.clear(); TICs.clear(); BPCs.clear(); polarities.clear(); }
    void append(const SpectrumRawMetadataMS &other)
    {
        scans.insert(scans.end(), other.scans.begin(), other.scans.end());
        times.insert(times.end(), other.times.begin(), other.times.end());
        TICs.insert(TICs.end(), other.TICs.begin(), other.TICs.end());
        BPCs.insert(BPCs.end(), other.BPCs.begin(), other.BPCs.end());
        polarities.insert(polarities.end(), other.polarities.begin(), other.polarities.end());
    }
};
struct SpectrumRawMetadataMSMS: public SpectrumRawMetadataMS
{
    std::vector<SpectrumRawTypes::IsolationRange> isolationRanges;
    std::vector<frameMSMSInfo> MSMSFrames;
    SpectrumRawMetadataMSMS(void) = default;
    SpectrumRawMetadataMSMS(size_t s, bool ims = false) : SpectrumRawMetadataMS(s),
        isolationRanges((!ims) ? s : 0), MSMSFrames((ims) ? s : 0) { }
    
    // NOTE: not virtual, to keep the class simple
    void clear(void) { SpectrumRawMetadataMS::clear(); isolationRanges.clear(); MSMSFrames.clear(); }
    void append(const SpectrumRawMetadataMSMS &other)
    {
        SpectrumRawMetadataMS::append(other);
        isolationRanges.insert(isolationRanges.end(), other.isolationRanges.begin(), other.isolationRanges.end());
        MSMSFrames.insert(MSMSFrames.end(), other.MSMSFrames.begin(), other.MSMSFrames.end());
    }
};
using SpectrumRawMetadata = std::pair<SpectrumRawMetadataMS, SpectrumRawMetadataMSMS>;

class SpectrumRaw
{
    std::vector<SpectrumRawTypes::Mass> mzs;
    std::vector<SpectrumRawTypes::Intensity> intensities;
    std::vector<SpectrumRawTypes::Mobility> mobilities;
    
public:
    SpectrumRaw(void) = default;
    SpectrumRaw(const std::vector<SpectrumRawTypes::Mass> &m,
                const std::vector<SpectrumRawTypes::Intensity> &i) : mzs(m), intensities(i) { }
    SpectrumRaw(size_t size, bool mobs = false) : mzs(size), intensities(size), mobilities((mobs) ? size : 0) { }
    
    const auto &getMZs(void) const { return mzs; }
    const auto &getIntensities(void) const { return intensities; }
    const auto &getMobilities(void) const { return mobilities; }
    
    bool hasMobilities(void) const { return !mobilities.empty(); }
    
    void append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten);
    void append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten, SpectrumRawTypes::Mobility mob);
    void append(const SpectrumRaw &sp);
    void insert(size_t i, SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten);
    void insert(size_t i, SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten, SpectrumRawTypes::Mobility mob);
    void setPeak(size_t i, SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten) { mzs[i] = mz; intensities[i] = inten; }
    void setPeak(size_t i, SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten,
                 SpectrumRawTypes::Mobility mob) { mzs[i] = mz; intensities[i] = inten; mobilities[i] = mob; }
    
    auto size(void) const { return mzs.size(); }
    bool empty(void) const { return mzs.empty(); }
    void clear(void) { mzs.clear(); intensities.clear(); mobilities.clear(); }
};

class SpectrumRawAveraged: public SpectrumRaw
{
    std::vector<SpectrumRawTypes::PeakAbundance> abundances, averagedPreviousAbundance;
    
public:
    SpectrumRawAveraged(void) = default;
    SpectrumRawAveraged(const std::vector<SpectrumRawTypes::Mass> &m,
                        const std::vector<SpectrumRawTypes::Intensity> &i) : SpectrumRaw(m, i) { }
    SpectrumRawAveraged(size_t size, bool ab = false, bool aab = false) : SpectrumRaw(size),
                                                                          abundances((ab) ? size : 0),
                                                                          averagedPreviousAbundance((aab) ? size : 0) { }
    
    const auto &getAbundances(void) const { return abundances; }
    const auto &getAvgPrevAbundances(void) const { return averagedPreviousAbundance; }
    
    using SpectrumRaw::append;
    void append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten, SpectrumRawTypes::PeakAbundance ab);
    void append(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten, SpectrumRawTypes::PeakAbundance ab,
                SpectrumRawTypes::PeakAbundance avgPrvAb);
    void append(const SpectrumRawAveraged &sp);
    
    void setAvgPrevAbundance(size_t i, SpectrumRawTypes::PeakAbundance ab) { averagedPreviousAbundance[i] = ab; }
};

struct SpectrumRawFilter
{
    SpectrumRawTypes::Intensity minIntensity = 0;
    NumRange<SpectrumRawTypes::Mass> mzRange;
    unsigned topMost = 0;
    bool withPrecursor = false, retainPrecursor = false;
    
    SpectrumRawFilter &setMinIntensity(unsigned i) { minIntensity = i; return *this; }
    SpectrumRawFilter &setMZRange(SpectrumRawTypes::Mass s, SpectrumRawTypes::Mass e) { mzRange.set(s, e); return *this; }
    SpectrumRawFilter &setTopMost(unsigned t) { topMost = t; return *this; }
    SpectrumRawFilter &setWithPrecursor(bool p) { withPrecursor = p; return *this; }
    SpectrumRawFilter &setRetainPrecursor(bool p) { retainPrecursor = p; return *this; }
};

struct SpectrumRawSelection
{
    SpectrumRawTypes::Scan index;
    std::vector<size_t> MSMSFrameIndices;
    SpectrumRawSelection(void) = default;
    SpectrumRawSelection(SpectrumRawTypes::Scan i) : index(i) { }
    bool operator==(const SpectrumRawSelection &o) const { return index == o.index && MSMSFrameIndices == o.MSMSFrameIndices; }
    bool operator!=(const SpectrumRawSelection &o) const { return index != o.index || MSMSFrameIndices != o.MSMSFrameIndices; }
};

std::vector<SpectrumRawSelection> getSpecRawSelections(const SpectrumRawMetadata &specMeta,
                                                       const NumRange<SpectrumRawTypes::Time> &timeRange,
                                                       SpectrumRawTypes::MSLevel MSLevel,
                                                       SpectrumRawTypes::Mass precursor,
                                                       SpectrumRawTypes::Intensity minBPIntensity = 0);
SpectrumRaw filterSpectrumRaw(const SpectrumRaw &spectrum, const SpectrumRawFilter &filter,
                              SpectrumRawTypes::Mass precursor);
SpectrumRaw filterIMSFrame(const SpectrumRaw &spectrum, const SpectrumRawFilter &filter,
                           SpectrumRawTypes::Mass precursor, const SpectrumRawTypes::MobilityRange &mobRange);

SpectrumRawAveraged averageSpectraRaw(const SpectrumRawAveraged &flattenedSpecs, const std::vector<size_t> &specIDs,
                                      clusterMethod method, SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance);
SpectrumRawAveraged averageSpectraRaw(const SpectrumRaw &flattenedSpecs, const std::vector<size_t> &specIDs,
                                      clusterMethod method, SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance);
SpectrumRawAveraged averageSpectraRaw(const std::vector<SpectrumRaw> &spectra, clusterMethod method,
                                      SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance);
SpectrumRawAveraged averageSpectraRaw(const std::vector<SpectrumRawAveraged> &spectra, clusterMethod method,
                                      SpectrumRawTypes::Mass window, bool averageIntensities,
                                      SpectrumRawTypes::Intensity minIntensity,
                                      SpectrumRawTypes::PeakAbundance minAbundance);
std::vector<size_t> frameSubSpecIDs(const SpectrumRaw &frame);

#endif
