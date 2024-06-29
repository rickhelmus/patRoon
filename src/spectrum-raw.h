#ifndef PATROON_SPECTRUM_RAW_H
#define PATROON_SPECTRUM_RAW_H

#include <stddef.h>
#include <vector>

namespace SpectrumRawTypes {

using Time = float;
using Mass = float;
using Intensity = float; // UNDONE: or unsigned?

}

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

#endif
