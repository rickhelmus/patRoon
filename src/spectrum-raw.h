#ifndef PATROON_SPECTRUM_RAW_H
#define PATROON_SPECTRUM_RAW_H

#include <vector>

class SpectrumRaw
{
public:
    using NumVecType = std::vector<double>;
    
private:
    // UNDONE: should intensities be double or int?
    float time;
    NumVecType mzs, intensities;
    
public:
    SpectrumRaw() = default;
    SpectrumRaw(float t, const NumVecType &m, const NumVecType &i) : time(t), mzs(m), intensities(i) { }
    SpectrumRaw(NumVecType::size_type size) : mzs(size), intensities(size) { }
    
    float getTime(void) const { return time; }
    const auto &getMZs(void) const { return mzs; }
    const auto &getIntensities(void) const { return intensities; }
    
    void setTime(float t) { time = t; }
    void append(double mz, double inten);
    void append(const SpectrumRaw &sp);
    void setPeak(NumVecType::size_type i, double mz, unsigned inten) { mzs[i] = mz; intensities[i] = inten; }
    
    auto size(void) const { return mzs.size(); }
    bool empty(void) const { return mzs.empty(); }
    void resize(NumVecType::size_type s) { mzs.resize(s); intensities.resize(s); }
    void clear(void) { mzs.clear(); intensities.clear(); }
};

#endif
