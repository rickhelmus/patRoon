#ifndef PATROON_EIM_RUNNING_H
#define PATROON_EIM_RUNNING_H

#include "spectrum-raw.h"
#include "utils.hpp"

#include <deque>
#include <map>
#include <vector>

class EIMRunning
{
public:
    struct EIM
    {
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        EIM(void) = default;
        EIM(const std::vector<SpectrumRawTypes::Mobility> &m, const std::vector<SpectrumRawTypes::Intensity> &i)
            : mobilities(m), intensities(i) { }
        EIM(std::vector<SpectrumRawTypes::Mobility> &&m, std::vector<SpectrumRawTypes::Intensity> &&i)
            : mobilities(std::move(m)), intensities(std::move(i)) { }
        EIM(size_t s) : mobilities(s), intensities(s) { }
    };
    
private:
    const size_t maxQueueSize;
    std::deque<EIM> EIMs;
    
    void maybePop(void)
    {
        if (EIMs.size() > maxQueueSize)
            EIMs.pop_front();
    }
    
public:
    EIMRunning(size_t sz) : maxQueueSize(sz) { }
    
    void add(const std::vector<SpectrumRawTypes::Mobility> &mobilities,
             const std::vector<SpectrumRawTypes::Intensity> &intensities)
    {
        EIMs.emplace_back(mobilities, intensities);
        maybePop();
    }
    
    void add(std::vector<SpectrumRawTypes::Mobility> &&mobilities,
             std::vector<SpectrumRawTypes::Intensity> &&intensities)
    {
        EIMs.emplace_back(mobilities, intensities);
        maybePop();
    }
    
    void addZero(void)
    {
        EIMs.emplace_back();
        maybePop();
    }
    
    EIM get(unsigned smoothWindow, SpectrumRawTypes::Mobility mobStart, SpectrumRawTypes::Mobility mobEnd) const
    {
        std::map<SpectrumRawTypes::Mobility, SpectrumRawTypes::Intensity> merged;
        for (const auto& eim : EIMs)
        {
            for (size_t i = 0; i < eim.mobilities.size(); ++i)
                merged[eim.mobilities[i]] += eim.intensities[i];
        }
        
        EIM ret(merged.size());
        size_t i = 0;
        for (auto it=merged.begin(); it!=merged.end(); ++it, ++i)
        {
            ret.mobilities[i] = it->first;
            ret.intensities[i] = it->second;
        }
        
        if (smoothWindow > 0 && ret.mobilities.size() >= 3)
        {
            // pad to get smoothing right
            
            SpectrumRawTypes::Mobility minDiff = 0;
            for (size_t j=1; j<ret.mobilities.size(); ++j)
            {
                const auto diff = ret.mobilities[j] - ret.mobilities[j-1];
                if (minDiff == 0 || diff < minDiff)
                    minDiff = diff;
            }
            if (minDiff > 0.0)
            {
                auto mobMinDiff = ret.mobilities.front() - mobStart;
                for (unsigned j=0; j<smoothWindow && mobMinDiff>0.0; ++j)
                {
                    const auto m = ret.mobilities.front() - minDiff;
                    ret.mobilities.insert(ret.mobilities.begin(), m);
                    ret.intensities.insert(ret.intensities.begin(), 0);
                    mobMinDiff -= minDiff;
                }
                auto mobMaxDiff = mobEnd - ret.mobilities.back();
                for (unsigned j=0; j<smoothWindow && mobMaxDiff>0.0; ++j)
                {
                    const auto m = ret.mobilities.back() + minDiff;
                    ret.mobilities.push_back(m);
                    ret.intensities.push_back(0);
                    mobMaxDiff -= minDiff;
                }
            }
            ret.intensities = movingAverage(ret.intensities, smoothWindow);
        }
        
        return ret;
    }
    
    void pop(void)
    {
        if (!EIMs.empty())
            EIMs.pop_front();
    }
    
    bool empty(void) const { return EIMs.empty(); }
    size_t size(void) const { return EIMs.size(); }
    size_t sizeNoZero(void) const
    {
        size_t count = 0;
        for (const auto& eim : EIMs)
        {
            if (!eim.mobilities.empty())
                ++count;
        }
        return count;
    }
};

#endif // PATROON_EIM_RUNNING_H
