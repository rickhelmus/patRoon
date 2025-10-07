#ifndef MSDATA_EIC_H
#define MSDATA_EIC_H

#include <deque>
#include <map>
#include <memory>
#include <vector>

#include "spectrum-raw.h"

enum class EICMode { SIMPLE, FULL, FULL_MZ, TEST };

template <typename T> class IMSFrameSummer
{
    struct Frame
    {
        std::vector<T> xvalues;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        Frame(void) = default;
        Frame(std::vector<T> &&x, std::vector<SpectrumRawTypes::Intensity> &&i)
            : xvalues(std::move(x)), intensities(std::move(i)) { }
        size_t size(void) const { return xvalues.size(); }
    };

    const size_t maxQueueSize;
    std::deque<Frame> frames;

    void maybePop(void)
    {
        if (frames.size() > maxQueueSize)
            frames.pop_front();
    }

public:
    using XInts = std::pair<std::vector<T>, std::vector<SpectrumRawTypes::Intensity>>;

    IMSFrameSummer(size_t sz) : maxQueueSize(sz) { }

    void add(std::vector<T> &&xvalues, std::vector<SpectrumRawTypes::Intensity> &&intensities)
    {
        frames.emplace_back(std::move(xvalues), std::move(intensities));
        maybePop();
    }

    void addZero(void)
    {
        frames.emplace_back();
        maybePop();
    }

    XInts get(void) const
    {
        std::map<T, SpectrumRawTypes::Intensity> merged;
        for (const auto &fr : frames)
        {
            for (size_t j=0; j<fr.size(); ++j)
            {
                merged[fr.xvalues[j]] += fr.intensities[j];
            }
        }

        std::vector<T> xv(merged.size());
        std::vector<SpectrumRawTypes::Intensity> ints(merged.size());
        size_t i = 0;
        for (auto it=merged.begin(); it!=merged.end(); ++it, ++i)
        {
            xv[i] = it->first;
            ints[i] = it->second;
        }

        return std::make_pair(std::move(xv), std::move(ints));
    }

    void pop(void)
    {
        if (!frames.empty())
            frames.pop_front();
    }

    void clear(void) { frames.clear(); }

    bool empty(void) const { return frames.empty(); }
    size_t size(void) const { return frames.size(); }
    size_t sizeNoZero(void) const
    {
        size_t count = 0;
        for (const auto& fr : frames)
        {
            if (!fr.xvalues.empty())
                ++count;
        }
        return count;
    }
    size_t maxSize(void) const { return maxQueueSize; }
    size_t flank(void) const { return (maxQueueSize - 1) / 2; }
};

class EIC
{
    struct EICPoint
    {
        SpectrumRawTypes::Mass mz = 0, mzMin = 0, mzMax = 0, mzBP = 0;
        SpectrumRawTypes::Intensity intensity = 0, intensityBP = 0;
        SpectrumRawTypes::Mobility mobMin = 0, mobMax = 0;
        std::vector<SpectrumRawTypes::Mass> allMZs;
        std::vector<SpectrumRawTypes::Mobility> allMobs;
        std::vector<SpectrumRawTypes::Intensity> allIntsMZs, allIntsMobs;
        void clear(void)
        {
            mz = mzMin = mzMax = mzBP = 0;
            intensity = intensityBP = 0;
            mobMin = mobMax = 0;
            allMZs.clear();
            allMobs.clear();
            allIntsMZs.clear();
            allIntsMobs.clear();
        }
    };
    
    std::vector<SpectrumRawTypes::Scan> scanInds;
    std::vector<SpectrumRawTypes::Mass> mzs, mzMins, mzMaxs;
    std::vector<SpectrumRawTypes::Intensity> intensities;
    std::vector<SpectrumRawTypes::Mobility> mobilities, mobMins, mobMaxs;
    std::vector<SpectrumRawTypes::Mass> mzsBP;
    std::vector<SpectrumRawTypes::Mobility> mobilitiesBP;
    
    // for storing m/z profiles and EIMs for each point
    std::vector<std::pair<std::vector<SpectrumRawTypes::Mass>, std::vector<SpectrumRawTypes::Intensity>>> mzProfiles;
    std::vector<std::pair<std::vector<SpectrumRawTypes::Mobility>, std::vector<SpectrumRawTypes::Intensity>>> EIMs;
    
    EICMode mode;
    bool withMob;
    EICPoint curPoint;
    SpectrumRawTypes::Intensity maxIntensity = 0;
    
    // validity checks
    SpectrumRawTypes::Intensity minAdjIntensity;
    SpectrumRawTypes::Time minAdjTime;
    unsigned minAdjPoints;
    bool enoughTimeAboveThr = false;
    bool enoughPointsAboveThr = false;
    SpectrumRawTypes::Time startTimeAboveThr = 0.0;
    unsigned adjPointsAboveThr = 0;
    
    // IMS related
    IMSFrameSummer<SpectrumRawTypes::Mass> mzSummer;
    IMSFrameSummer<SpectrumRawTypes::Mobility> mobSummer;
    unsigned smoothWindowMZ, smoothWindowMob;
    SpectrumRawTypes::Mass smoothExtMZ;
    SpectrumRawTypes::Mobility smoothExtMob;
    SpectrumRawTypes::Mass mzStart, mzEnd;
    SpectrumRawTypes::Mobility mobStart, mobEnd;
    bool saveMZProfiles, saveEIMs;
    
    bool pointMZInRange(SpectrumRawTypes::Mass mz) const
    {
        return numberGTE(mz, mzStart) && (mzEnd == 0.0 || numberLTE(mz, mzEnd));
    }
    bool pointMobInRange(SpectrumRawTypes::Mobility mob) const
    {
        return (!withMob || (numberGTE(mob, mobStart) && (mobEnd == 0.0 || numberLTE(mob, mobEnd))));
    }
    
    size_t findEICIndexFromScan(SpectrumRawTypes::Scan scanInd) const;
    void setSummedFrameMZ(SpectrumRawTypes::Scan scanInd);
    void setSummedFrameMob(SpectrumRawTypes::Scan scanInd);
    void updateFrameSummer(void);
    void commitPoints(SpectrumRawTypes::Scan curScanInd);

public:
    EIC(EICMode em, bool wm, SpectrumRawTypes::Intensity minAdjI, SpectrumRawTypes::Time minAdjT, unsigned minAdjP,
        size_t sumFramesMZ, size_t sumFramesMob, unsigned smoMZ, unsigned smoMob, SpectrumRawTypes::Mass smoExtMZ,
        SpectrumRawTypes::Mobility smoExtMob, bool svMZPs,
        bool svEIMs) : mode(em), withMob(wm), minAdjIntensity(minAdjI), minAdjTime(minAdjT), minAdjPoints(minAdjP),
        mzSummer((withMob && (em == EICMode::FULL || em == EICMode::FULL_MZ)) ? sumFramesMZ : 0),
        mobSummer((withMob && mode == EICMode::FULL) ? sumFramesMob : 0), smoothWindowMZ(smoMZ),
        smoothWindowMob(smoMob), smoothExtMZ(smoExtMZ), smoothExtMob(smoExtMob),
        saveMZProfiles(svMZPs), saveEIMs(svEIMs) { }

    void setBoundaries(SpectrumRawTypes::Mass mzS, SpectrumRawTypes::Mass mzE,
                       SpectrumRawTypes::Mobility mobS, SpectrumRawTypes::Mobility mobE)
    {
        mzStart = mzS; mzEnd = mzE;
        mobStart = mobS; mobEnd = mobE;
    }
    
    void addPoint(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten);
    void addPoint(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Mobility mob, SpectrumRawTypes::Intensity inten);
    void commit(SpectrumRawTypes::Scan curScanInd, SpectrumRawTypes::Time curTime);
    void finalize(void);
    void pad(SpectrumRawTypes::Scan scanStart, SpectrumRawTypes::Scan scanEnd);
    
    bool verifyAdjancency(void) const
    {
        if (minAdjIntensity != 0.0)
        {
            if (minAdjTime > 0.0 && !enoughTimeAboveThr)
                return false;
            if (minAdjPoints > 0 && !enoughPointsAboveThr)
                return false;
        }
        return true;
    }
    
    void clear(void)
    {
        scanInds.clear();
        mzs.clear();
        mzMins.clear();
        mzMaxs.clear();
        intensities.clear();
        mobMins.clear();
        mobMaxs.clear();
        mzsBP.clear();
        mobilities.clear();
        mobilitiesBP.clear();
        mzProfiles.clear();
        EIMs.clear();
        curPoint.clear();
        maxIntensity = 0;
        enoughTimeAboveThr = false;
        enoughPointsAboveThr = false;
        startTimeAboveThr = 0.0;
        adjPointsAboveThr = 0;
        mzSummer.clear();
        mobSummer.clear();
    }
    
    size_t size(void) const { return scanInds.size(); }
    bool empty(void) const { return scanInds.empty(); }
    
    SpectrumRawTypes::Intensity getMaxIntensity(void) const { return maxIntensity; }
    
    auto getScanIndices(void) const { return scanInds; }
    auto getMZs(void) const { return mzs; }
    auto getMZMins(void) const { return mzMins; }
    auto getMZMaxs(void) const { return mzMaxs; }
    auto getIntensities(void) const { return intensities; }
    auto getMobMins(void) const { return mobMins; }
    auto getMobMaxs(void) const { return mobMaxs; }
    auto getMobilities(void) const { return mobilities; }
    auto getMZsBP(void) const { return mzsBP; }
    auto getMobilitiesBP(void) const { return mobilitiesBP; }
    auto getMZProfiles(void) const { return mzProfiles; }
    auto getEIMs(void) const { return EIMs; }
};


#endif // MSDATA_EIC_H
