#ifndef MSDATA_EIC_H
#define MSDATA_EIC_H

#include <deque>
#include <vector>

#include "spectrum-raw.h"

enum class EICMode { SIMPLE, FULL, FULL_MZ, TEST };

class IMSFrameSummer
{
    struct Frame
    {
        std::vector<SpectrumRawTypes::Mass> mzs;
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        Frame(void) = default;
        Frame(std::vector<SpectrumRawTypes::Mass> &&m, std::vector<SpectrumRawTypes::Intensity> &&i)
            : mzs(std::move(m)), intensities(std::move(i)) { }
        Frame(std::vector<SpectrumRawTypes::Mass> &&m, std::vector<SpectrumRawTypes::Mobility> &&mob,
              std::vector<SpectrumRawTypes::Intensity> &&i)
            : mzs(std::move(m)), mobilities(mob), intensities(std::move(i)) { }
        size_t size(void) const { return mzs.size(); }
    };
    
    const size_t maxQueueSize;
    std::deque<Frame> frames;
    
    void maybePop(void)
    {
        if (frames.size() > maxQueueSize)
            frames.pop_front();
    }
    
    template <typename T, typename U>
    std::pair<std::vector<T>, std::vector<SpectrumRawTypes::Intensity>> get(U Frame::* memberPtr) const
    {
        std::map<T, SpectrumRawTypes::Intensity> merged;
        for (const auto &fr : frames)
        {
            for (size_t i=0; i<fr.size(); ++i)
                merged[(fr.*memberPtr)[i]] += fr.intensities[i];
        }
        
        std::vector<T> xvalues(merged.size());
        std::vector<SpectrumRawTypes::Intensity> ints(merged.size());
        size_t i = 0;
        for (auto it=merged.begin(); it!=merged.end(); ++it, ++i)
        {
            xvalues[i] = it->first;
            ints[i] = it->second;
        }
        
        return std::make_pair(std::move(xvalues), std::move(ints));
    }
    
public:
    using MassInts = std::pair<std::vector<SpectrumRawTypes::Mass>, std::vector<SpectrumRawTypes::Intensity>>;
    using MobInts = std::pair<std::vector<SpectrumRawTypes::Mobility>, std::vector<SpectrumRawTypes::Intensity>>;
    
    IMSFrameSummer(size_t sz) : maxQueueSize(sz) { }
    
    void add(std::vector<SpectrumRawTypes::Mass> &&mzs, std::vector<SpectrumRawTypes::Intensity> &&intensities)
    {
        frames.emplace_back(std::move(mzs), std::move(intensities));
        maybePop();
    }
    
    void add(std::vector<SpectrumRawTypes::Mass> &&mzs, std::vector<SpectrumRawTypes::Mobility> &&mobilities,
             std::vector<SpectrumRawTypes::Intensity> &&intensities)
    {
        frames.emplace_back(std::move(mzs), std::move(mobilities), std::move(intensities));
        maybePop();
    }
    
    void addZero(void)
    {
        frames.emplace_back();
        maybePop();
    }
    
    MassInts getMZs(void) const
    {
        return get<SpectrumRawTypes::Mass>(&Frame::mzs);
    }
    
    MobInts getMobilities(void) const
    {
        return get<SpectrumRawTypes::Mobility>(&Frame::mobilities);
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
            if (!fr.mobilities.empty())
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
        std::vector<SpectrumRawTypes::Intensity> allInts;
        void clear(void)
        {
            mz = mzMin = mzMax = mzBP = 0;
            intensity = intensityBP = 0;
            mobMin = mobMax = 0;
            allMZs.clear();
            allMobs.clear();
            allInts.clear();
        }
    };
    
    std::vector<SpectrumRawTypes::Scan> scanInds;
    std::vector<SpectrumRawTypes::Mass> mzs, mzMins, mzMaxs;
    std::vector<SpectrumRawTypes::Intensity> intensities;
    std::vector<SpectrumRawTypes::Mobility> mobilities, mobMins, mobMaxs;
    std::vector<SpectrumRawTypes::Mass> mzsBP;
    std::vector<SpectrumRawTypes::Mobility> mobilitiesBP;
    std::vector<SpectrumRawTypes::Intensity> intensitiesBP;
    
    // for storing m/z profiles and EIMs for each point
    std::vector<std::pair<std::vector<SpectrumRawTypes::Mass>, std::vector<SpectrumRawTypes::Intensity>>> MZProfiles;
    std::vector<std::pair<std::vector<SpectrumRawTypes::Mobility>, std::vector<SpectrumRawTypes::Intensity>>> EIMs;
    
    EICMode mode;
    bool withMob;
    EICPoint curPoint;
    
    // validity checks
    SpectrumRawTypes::Intensity minEICAdjIntensity = 0.0;
    SpectrumRawTypes::Time minEICAdjTime = 0.0;
    unsigned minEICAdjPoints = 0;
    bool enoughTimeAboveThr = false;
    bool enoughPointsAboveThr = false;
    SpectrumRawTypes::Time startTimeAboveThr = 0.0;
    unsigned adjPointsAboveThr = 0;
    
    // IMS related
    IMSFrameSummer frameSummer;
    unsigned smoothWindowMZ, smoothWindowMob;
    SpectrumRawTypes::Mass smoothExtMZ;
    SpectrumRawTypes::Mobility smoothExtMob;
    SpectrumRawTypes::Mass mzStart, mzEnd;
    SpectrumRawTypes::Mobility mobStart, mobEnd;
    bool saveMZProfiles, saveEIMs;
    
    void setSummedFrame(SpectrumRawTypes::Scan scanInd);
    void updateFrameSummer(SpectrumRawTypes::Scan curScanInd);
    void commitPoints(SpectrumRawTypes::Scan curScanInd);

public:
    EIC(EICMode em, bool wm, size_t sumFrames, unsigned smoMZ, unsigned smoMob, SpectrumRawTypes::Mass smoExtMZ,
        SpectrumRawTypes::Mobility smoExtMob, bool svMZPs,
        bool svEIMs) : mode(em), withMob(wm), frameSummer((withMob && (mode == EICMode::FULL || mode == EICMode::FULL_MZ)) ? sumFrames : 0),
        smoothWindowMZ(smoMZ), smoothWindowMob(smoMob), smoothExtMZ(smoExtMZ), smoothExtMob(smoExtMob),
        saveMZProfiles(svMZPs), saveEIMs(svEIMs) { }
    
    void setBoundaries(SpectrumRawTypes::Mass mzS, SpectrumRawTypes::Mass mzE,
                       SpectrumRawTypes::Mobility mobS, SpectrumRawTypes::Mobility mobE)
    {
        mzStart = mzS; mzEnd = mzE;
        mobStart = mobS; mobEnd = mobE;
    }
    
    void addPoint(SpectrumRawTypes::Intensity inten);
    void addPoint(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten);
    void addPoint(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Mobility mob, SpectrumRawTypes::Intensity inten);
    void commit(SpectrumRawTypes::Scan curScanInd, SpectrumRawTypes::Time curTime);
    void finalize(void);
    void pad(SpectrumRawTypes::Scan scanStart, SpectrumRawTypes::Scan scanEnd);
    
    bool valid(void) const
    {
        if (minEICAdjTime > 0.0 && !enoughTimeAboveThr)
            return false;
        if (minEICAdjPoints > 0 && !enoughPointsAboveThr)
            return false;
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
        intensitiesBP.clear();
        mobilities.clear();
        mobilitiesBP.clear();
        MZProfiles.clear();
        EIMs.clear();
        curPoint.clear();
        enoughTimeAboveThr = false;
        enoughPointsAboveThr = false;
        startTimeAboveThr = 0.0;
        adjPointsAboveThr = 0;
        frameSummer.clear();
    }
    
    size_t size(void) const { return scanInds.size(); }
    bool empty(void) const { return scanInds.empty(); }
    
    /* re-use eic in OpenMP loop
     * set capacity in constructor to max size
     * store as vector of structs
    */
};


#endif // MSDATA_EIC_H
