/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

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
        SpectrumRawTypes::Time time = 0.0;
        std::vector<T> xvalues;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        Frame(void) = default;
        Frame(SpectrumRawTypes::Time t, std::vector<T> &&x, std::vector<SpectrumRawTypes::Intensity> &&i)
            : time(t), xvalues(std::move(x)), intensities(std::move(i)) { }
        size_t size(void) const { return xvalues.size(); }
    };

    std::deque<Frame> frames;

public:
    using XInts = std::pair<std::vector<T>, std::vector<SpectrumRawTypes::Intensity>>;

    void add(SpectrumRawTypes::Time time, std::vector<T> &&xvalues,
             std::vector<SpectrumRawTypes::Intensity> &&intensities)
    {
        frames.emplace_back(time, std::move(xvalues), std::move(intensities));
    }

    XInts get(void) const
    {
        std::map<T, SpectrumRawTypes::Intensity> merged;
        for (const auto &fr : frames)
        {
            // Rcpp::Rcout << "  IMS frame: " << fr.size() << " points at time " << fr.time << std::endl;
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
    SpectrumRawTypes::Time frontTime(void) const { return frames.empty() ? 0.0 : frames.front().time; }
    SpectrumRawTypes::Time backTime(void) const { return frames.empty() ? 0.0 : frames.back().time; }
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
    
    const std::vector<SpectrumRawTypes::Time> &EICTimes;
    const EICMode mode;
    const bool withMob;
    EICPoint curPoint;
    SpectrumRawTypes::Intensity maxIntensity = 0;
    
    // validity checks
    const SpectrumRawTypes::Intensity minAdjIntensity;
    const SpectrumRawTypes::Time minAdjTime;
    const unsigned minAdjPoints;
    bool enoughTimeAboveThr = false;
    bool enoughPointsAboveThr = false;
    SpectrumRawTypes::Time startTimeAboveThr = 0.0;
    unsigned adjPointsAboveThr = 0;
    
    // IMS related
    IMSFrameSummer<SpectrumRawTypes::Mass> mzSummer;
    IMSFrameSummer<SpectrumRawTypes::Mobility> mobSummer;
    const SpectrumRawTypes::Time sumWindowMZ, sumWindowMob;
    size_t updateSumMZIndex = 0, updateSumMobIndex = 0;
    const unsigned smoothWindowMZ, smoothWindowMob;
    const SpectrumRawTypes::Mass smoothExtMZ;
    const SpectrumRawTypes::Mobility smoothExtMob;
    SpectrumRawTypes::Mass mzStart, mzEnd;
    SpectrumRawTypes::Mobility mobStart, mobEnd;
    const bool saveMZProfiles, saveEIMs;
    
    bool pointMZInRange(SpectrumRawTypes::Mass mz) const
    {
        return numberGTE(mz, mzStart) && (mzEnd == 0.0 || numberLTE(mz, mzEnd));
    }
    bool pointMobInRange(SpectrumRawTypes::Mobility mob) const
    {
        return (!withMob || (numberGTE(mob, mobStart) && (mobEnd == 0.0 || numberLTE(mob, mobEnd))));
    }
    
    SpectrumRawTypes::Time getUpdateSumMZTime(void) const
    {
        return (updateSumMZIndex >= scanInds.size()) ? 0.0 : EICTimes[scanInds[updateSumMZIndex]];
    }
    SpectrumRawTypes::Time getUpdateSumMobTime(void) const
    {
        return (updateSumMobIndex >= scanInds.size()) ? 0.0 : EICTimes[scanInds[updateSumMobIndex]];
    }
    
    void setSummedFrameMZ(void);
    void setSummedFrameMob(void);
    void updateFrameSummer(SpectrumRawTypes::Time curTime);
    void commitPoints(SpectrumRawTypes::Scan curScanInd);

public:
    EIC(const std::vector<SpectrumRawTypes::Time> &times, EICMode em, bool wm, SpectrumRawTypes::Intensity minAdjI,
        SpectrumRawTypes::Time minAdjT, unsigned minAdjP, SpectrumRawTypes::Time sumMZ, SpectrumRawTypes::Time sumMob,
        unsigned smoMZ, unsigned smoMob, SpectrumRawTypes::Mass smoExtMZ, SpectrumRawTypes::Mobility smoExtMob,
        bool svMZPs, bool svEIMs) : EICTimes(times), mode(em), withMob(wm), minAdjIntensity(minAdjI),
        minAdjTime(minAdjT), minAdjPoints(minAdjP), sumWindowMZ(sumMZ), sumWindowMob(sumMob), smoothWindowMZ(smoMZ),
        smoothWindowMob(smoMob), smoothExtMZ(smoExtMZ), smoothExtMob(smoExtMob), saveMZProfiles(svMZPs),
        saveEIMs(svEIMs) { }

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
    
    void clear(bool keepMZProfiles = false, bool keepEIMs = false)
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
        curPoint.clear();
        maxIntensity = 0;
        enoughTimeAboveThr = false;
        enoughPointsAboveThr = false;
        startTimeAboveThr = 0.0;
        adjPointsAboveThr = 0;
        mzSummer.clear();
        mobSummer.clear();
        updateSumMZIndex = updateSumMobIndex = 0;
        if (!keepMZProfiles)
            mzProfiles.clear();
        if (!keepEIMs)
            EIMs.clear();
    }
    
    size_t size(void) const { return scanInds.size(); }
    bool empty(void) const { return scanInds.empty(); }
    
    SpectrumRawTypes::Intensity getMaxIntensity(void) const { return maxIntensity; }
    
    const auto &getScanIndices(void) const { return scanInds; }
    const auto &getMZs(void) const { return mzs; }
    const auto &getMZMins(void) const { return mzMins; }
    const auto &getMZMaxs(void) const { return mzMaxs; }
    const auto &getIntensities(void) const & { return intensities; }
    const auto &getMobMins(void) const { return mobMins; }
    const auto &getMobMaxs(void) const { return mobMaxs; }
    const auto &getMobilities(void) const { return mobilities; }
    const auto &getMZsBP(void) const { return mzsBP; }
    const auto &getMobilitiesBP(void) const { return mobilitiesBP; }
    const auto &getMZProfiles(void) const { return mzProfiles; }
    const auto &getEIMs(void) const { return EIMs; }
};

class SimpleEIC
{
    std::vector<SpectrumRawTypes::Time> times;
    std::vector<SpectrumRawTypes::Intensity> intensities;
    
public:
    SimpleEIC(void) = default;
    SimpleEIC(const EIC &eic, const std::vector<SpectrumRawTypes::Time> &scanTimes) : times(eic.size()), intensities(eic.getIntensities())
    {
        for (size_t i=0; i<eic.size(); ++i)
            times[i] = scanTimes[eic.getScanIndices()[i]];
    }
    
    size_t size(void) const { return times.size(); }
    bool empty(void) const { return times.empty(); }
    
    void fillGaps(SpectrumRawTypes::Time medRTDiff, double gapFactor, bool pad, SpectrumRawTypes::Time startTime,
                  SpectrumRawTypes::Time endTime);
    
    const auto &getTimes(void) const { return times; }
    const auto &getIntensities(void) const { return intensities; }
};


#endif // MSDATA_EIC_H
