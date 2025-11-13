/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <Rcpp.h>

#include <algorithm>
#include <chrono>
#include <iomanip>

#include "msdata.hpp"
#include "msdata-eic.h"

// #define LOG_EIC

namespace {

class SimpleTimer
{
    using clock = std::chrono::steady_clock;
    
    clock::time_point startTime;
    bool running, enabled;
    
public:
    SimpleTimer(bool e = true) : running(false), enabled(e) {}
    
    void start(const std::string &msg = "")
    {
        if (!enabled)
            return;
        
        if (!msg.empty())
            Rcpp::Rcout << msg << ": ";
        startTime = clock::now();
        running = true;
    }
    
    double stop(void)
    {
        if (!running || !enabled)
            return 0.0;
        const auto end = clock::now();
        running = false;
        const double ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - startTime).count();
        
        if (ms >= 1000.0)
            Rcpp::Rcout << std::fixed << std::setprecision(3) << (ms / 1000.0) << " s\n";
        else
            Rcpp::Rcout << std::fixed << std::setprecision(3) << ms << " ms\n";
        
        return ms;
    }
};
    
SpectrumRawTypes::Time medianRTDiff(const std::vector<SpectrumRawTypes::Time> &times)
{
    if (times.size() < 2)
        return 0.0;
    std::vector<SpectrumRawTypes::Time> diffs(times.size() - 1, 0);
    for (size_t i=1; i<times.size(); ++i)
        diffs[i - 1] = times[i] - times[i - 1];
    return median(diffs);
}

std::pair<std::vector<size_t>, std::vector<SpectrumRawTypes::Time>> fillRTGaps(const std::vector<SpectrumRawTypes::Time> &EICTimes,
                                                                               SpectrumRawTypes::Time medRTDiff,
                                                                               double gapFactor)
{
    // NOTE: Bruker TIMS data (and maybe others?) seem to omit zero intensity scans, leading to time gaps. Since
    // this leads to incorrect EICs, we pad here. To detect gaps, we take the median difference between scans
    // and assume a gap is at least X higher than that.
    
    
    if (EICTimes.empty())
        return std::make_pair(std::vector<size_t>(), std::vector<SpectrumRawTypes::Time>());
    if (EICTimes.size() == 1)
        return std::make_pair(std::vector<size_t>(1, 0), EICTimes);
    
    std::vector<size_t> filledInds;
    std::vector<SpectrumRawTypes::Time> filledTimes;
    const size_t fakeIndex = EICTimes.size();
    
    for (size_t i=0; i<EICTimes.size(); ++i)
    {
        filledInds.push_back(i);
        filledTimes.push_back(EICTimes[i]);
        if (i == (EICTimes.size() - 1))
            break;
        
        const auto diff = EICTimes[i + 1] - EICTimes[i];
        if (diff > (medRTDiff * gapFactor)) // add dummy point after current
        {
            filledInds.push_back(fakeIndex);
            filledTimes.push_back(EICTimes[i] + medRTDiff);
            if (diff > (medRTDiff * 2.0)) // add dummy point before next
            {
                filledInds.push_back(fakeIndex);
                filledTimes.push_back(EICTimes[i + 1] - medRTDiff);
            }
        }
    }
    
    return std::make_pair(std::move(filledInds), std::move(filledTimes));
}

template <typename T>
void smoothSummedFrameData(std::pair<std::vector<T>, std::vector<SpectrumRawTypes::Intensity>> &data,
                           unsigned window, T start, T end)
{
    if (data.first.size() < 3 || window < 1)
        return;
    
    // pad to get smoothing right
    
    T minDiff = 0;
    for (size_t i=1; i<data.first.size(); ++i)
    {
        const auto diff = data.first[i] - data.first[i-1];
        if (minDiff == 0 || diff < minDiff)
            minDiff = diff;
    }
    if (minDiff > 0.0)
    {
        auto startDiff = data.first.front() - start;
        for (unsigned i=0; i<window && startDiff>0.0; ++i)
        {
            const auto m = std::max(data.first.front() - minDiff, start);
            data.first.insert(data.first.begin(), m);
            data.second.insert(data.second.begin(), 0);
            startDiff -= minDiff;
            // Rcpp::Rcout << "pad start: " << m << "/" << startDiff << "\n";
        }
        auto endDiff = end - data.first.back();
        for (unsigned i=0; i<window && endDiff>0.0; ++i)
        {
            const auto m = std::min(data.first.back() + minDiff, end);
            data.first.push_back(m);
            data.second.push_back(0);
            endDiff -= minDiff;
            // Rcpp::Rcout << "pad end: " << m << "/" << endDiff << "\n";
        }
    }
    
    data.second = movingAverage(data.second, window);
}

template <typename T> std::pair<T, T>
summedFrameProps(const std::pair<std::vector<T>, std::vector<SpectrumRawTypes::Intensity>> &data,
                 T start, T end)
{
    T wm = 0.0, bp = 0.0;
    SpectrumRawTypes::Intensity totInten = 0, maxInten = 0;
    for (size_t i=0; i<data.first.size(); ++i)
    {
        // omit smoothing extension range
        // NOTE: assume data is sorted
        if (data.first[i] < start)
            continue;
        if (end != 0.0 && data.first[i] > end)
            break;
        
        // Rcpp::Rcout << "summedFrameProps: " << data.first[i] << "/" << data.second[i] << "\n";
        
        wm += (data.first[i] * data.second[i]);
        totInten += data.second[i];
        if (data.second[i] > maxInten)
        {
            maxInten = data.second[i];
            bp = data.first[i];
        }
    }
    
    if (totInten > 0)
        wm /= totInten;
    
    return std::make_pair(wm, bp);
}

template<typename T> Rcpp::List convertEICProfilesToR(const T &profiles, const char *name)
{
    Rcpp::List ret(profiles.size());
    for (size_t i=0; i<profiles.size(); ++i)
    {
        Rcpp::NumericMatrix mat(profiles[i].first.size(), 2);
        for (size_t j=0; j<profiles[i].first.size(); ++j)
        {
            mat(j, 0) = profiles[i].first[j];
            mat(j, 1) = profiles[i].second[j];
        }
        Rcpp::colnames(mat) = Rcpp::CharacterVector::create(name, "intensity");
        ret[i] = mat;
    }
    return ret;
}

}


void EIC::setSummedFrameMZ()
{
    if (empty() || updateSumMZIndex >= mzs.size())
        return;
    
    const auto updateTime = getUpdateSumMZTime();
    
#ifdef LOG_EIC
    Rcpp::Rcout << "EIC set MZ: " << updateSumMZIndex << " @ " << updateTime << std::endl;
#endif
    
    while (!mzSummer.empty() && (updateTime - mzSummer.frontTime()) > sumWindowMZ)
        mzSummer.pop();
    
    auto summedMZFrame = mzSummer.get();
#ifdef LOG_EIC
    Rcpp::Rcout << "  EIC mz frame: " << summedMZFrame.first.size() << "/" << summedMZFrame.second.size() << std::endl;
    for (size_t i=0; i<std::min(summedMZFrame.first.size(), size_t(5)); ++i)
        Rcpp::Rcout << "    " << summedMZFrame.first[i] << "/" << summedMZFrame.second[i] << std::endl;
#endif
    if (smoothWindowMZ > 0)
        smoothSummedFrameData(summedMZFrame, smoothWindowMZ, mzStart - smoothExtMZ, mzEnd + smoothExtMZ);
    const auto propsMZ = summedFrameProps(summedMZFrame, mzStart, mzEnd);
#ifdef LOG_EIC
    Rcpp::Rcout << "  EIC mz props: " << propsMZ.first << "/" << propsMZ.second << std::endl;
#endif
    mzs[updateSumMZIndex] = propsMZ.first; mzsBP[updateSumMZIndex] = propsMZ.second;
    if (saveMZProfiles)
        mzProfiles[updateSumMZIndex] = std::move(summedMZFrame);
    
    ++updateSumMZIndex;
}

void EIC::setSummedFrameMob()
{
    if (empty() || updateSumMobIndex >= mobilities.size())
        return;
    
    const auto updateTime = getUpdateSumMobTime();
    
#ifdef LOG_EIC
    Rcpp::Rcout << "EIC set Mob: " << updateSumMobIndex << " @ " << updateTime << std::endl;
#endif
    
    while (!mobSummer.empty() && (updateTime - mobSummer.frontTime()) > sumWindowMob)
        mobSummer.pop();
    
    auto summedMobFrame = mobSummer.get();
#ifdef LOG_EIC
    Rcpp::Rcout << "  EIC mob frame: " << summedMobFrame.first.size() << "/" << summedMobFrame.second.size() << std::endl;
    for (size_t i=0; i<std::min(summedMobFrame.first.size(), size_t(5)); ++i)
        Rcpp::Rcout << "    " << summedMobFrame.first[i] << "/" << summedMobFrame.second[i] << std::endl;
#endif
    if (smoothWindowMob > 0)
        smoothSummedFrameData(summedMobFrame, smoothWindowMob, mobStart - smoothExtMob, mobEnd + smoothExtMob);
    const auto propsMob = summedFrameProps(summedMobFrame, mobStart, mobEnd);
#ifdef LOG_EIC
    Rcpp::Rcout << "  EIC mob props: " << propsMob.first << "/" << propsMob.second << std::endl;
#endif
    mobilities[updateSumMobIndex] = propsMob.first; mobilitiesBP[updateSumMobIndex] = propsMob.second;
    if (saveEIMs)
        EIMs[updateSumMobIndex] = std::move(summedMobFrame);
    
    ++updateSumMobIndex;
}

void EIC::updateFrameSummer(SpectrumRawTypes::Time curTime)
{
    if (empty())
        return;

    while (!mzSummer.empty() && (curTime - getUpdateSumMZTime()) > sumWindowMZ)
        setSummedFrameMZ();
    mzSummer.add(curTime, std::move(curPoint.allMZs), std::move(curPoint.allIntsMZs));
    
    if (mode == EICMode::FULL)
    {
        while (!mobSummer.empty() && (curTime - getUpdateSumMobTime()) > sumWindowMob)
            setSummedFrameMob();
        mobSummer.add(curTime, std::move(curPoint.allMobs), std::move(curPoint.allIntsMobs));
    }
}

void EIC::commitPoints(SpectrumRawTypes::Scan curScanInd)
{
    scanInds.push_back(curScanInd);
    mzs.push_back(curPoint.mz);
    intensities.push_back(curPoint.intensity);
    if (mode == EICMode::FULL || mode == EICMode::FULL_MZ)
    {
        mzMins.push_back(curPoint.mzMin);
        mzMaxs.push_back(curPoint.mzMax);
        mzsBP.push_back(curPoint.mzBP);
    }
    if (withMob && saveMZProfiles && (mode == EICMode::FULL || mode == EICMode::FULL_MZ))
        mzProfiles.emplace_back(std::make_pair(std::vector<SpectrumRawTypes::Mass>(), std::vector<SpectrumRawTypes::Intensity>()));
    if (withMob && mode == EICMode::FULL)
    {
        mobMins.push_back(curPoint.mobMin);
        mobMaxs.push_back(curPoint.mobMax);
        // will be set by updateFrameSummer()
        mobilities.push_back(0);
        mobilitiesBP.push_back(0);
        if (saveEIMs)
            EIMs.emplace_back(std::make_pair(std::vector<SpectrumRawTypes::Mobility>(), std::vector<SpectrumRawTypes::Intensity>()));
    }
}

void EIC::addPoint(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten)
{
    if (inten == 0.0)
        return;
    
    curPoint.intensity += inten;
    
    if ((mode == EICMode::FULL || mode == EICMode::FULL_MZ))
    {
        if (curPoint.mzMin == 0.0 || mz < curPoint.mzMin)
            curPoint.mzMin = mz;
        if (mz > curPoint.mzMax)
            curPoint.mzMax = mz;
        if (!withMob)
        {
            curPoint.mz += mz * inten;
            if (inten > curPoint.intensityBP)
            {
                curPoint.intensityBP = inten;
                curPoint.mzBP = mz;
            }
        }
    }
}

void EIC::addPoint(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Mobility mob, SpectrumRawTypes::Intensity inten)
{
    if (inten == 0.0)
        return;
    
       // Rcpp::Rcout << "addPoint mob: " << mz << "/" << mob << "\t" << inten << "\t" << curPoint.mzMin << "/" << curPoint.mzMax << "\t" << std::endl;
    
    if (mode == EICMode::FULL || mode == EICMode::FULL_MZ)
    {
        // add data for frame summers
        // NOTE: point is outside mz/mob range then it is a smoothing extension point
        
        if (pointMobInRange(mob))
        {
            curPoint.allMZs.push_back(mz);
            curPoint.allIntsMZs.push_back(inten);
        }
        if (mode == EICMode::FULL && pointMZInRange(mz)) // NOTE: if not then in mz extension range
        {
            curPoint.allMobs.push_back(mob);
            curPoint.allIntsMobs.push_back(inten);
        }
        
        // Rcpp::Rcout << "  addPoint mob: " << mz << "/" << mob << "/" << inten  << std::endl;
    }
    
    if (!pointMZInRange(mz) || !pointMobInRange(mob))
        return; // out of range (but may be in extension range, as that should've been checked before calling addPoint())

    addPoint(mz, inten);
    
    if (mode == EICMode::FULL)
    {
        if (curPoint.mobMin == 0.0 || mob < curPoint.mobMin)
            curPoint.mobMin = mob;
        if (mob > curPoint.mobMax)
            curPoint.mobMax = mob;
    }
}

void EIC::commit(SpectrumRawTypes::Scan curScanInd, SpectrumRawTypes::Time curTime)
{
    const auto prevScanInd = (empty()) ? 0 : scanInds.back();
    
    if (curPoint.intensity > 0)
    {
        commitPoints(curScanInd);
        if (curPoint.intensity > maxIntensity)
            maxIntensity = curPoint.intensity;
        if (mode == EICMode::FULL || mode == EICMode::FULL_MZ)
        {
            if (withMob)
                updateFrameSummer(curTime);
            else
                curPoint.mz /= curPoint.intensity;
        }
    }
    
    if (minAdjIntensity > 0.0 && ((!enoughTimeAboveThr && minAdjTime > 0.0) ||
        (!enoughPointsAboveThr && minAdjPoints > 0)))
    {
        // if first point, gap in scans, or below intensity threshold, reset
        if (curPoint.intensity == 0.0 || scanInds.size() == 1 || prevScanInd != (curScanInd - 1) ||
            !numberGTE(curPoint.intensity, minAdjIntensity))
        {
            // Rcpp::Rcout << "Stop! EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
            startTimeAboveThr = 0.0;
            adjPointsAboveThr = 0;   
        }
        else
        {
            if (minAdjTime > 0.0)
            {
                if (startTimeAboveThr == 0.0)
                    startTimeAboveThr = curTime;
                else if (numberGTE(curTime - startTimeAboveThr, minAdjTime))
                    enoughTimeAboveThr = true;
            }
            if (minAdjPoints > 0)
            {
                ++adjPointsAboveThr;
                // Rcpp::Rcout << "above: EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
                if (adjPointsAboveThr >= minAdjPoints)
                    enoughPointsAboveThr = true;
            }
        }
    }
    
    curPoint.clear();
    // Rcpp::Rcout << "Cleared!\n";
}

void EIC::pad(SpectrumRawTypes::Scan scanStart, SpectrumRawTypes::Scan scanEnd)
{
    if (mode != EICMode::SIMPLE)
        Rcpp::stop("Padding only supported for simple EICs");
    
    if (empty())
    {
        scanInds = { scanStart, scanEnd };
        intensities.resize(2, 0.0);
        return;
    }

    std::vector<SpectrumRawTypes::Scan> paddedScans;
    std::vector<SpectrumRawTypes::Intensity> paddedInts;
    
    for (auto it=scanInds.begin(); it!=scanInds.end(); )
    {
        if (*it > 0)
        {
            // add zero point if previous was missing in EIC and within the RT range
            const auto prevScanInd = *it - 1;
            if ((paddedScans.empty() || prevScanInd != paddedScans.back()) && prevScanInd >= scanStart)
                
            {
                paddedScans.push_back(prevScanInd);
                paddedInts.push_back(0.0);
            }
        }
        
        // add EIC point
        paddedScans.push_back(*it);
        paddedInts.push_back(intensities[std::distance(scanInds.begin(), it)]);
        
        // add zero point if next is missing in EIC
        const auto nextScanInd = *it + 1;
        ++it;
        if ((it == scanInds.end() || *it != nextScanInd) && nextScanInd <= scanEnd)
        {
            paddedScans.push_back(nextScanInd);
            paddedInts.push_back(0.0);
        }
    }
    
    // add start/stop points if needed
    if (paddedScans.empty() || paddedScans.front() > scanStart)
    {
        paddedScans.insert(paddedScans.begin(), scanStart);
        paddedInts.insert(paddedInts.begin(), 0.0);
    }
    if (paddedScans.empty() || paddedScans.back() < scanEnd)
    {
        paddedScans.push_back(scanEnd);
        paddedInts.push_back(0.0);
    }
    
    scanInds = std::move(paddedScans);
    intensities = std::move(paddedInts);
}

void EIC::finalize()
{
    if (empty())
        return;

    // fill in left-over trailing scans
    while (!mzSummer.empty() && updateSumMZIndex < mzs.size())
        setSummedFrameMZ();
    while (!mobSummer.empty() && updateSumMobIndex < mobilities.size())
        setSummedFrameMob();
}        

void SimpleEIC::fillGaps(SpectrumRawTypes::Time medRTDiff, double gapFactor, bool pad,
                         SpectrumRawTypes::Time timeStart, SpectrumRawTypes::Time timeEnd)
{
    if (times.size() < 2)
        return;
    
    auto filled = fillRTGaps(times, medRTDiff, gapFactor);
    auto &filledTimes = filled.second;
    std::vector<SpectrumRawTypes::Intensity> filledInts(filled.first.size(), 0.0);
    // fill in intensities of existing time points, leave gaps at zero intensity
    for (size_t i=0; i<filled.first.size(); ++i)
    {
        if (filled.first[i] < times.size())
            filledInts[i] = intensities[filled.first[i]];
    }
    
    // check if padding is needed, eg in case front/back fell in gap
    if (pad)
    {
        // add two points in front/back: one just before/after the first/last and the other at the start/end
        // time range
        auto diff = filledTimes.front() - timeStart;
        if (diff > (medRTDiff * gapFactor))
        {
            filledTimes.insert(filledTimes.begin(), filledTimes.front() - medRTDiff);
            filledInts.insert(filledInts.begin(), 0.0);
            if (diff >= (medRTDiff * 2.0))
            {
                filledTimes.insert(filledTimes.begin(), timeStart);
                filledInts.insert(filledInts.begin(), 0.0);
            }
        }
        diff = timeEnd - filledTimes.back();
        if (diff > (medRTDiff * gapFactor))
        {
            filledTimes.push_back(filledTimes.back() + medRTDiff);
            filledInts.push_back(0.0);
            if (diff >= (medRTDiff * 2.0))
            {
                filledTimes.push_back(timeEnd);
                filledInts.push_back(0.0);
            }
        }
    }
    
    times = std::move(filledTimes);
    intensities = std::move(filledInts);
}

// [[Rcpp::export]]
Rcpp::List getEICList(const MSReadBackend &backend, const std::vector<SpectrumRawTypes::Mass> &startMZs,
                      const std::vector<SpectrumRawTypes::Mass> &endMZs,
                      const std::vector<SpectrumRawTypes::Time> &startTimes,
                      const std::vector<SpectrumRawTypes::Time> &endTimes,
                      const std::vector<SpectrumRawTypes::Mobility> &startMobs,
                      const std::vector<SpectrumRawTypes::Mobility> &endMobs,
                      SpectrumRawTypes::Time gapFactor, SpectrumRawTypes::Mass mzExpIMSWindow,
                      SpectrumRawTypes::Intensity minIntensityIMS, const std::string &mode = "simple",
                      SpectrumRawTypes::Time sumWindowMZ = 0, SpectrumRawTypes::Time sumWindowMob = 0,
                      unsigned smoothWindowMZ = 3, unsigned smoothWindowMob = 3, SpectrumRawTypes::Mass smoothExtMZ = 0,
                      SpectrumRawTypes::Mobility smoothExtMob = 0, bool saveMZProfiles = false, bool saveEIMs = false,
                      bool pad = false, SpectrumRawTypes::Intensity minEICIntensity = 0,
                      SpectrumRawTypes::Time minEICAdjTime = 0, unsigned minEICAdjPoints = 0,
                      SpectrumRawTypes::Intensity minEICAdjIntensity = 0, unsigned topMost = 0)
{
    // NOTE: startTimes/endTimes may be length one vectors, in which case they are used for all EICs
    
    struct AllPeaks
    {
        std::vector<SpectrumRawTypes::Scan> indices;
        std::vector<SpectrumRawTypes::Mass> mzs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        std::vector<SpectrumRawTypes::Mobility> mobilities;
        AllPeaks() = default;
        AllPeaks(size_t size, bool withMob) : indices(size), mzs(size), intensities(size), mobilities((withMob) ? size : 0) { }
        void clear(void)
        {
            indices.clear();
            mzs.clear();
            intensities.clear();
            mobilities.clear();
        }
    };
    
    const auto EICCount = startMZs.size();
    if (EICCount == 0)
        return Rcpp::List();
    
    const auto eicMode = (mode == "simple") ? EICMode::SIMPLE : (mode == "full") ? EICMode::FULL : (mode == "full_mz") ? EICMode::FULL_MZ : EICMode::TEST;
    
    const auto &specMeta = backend.getSpecMetadata();
    
    if (smoothWindowMZ == 0 || (eicMode != EICMode::FULL && eicMode != EICMode::FULL_MZ))
        smoothExtMZ = 0;
    if (smoothWindowMob == 0 || eicMode != EICMode::FULL)
        smoothExtMob = 0;
    
    const auto minMZ = (*(std::min_element(startMZs.begin(), startMZs.end()))) - smoothExtMZ;
    const auto maxMZ = (*(std::max_element(endMZs.begin(), endMZs.end()))) + smoothExtMZ;
    const auto minMob = (startMobs.empty()) ? 0 : (*(std::min_element(startMobs.begin(), startMobs.end())) - smoothExtMob);
    const auto maxMob = (endMobs.empty()) ? 0 : (*(std::max_element(endMobs.begin(), endMobs.end())) + smoothExtMob);
    
    const auto sfunc = [&](const SpectrumRaw &spec, const SpectrumRawSelection &ssel, size_t)
    {
        if (minMZ == 0.0 && maxMZ == 0.0 && minMob == 0.0 && maxMob == 0.0)
            return spec;
        
        const bool hasMob = spec.hasMobilities();
        
        // NOTE: specs are mz sorted
        const auto startIt = std::lower_bound(spec.getMZs().begin(), spec.getMZs().end(), minMZ);
        
        SpectrumRaw specf;
        for (size_t i=std::distance(spec.getMZs().begin(), startIt); i<spec.size(); ++i)
        {
            if (spec.getMZs()[i] > maxMZ)
                break;
            
            if (hasMob)
            {
                if (spec.getMobilities()[i] < minMob || (maxMob != 0.0 && spec.getMobilities()[i] > maxMob))
                    continue;
                specf.append(spec.getMZs()[i], spec.getIntensities()[i], spec.getMobilities()[i]);
            }
            else
                specf.append(spec.getMZs()[i], spec.getIntensities()[i]);
        }
        return specf;
    };
    
    SimpleTimer timer(true);
    
    // get scan selections: if EICCount is small, we minimize them to save file reads. For larger counts, we just read
    // the full min/max time range, as figuring out the specifics will take relatively long
    std::vector<std::vector<SpectrumRawSelection>> scanSels(1);
    
    if (EICCount < 1000)
    {
        timer.start("Collecting all scans (selectively)");
        std::set<SpectrumRawTypes::Scan> allScans;
        for (size_t i=0; i<startTimes.size(); ++i)
        {
            const auto sels = getSpecRawSelections(specMeta, makeNumRange(startTimes[i], endTimes[i]),
                                                   SpectrumRawTypes::MSLevel::MS1, 0);
            for (const auto &sel : sels)
                allScans.insert(sel.index);
        }
        for (const auto &scan : allScans)
            scanSels[0].emplace_back(scan);
        timer.stop();
    }
    else
    {
        timer.start("Collecting all scans (range)");
        scanSels[0] = getSpecRawSelections(specMeta, makeNumRange(*std::min_element(startTimes.begin(), startTimes.end()),
                                                                  *std::max_element(endTimes.begin(), endTimes.end())),
                                                                  SpectrumRawTypes::MSLevel::MS1, 0);
        timer.stop();
    }
    
    timer.start("Loading all MS spectra");
    auto allSpectra = applyMSData<SpectrumRaw>(backend, SpectrumRawTypes::MSLevel::MS1, scanSels, sfunc,
                                               minIntensityIMS, SpectrumRawTypes::MSSortType::MZ);
    timer.stop();
    if (allSpectra.empty())
        return Rcpp::List();
    
    timer.start("Preparing EICs (AllPeaks)");
    const auto totalPeaks = std::accumulate(allSpectra[0].begin(), allSpectra[0].end(), 0,
                                            [](size_t sum, const SpectrumRaw &spec) { return sum + spec.size(); });
    
    bool anySpecHasMob = false; // NOTE: assume all MS spectra have or have not IMS data (UNDONE?)
    for (size_t i=0; i<allSpectra[0].size(); ++i)
    {
        if (allSpectra[0][i].getMZs().empty())
            continue; // no peaks in this spectrum, need some to see if there are mobilities
        anySpecHasMob = allSpectra[0][i].hasMobilities();
        break;
    }
    
    if (!anySpecHasMob)
    {
        smoothExtMZ = 0;
        smoothExtMob = 0;
    }
    
    AllPeaks allPeaks(totalPeaks, anySpecHasMob);
    for (size_t i=0, peaksi=0; i<allSpectra[0].size(); ++i)
    {
        if (allSpectra[0][i].getMZs().empty())
            continue; // no peaks in this spectrum
        
        std::fill(allPeaks.indices.begin() + peaksi, allPeaks.indices.begin() + peaksi + allSpectra[0][i].size(),
                  scanSels[0][i].index);
        std::copy(allSpectra[0][i].getMZs().begin(), allSpectra[0][i].getMZs().end(),
                  allPeaks.mzs.begin() + peaksi);
        std::copy(allSpectra[0][i].getIntensities().begin(), allSpectra[0][i].getIntensities().end(),
                  allPeaks.intensities.begin() + peaksi);
        if (anySpecHasMob)
        {
            std::copy(allSpectra[0][i].getMobilities().begin(), allSpectra[0][i].getMobilities().end(),
                      allPeaks.mobilities.begin() + peaksi);
        }
        
        peaksi += allSpectra[0][i].size();
        allSpectra[0][i].clear();
    }
    timer.stop();
    
    timer.start("Sorting all peaks");
    const auto sortedInds = getSortedInds(allPeaks.mzs);
    
    AllPeaks allPeaksSorted(allPeaks.indices.size(), anySpecHasMob);
    for (size_t i=0; i<sortedInds.size(); ++i)
    {
        const auto j = sortedInds[i];
        allPeaksSorted.indices[i] = allPeaks.indices[j];
        allPeaksSorted.mzs[i] = allPeaks.mzs[j];
        allPeaksSorted.intensities[i] = allPeaks.intensities[j];
        if (anySpecHasMob)
            allPeaksSorted.mobilities[i] = allPeaks.mobilities[j];
    }
    allPeaks.clear();
    timer.stop();
    
    timer.start("Processing EICs");
    std::vector<EIC> allEICs;
    std::vector<SpectrumRawTypes::Intensity> allEICMaxIntensities(EICCount);
    
    allEICs.reserve(EICCount);
    for (size_t i=0; i<EICCount; ++i)
        allEICs.emplace_back(specMeta.first.times, eicMode, anySpecHasMob, minEICAdjIntensity, minEICAdjTime,
                             minEICAdjPoints, sumWindowMZ, sumWindowMob, smoothWindowMZ, smoothWindowMob, smoothExtMZ,
                             smoothExtMob, saveMZProfiles, saveEIMs);
    
    #pragma omp parallel for
    for (size_t i=0; i<EICCount; ++i)
    {
        EIC &eic = allEICs[i];
        
        const auto mzStart = startMZs[i] - ((anySpecHasMob) ? mzExpIMSWindow : 0.0);
        const auto mzEnd = endMZs[i] + ((anySpecHasMob) ? mzExpIMSWindow : 0.0);
        const auto mzExtStart = (mzStart > 0.0) ? (mzStart - smoothExtMZ) : 0.0;
        const auto mzExtEnd = (mzEnd > 0.0) ? mzEnd + smoothExtMZ : 0.0;
        const auto timeStart = (startTimes.size() == 1) ? startTimes[0] : startTimes[i];
        const auto timeEnd = (endTimes.size() == 1) ? endTimes[0] : endTimes[i];
        const auto mobStart = (anySpecHasMob) ? startMobs[i] : 0.0, mobEnd = (anySpecHasMob) ? endMobs[i] : 0.0;
        const auto mobExtStart = (mobStart > 0.0) ? (mobStart - smoothExtMob) : 0.0;
        const auto mobExtEnd = (mobEnd > 0.0) ? (mobEnd + smoothExtMob) : 0.0;
        
        const auto itStart = std::lower_bound(allPeaksSorted.mzs.cbegin(), allPeaksSorted.mzs.cend(), mzStart);
        if (itStart == allPeaksSorted.mzs.cend() || *itStart > mzEnd)
            continue;
        const auto itEnd = std::prev(std::upper_bound(itStart, allPeaksSorted.mzs.cend(), mzEnd));
        const auto startInd = std::distance(allPeaksSorted.mzs.cbegin(), itStart);
        auto endInd = std::distance(allPeaksSorted.mzs.cbegin(), itEnd);
        if (startInd > endInd)
            endInd = startInd; // only one peak at the end
        const auto sortedInds = getSortedInds(allPeaksSorted.indices.cbegin() + startInd,
                                              allPeaksSorted.indices.cbegin() + endInd);
        
        eic.setBoundaries(mzStart, mzEnd, mobStart, mobEnd);
        
        size_t curScanInd = allPeaksSorted.indices.size();
        bool init = true;
        for (size_t j=0; ; ++j)
        {
            const bool ended = j == sortedInds.size();
            if (ended && init)
                break; // nothing was done
            
            const auto allPeaksSortedInd = (ended) ? 0 : (startInd + sortedInds[j]);
            const auto scanInd = allPeaksSorted.indices[allPeaksSortedInd];
            
            //Rcpp::Rcout << "EIC: " << j << "/" << sortedInds.size() << "/" << specMeta.first.scans[curScanInd] << "/" << curScanInd << "/" << init << "/" << ended << "/" << startInd << "/" << endInd << "\n";
            
            if (ended || (!init && scanInd != curScanInd))
            {
                const auto curTime = specMeta.first.times[curScanInd];
                eic.commit(curScanInd, curTime);
                if (ended)
                    break;
                curScanInd = scanInd;
            }
            
            const auto time = specMeta.first.times[scanInd];
            if (!timeGTE(time, timeStart) || (timeEnd != 0.0 && !timeLTE(time, timeEnd)))
                continue;
            
            //Rcpp::Rcout << "EIC: " << allPeaksSortedInd << "/" << j << "/" << specMeta.first.times[allPeaksSorted.indices[allPeaksSortedInd]] << "/" << mz << "/" << mzStart << "/" << mzEnd << "\n";
            const auto mz = allPeaksSorted.mzs[allPeaksSortedInd];
            const auto mob = (anySpecHasMob) ? allPeaksSorted.mobilities[allPeaksSortedInd] : 0;
            const auto inten = allPeaksSorted.intensities[allPeaksSortedInd];
            
            // skip check if outside smooth window.
            // NOTE: xxxExtStart/xxxExtEnd == xxxStart/xxxEnd if no smoothing is performed
            if (mz < mzExtStart || (mzExtEnd != 0.0 && mz > mzExtEnd) ||
                (anySpecHasMob && (mob < mobExtStart || (mobExtEnd != 0.0 && mob > mobExtEnd))))
                continue;
            
            if (init)
            {
                curScanInd = scanInd;
                init = false;
            }
            
            if (inten > 0)
            {
                if (anySpecHasMob)
                    eic.addPoint(mz, mob, inten);
                else
                    eic.addPoint(mz, inten);
            }
        }
        
        if (!eic.verifyAdjancency() || (minEICIntensity != 0.0 && eic.getMaxIntensity() < minEICIntensity))
        {
            eic.clear();
            continue;   
        }

        if (eicMode == EICMode::SIMPLE && pad)
        {
            // figure out start/stop scans
            auto startIt = std::lower_bound(specMeta.first.times.begin(), specMeta.first.times.end(), timeStart);
            if (startIt == specMeta.first.times.end())
                startIt = specMeta.first.times.begin();
            auto endIt = (timeEnd == 0.0) ? specMeta.first.times.end() : std::lower_bound(startIt, specMeta.first.times.end(), timeEnd);
            if (endIt == specMeta.first.times.end() || (endIt != specMeta.first.times.begin() && *endIt > timeEnd))
                --endIt;
            eic.pad(std::distance(specMeta.first.times.begin(), startIt),
                    std::distance(specMeta.first.times.begin(), endIt));
        }
        
        eic.finalize();
        
        allEICMaxIntensities[i] = eic.getMaxIntensity();
    }
    
    allPeaksSorted.clear(); // no need for this anymore
    timer.stop();
    
    SpectrumRawTypes::Intensity minMaxIntens = 0.0;
    if (topMost > 0 && allEICMaxIntensities.size() > topMost)
    {
        timer.start("Determining EIC topMost...");
        auto sortedAEMI = allEICMaxIntensities;
        std::nth_element(sortedAEMI.begin(), sortedAEMI.begin() + (topMost-1), sortedAEMI.end(),
                         std::greater<SpectrumRawTypes::Intensity>());
        minMaxIntens = sortedAEMI[topMost-1];
        timer.stop();
    }
    
    timer.start("Converting EICs to R data");
    const auto medRTDiff = medianRTDiff(specMeta.first.times);
    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        auto &eic = allEICs[i];
        
        if (allEICMaxIntensities[i] < minMaxIntens)
        {
            // just clear the EIC so an empty will be added below. Note that we would clear the EIC at the end anyway.
            eic.clear();
        }
        
        // NOTE: constructing the output list directly with create() is much faster than assigning items afterwards.
        // Hence, we do quite some if statements here...
        
        if (eicMode == EICMode::SIMPLE)
        {
            SimpleEIC simpleEIC(eic, specMeta.first.times);
            if (gapFactor > 0.0)
            {
                const auto timeStart = (startTimes.size() == 1) ? startTimes[0] : startTimes[i];
                const auto timeEnd = (endTimes.size() == 1) ? endTimes[0] : endTimes[i];
                simpleEIC.fillGaps(medRTDiff, gapFactor, pad, timeStart,
                                   (timeEnd == 0.0) ? specMeta.first.times.back() : timeEnd);
            }
            
            auto mat = Rcpp::NumericMatrix(simpleEIC.size(), 2);
            for (size_t i=0; i<simpleEIC.size(); ++i)
            {
                mat(i, 0) = simpleEIC.getTimes()[i];
                mat(i, 1) = simpleEIC.getIntensities()[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity");
            ret[i] = mat;
        }
        else if (eicMode == EICMode::FULL && anySpecHasMob)
        {
            auto mat = Rcpp::NumericMatrix(eic.size(), 10);
            for (size_t i=0; i<eic.size(); ++i)
            {
                mat(i, 0) = specMeta.first.times[eic.getScanIndices()[i]];
                mat(i, 1) = eic.getIntensities()[i];
                mat(i, 2) = eic.getMZs()[i];
                mat(i, 3) = eic.getMZsBP()[i];
                mat(i, 4) = eic.getMZMins()[i];
                mat(i, 5) = eic.getMZMaxs()[i];
                mat(i, 6) = eic.getMobilities()[i];
                mat(i, 7) = eic.getMobMins()[i];
                mat(i, 8) = eic.getMobMaxs()[i];
                mat(i, 9) = eic.getMobilitiesBP()[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity", "mz", "mzBP", "mzmin", "mzmax",
                           "mobility", "mobmin", "mobmax", "mobilityBP");
            ret[i] = mat;
        }
        else if (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ)
        {
            auto mat = Rcpp::NumericMatrix(eic.size(), 6);
            for (size_t i=0; i<eic.size(); ++i)
            {
                mat(i, 0) = specMeta.first.times[eic.getScanIndices()[i]];
                mat(i, 1) = eic.getIntensities()[i];
                mat(i, 2) = eic.getMZs()[i];
                mat(i, 3) = eic.getMZsBP()[i];
                mat(i, 4) = eic.getMZMins()[i];
                mat(i, 5) = eic.getMZMaxs()[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity", "mz", "mzBP", "mzmin", "mzmax");
            ret[i] = mat;
        }
        else // if (eicMode == EICMode::TEST)
        {
            ret[i] = !eic.empty();
        }
        
        eic.clear(saveMZProfiles, saveEIMs); // free memory as EICs may consume a lot
    }
    
    if (gapFactor > 0.0 && eicMode != EICMode::TEST)
    {
        // For padding we just have to add additional time points, EIC decompression will assume these are zero
        // intensity points.
        ret.attr("allXValues") = (fillRTGaps(specMeta.first.times, medRTDiff, gapFactor)).second;
    }
    
    if (anySpecHasMob && (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ) && saveMZProfiles)
    {
        Rcpp::List mzProfiles(EICCount);
        for (size_t i=0; i<EICCount; ++i)
        {
            const auto &eic = allEICs[i];
            if (eic.getMZProfiles().empty())
                continue;
            mzProfiles[i] = convertEICProfilesToR(eic.getMZProfiles(), "mz");
        }
        ret.attr("mzProfiles") = mzProfiles;
    }
    
    if (anySpecHasMob && eicMode == EICMode::FULL && saveEIMs)
    {
        Rcpp::List allEIMs(EICCount);
        for (size_t i=0; i<EICCount; ++i)
        {
            const auto &eic = allEICs[i];
            if (eic.getEIMs().empty())
                continue;
            allEIMs[i] = convertEICProfilesToR(eic.getEIMs(), "mobility");
        }
        ret.attr("EIMs") = allEIMs;
    }
    
    timer.stop();
    
    return ret;
}

// [[Rcpp::export]]
std::vector<SpectrumRawTypes::Intensity> doFillEIXIntensities(const std::vector<SpectrumRawTypes::Time> &allXValues,
                                                              const std::vector<SpectrumRawTypes::Time> &xvalues,
                                                              const std::vector<SpectrumRawTypes::Intensity> &intensities)
{
    return fillEIXIntensities(allXValues, xvalues, intensities);
}
