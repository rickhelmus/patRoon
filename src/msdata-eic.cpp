#include <Rcpp.h>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <map>

#include "msdata.hpp"
#include "msdata-eic.h"

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

template <typename T>
void smoothSummedFrameData(std::pair<std::vector<T>, std::vector<SpectrumRawTypes::Intensity>> &data,
                           unsigned window, T start, T end)
{
    if (data.first.size() < 3 || window < 1)
        return;
    
    // pad to get smoothing right
    
    T minDiff = 0;
    for (size_t j=1; j<data.first.size(); ++j)
    {
        const auto diff = data.first[j] - data.first[j-1];
        if (minDiff == 0 || diff < minDiff)
            minDiff = diff;
    }
    if (minDiff > 0.0)
    {
        auto startDiff = data.first.front() - start;
        for (unsigned j=0; j<window && startDiff>0.0; ++j)
        {
            const auto m = data.first.front() - minDiff;
            data.first.insert(data.first.begin(), m);
            data.second.insert(data.second.begin(), 0);
            startDiff -= minDiff;
        }
        auto endDiff = end - data.first.back();
        for (unsigned j=0; j<window && endDiff>0.0; ++j)
        {
            const auto m = data.first.back() + minDiff;
            data.first.push_back(m);
            data.second.push_back(0);
            endDiff -= minDiff;
        }
    }
    
    data.second = movingAverage(data.second, window);

    return std::make_pair(std::move(data.first), std::move(data.second));
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

std::pair<std::vector<SpectrumRawTypes::Time>, std::vector<SpectrumRawTypes::Intensity>>
    padEIC(const std::vector<SpectrumRawTypes::Time> &allTimes,
           SpectrumRawTypes::Time startTime, SpectrumRawTypes::Time endTime,
           const std::vector<SpectrumRawTypes::Time> &times,
           const std::vector<SpectrumRawTypes::Intensity> &intensities)
{
    std::vector<SpectrumRawTypes::Time> outTimes;
    std::vector<SpectrumRawTypes::Intensity> outIntensities;
    
    if (allTimes.empty())
        return std::make_pair(outTimes, outIntensities);
    if (times.empty())
    {
        outTimes = { allTimes.front(), allTimes.back() };
        outIntensities.resize(2, 0.0);
        return std::make_pair(outTimes, outIntensities);
    }
    
    // snap start/end ranges to allTimes
    auto startTimeIt = std::lower_bound(allTimes.begin(), allTimes.end(), startTime);
    if (startTimeIt == allTimes.end())
        startTimeIt = std::prev(allTimes.end());
    startTime = *startTimeIt;
    auto endTimeIt = std::lower_bound(startTimeIt, allTimes.end(), endTime);
    if (endTimeIt == allTimes.end() || (endTimeIt != allTimes.begin() && *endTimeIt > endTime))
        --endTimeIt;
    endTime = *endTimeIt;
    
    auto it = times.begin();
    auto allIt = allTimes.begin();
    
    while (it != times.end() && allIt != allTimes.end())
    {
        // corresponding iterator in allTimes
        const auto allItCorresp = std::lower_bound(allIt, allTimes.end(), *it);
        
        // add zero point if previous was missing in EIC and within the RT range
        if (allItCorresp != allTimes.begin())
        {
            const auto aicPrv = std::prev(allItCorresp);
            if ((it == times.begin() || (!compareTol(*(std::prev(it)), *aicPrv) && !compareTol(*aicPrv, outTimes.back()))) &&
                numberGTE(*aicPrv, startTime))
                
            {
                outTimes.push_back(*aicPrv);
                outIntensities.push_back(0.0);
            }
        }
        
        // add EIC point
        outTimes.push_back(*it);
        outIntensities.push_back(intensities[std::distance(times.begin(), it)]);
        
        allIt = std::next(allItCorresp);
        if (allIt == allTimes.end())
            break;
        
        // add zero point if next is missing in EIC
        ++it;
        if ((it == times.end() || !compareTol(*it, *allIt)) && numberLTE(*allIt, endTime))
        {
            outTimes.push_back(*allIt);
            outIntensities.push_back(0.0);
            ++allIt;
        }
    }
    
    // add final point if needed
    if (allIt != allTimes.end() && numberLTE(allTimes.back(), endTime))
    {
        outTimes.push_back(allTimes.back());
        outIntensities.push_back(0.0);
    }
    
    // add start/stop points if needed
    if (outTimes.empty() || outTimes.front() > startTime)
    {
        outTimes.insert(outTimes.begin(), startTime);
        outIntensities.insert(outIntensities.begin(), 0.0);
    }
    if (outTimes.empty() || outTimes.back() < endTime)
    {
        outTimes.push_back(endTime);
        outIntensities.push_back(0.0);
    }
    
    return std::make_pair(outTimes, outIntensities);
}

}

void EIC::setSummedFrame(SpectrumRawTypes::Scan scanInd)
{
    size_t ind = scanInds.size() - 1;
    for (; ; --ind)
    {
        if (scanInds[ind] == scanInd)
            break; // found it
        if (scanInds[ind] < scanInd || ind == 0)
            return; // not actually present
    }
    
    auto summedMZFrame = frameSummer.getMZs();
    if (smoothWindowMZ > 0)
        smoothSummedFrameData(summedMZFrame, smoothWindowMZ, mzStart - smoothExtMZ, mzEnd + smoothExtMZ);
    const auto propsMZ = summedFrameProps(summedMZFrame, mzStart, mzEnd);
    mzs[ind] = propsMZ.first; mzsBP[ind] = propsMZ.second;
    if (saveMZProfiles)
        MZProfiles[ind] = std::move(summedMZFrame);
    
    if (mode == EICMode::FULL)
    {
        auto summedMobFrame = frameSummer.getMobilities();
        if (smoothWindowMob > 0)
            smoothSummedFrameData(summedMobFrame, smoothWindowMob, mobStart - smoothExtMob, mobEnd + smoothExtMob);
        const auto propsMob = summedFrameProps(summedMobFrame, mobStart, mobEnd);
        mobilities[ind] = propsMob.first; mobilitiesBP[ind] = propsMob.second;
        if (saveEIMs)
            EIMs[ind] = std::move(summedMobFrame);
    }
}

void EIC::updateFrameSummer(SpectrumRawTypes::Scan curScanInd)
{
    const auto prvScanInd = (empty()) ? 0 : scanInds.back();
    SpectrumRawTypes::Scan sc = (empty()) ? curScanInd : (prvScanInd + 1);
    
#ifdef LOG_EIM
    Rcpp::Rcout << "EIC sum: " << i << "/" << j << "/" << sc << "/" << scanInd
                << "/" << frameSummer.size() << "/" <<  specMeta.first.times[sc] << "/" << curTime << std::endl;
#endif
    
    while (true)
    {
        if (sc != curScanInd)
        {
            // add zero-intensity points for missing scans
            frameSummer.addZero();
#ifdef LOG_EIM
            Rcpp::Rcout << "added zero:" << i << "/" << j << "/" << sc << "/" << curScanInd << std::endl;
#endif
        }
        else
        {
            if (mode == EICMode::FULL)
                frameSummer.add(std::move(curPoint.allMZs), std::move(curPoint.allMobs),
                                std::move(curPoint.allInts));
            else
                frameSummer.add(std::move(curPoint.allMZs), std::move(curPoint.allInts));
#ifdef LOG_EIM
            Rcpp::Rcout << "added points:" << i << "/" << j << "/" << curPointMZs.size() << "/" << curScanInd << std::endl;
#endif
        }
        
        // sumEIMs == 3 --> flank == 1 --> set position 1 back from current
        // sumEIMs == 5 --> flank == 2 --> set position 2 back from current
        // etc
        
        if (!empty() && sc >= frameSummer.flank())
            setSummedFrame(sc - frameSummer.flank());
#ifdef LOG_EIM
        Rcpp::Rcout << "EIC update: " << specMeta.first.times[sc] << "/" << j << "/" << sc << "/" << curScanInd
                    << "/" << prvScanInd << "/" << frameSummer.size() << "/" << frameSummer.sizeNoZero() << std::endl;
#endif
        
        if (sc == curScanInd)
            break;
        
        if (sc >= (prvScanInd + frameSummer.maxSize()))
            sc = curScanInd; // jump to current: doesn't make sense to add more zeroes
        else
            ++sc;
    }
}

void EIC::commitPoints(SpectrumRawTypes::Scan curScanInd)
{
    scanInds.push_back(curScanInd);
    mzs.push_back(curPoint.mz);
    mzMins.push_back(curPoint.mzMin);
    mzMaxs.push_back(curPoint.mzMax);
    intensities.push_back(curPoint.intensity);
    mzsBP.push_back(curPoint.mzBP);
    intensitiesBP.push_back(curPoint.intensityBP);
    if (withMob)
    {
        mobMins.push_back(curPoint.mobMin);
        mobMaxs.push_back(curPoint.mobMax);
        // will be set by updateSummedEIMs()
        mobilities.push_back(0);
        mobilitiesBP.push_back(0);
        if (saveMZProfiles)
            MZProfiles.emplace_back(std::make_pair(std::vector<SpectrumRawTypes::Mass>(), std::vector<SpectrumRawTypes::Intensity>()));
        if (saveEIMs)
            EIMs.emplace_back(std::make_pair(std::vector<SpectrumRawTypes::Mobility>(), std::vector<SpectrumRawTypes::Intensity>()));
    }
}

void EIC::addPoint(SpectrumRawTypes::Intensity inten)
{
    curPoint.intensity += inten;
}

void EIC::addPoint(SpectrumRawTypes::Mass mz, SpectrumRawTypes::Intensity inten)
{
    addPoint(inten);
    if (inten > 0 && (mode == EICMode::FULL || mode == EICMode::FULL_MZ))
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
    // UNDONE: need (to check) EICMode?
    if (mode == EICMode::FULL || mode == EICMode::FULL_MZ)
    {
        // add data for frame summer
        curPoint.allMZs.push_back(mz);
        if (mode == EICMode::FULL)
            curPoint.allMobs.push_back(mob);
        curPoint.allInts.push_back(inten);
    }
    
    addPoint(mz, inten);
    
    if (inten > 0 && mode == EICMode::FULL)
    {
        if (curPoint.mobMin == 0.0 || mob < curPoint.mobMin)
            curPoint.mobMin = mob;
        if (mob > curPoint.mobMax)
            curPoint.mobMax = mob;
    }
    // Rcpp::Rcout << "EIC: " << curPoint.time << "/" << j << "\t" << std::fixed << mz << "\t" << inten << "\t" << curPoint.mzMin << "/" << curPoint.mzMax << "\t" << mob << std::endl;
}

void EIC::commit(SpectrumRawTypes::Scan curScanInd, SpectrumRawTypes::Time curTime)
{
    const auto prevScanInd = (empty()) ? 0 : scanInds.back();
    
    if (curPoint.intensity > 0)
    {
        if (mode == EICMode::FULL || mode == EICMode::FULL_MZ)
        {
            if (withMob)
                updateFrameSummer(curScanInd);
            else
                curPoint.mz /= curPoint.intensity;
        }
        commitPoints(curScanInd);
    }
    
    if ((!enoughTimeAboveThr && minEICAdjTime > 0.0) || (!enoughPointsAboveThr && minEICAdjPoints > 0))
    {
        // if first point, gap in scans, or below intensity threshold, reset
        if (scanInds.size() == 1 || prevScanInd != (curScanInd - 1) ||
            !numberGTE(curPoint.intensity, minEICAdjIntensity))
        {
            // Rcpp::Rcout << "Stop! EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
            startTimeAboveThr = 0.0;
            adjPointsAboveThr = 0;   
        }
        else
        {
            if (minEICAdjTime > 0.0)
            {
                if (startTimeAboveThr == 0.0)
                    startTimeAboveThr = curTime;
                else if (numberGTE(curTime - startTimeAboveThr, minEICAdjTime))
                    enoughTimeAboveThr = true;
            }
            if (minEICAdjPoints > 0)
            {
                ++adjPointsAboveThr;
                // Rcpp::Rcout << "above: EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
                if (adjPointsAboveThr >= minEICAdjPoints)
                    enoughPointsAboveThr = true;
            }
        }
    }
    
    curPoint.clear();
}

void EIC::pad(SpectrumRawTypes::Scan scanStart, SpectrumRawTypes::Scan scanEnd)
{
    if (mode != EICMode::SIMPLE)
        Rcpp::stop("Padding only supported for simple EICs");
    
    if (empty())
    {
        scanInds = { scanStart, scanEnd };
        intensities.resize(2, 0.0);
    }
    
    std::vector<SpectrumRawTypes::Scan> paddedScans;
    std::vector<SpectrumRawTypes::Intensity> paddedInts;
    
    auto it = scanInds.begin();
    
    while (it != scanInds.end())
    {
        // add zero point if previous was missing in EIC and within the RT range
        const auto prevScanInd = *it - 1;
        if ((paddedScans.empty() || prevScanInd != paddedScans.back()) && prevScanInd >= scanStart)
            
        {
            paddedScans.push_back(prevScanInd);
            paddedInts.push_back(0.0);
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
    if (!empty())
    {
        // fill in left-over trailing scans
        // NOTE: at this point, prvScanInd refers to the last filled in scan
        const auto lastScanInd = scanInds.back();
        for (SpectrumRawTypes::Scan s=lastScanInd+1-std::min(frameSummer.flank(), frameSummer.size());
             !frameSummer.empty() && s<=lastScanInd; ++s)
        {
#ifdef LOG_EIM
            Rcpp::Rcout << "EIM purge: " << s << "/" << curScanInd << "/" << specMeta.first.times[s] << "/"
                        << frameSummer.size() << std::endl;
#endif
            setSummedFrame(s);
            frameSummer.pop();
        }
    }
    
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
                      unsigned sumEIMs = 1, unsigned smoothWindowMZ = 3, unsigned smoothWindowMob = 3,
                      SpectrumRawTypes::Mass smoothExtMZ = 0, SpectrumRawTypes::Mobility smoothExtMob = 0,
                      bool saveEIMs = false, bool pad = false, SpectrumRawTypes::Intensity minEICIntensity = 0,
                      SpectrumRawTypes::Time minEICAdjTime = 0, unsigned minEICAdjPoints = 0,
                      SpectrumRawTypes::Intensity minEICAdjIntensity = 0, unsigned topMost = 0)
{
    // NOTE: startTimes/endTimes may be length one vectors, in which case they are used for all EICs
    
    // #define LOG_EIM
    
    struct EICPoint
    {
        SpectrumRawTypes::Time time = 0.0;
        SpectrumRawTypes::Mass mz = 0.0, mzMin = 0.0, mzMax = 0.0;
        SpectrumRawTypes::Intensity intensity = 0;
        SpectrumRawTypes::Mobility mobMin = 0.0, mobMax = 0.0;
        SpectrumRawTypes::Mass mzBP = 0.0;
        SpectrumRawTypes::Intensity intensityBP = 0.0;
    };
    
    struct EIC
    {
        std::vector<SpectrumRawTypes::Time> times;
        std::vector<SpectrumRawTypes::Mass> mzs, mzMins, mzMaxs;
        std::vector<SpectrumRawTypes::Intensity> intensities;
        std::vector<SpectrumRawTypes::Mobility> mobilities, mobMins, mobMaxs;
        std::vector<SpectrumRawTypes::Mass> mzsBP;
        std::vector<SpectrumRawTypes::Mobility> mobilitiesBP;
        std::vector<SpectrumRawTypes::Intensity> intensitiesBP;
        std::vector<std::pair<std::vector<SpectrumRawTypes::Mobility>, std::vector<SpectrumRawTypes::Intensity>>> EIMs; // for summed EIMs
        void addPoint(const EICPoint &p, bool withMob, bool saveEIM)
        {
            times.push_back(p.time);
            mzs.push_back(p.mz);
            mzMins.push_back(p.mzMin);
            mzMaxs.push_back(p.mzMax);
            intensities.push_back(p.intensity);
            mzsBP.push_back(p.mzBP);
            intensitiesBP.push_back(p.intensityBP);
            if (withMob)
            {
                mobMins.push_back(p.mobMin);
                mobMaxs.push_back(p.mobMax);
                // will be set by updateSummedEIMs()
                mobilities.push_back(0);
                mobilitiesBP.push_back(0);
                if (saveEIM)
                    EIMs.emplace_back(std::make_pair(std::vector<SpectrumRawTypes::Mobility>(), std::vector<SpectrumRawTypes::Intensity>()));
            }
        }
        void addPoint(SpectrumRawTypes::Time t, SpectrumRawTypes::Intensity i) { times.push_back(t); intensities.push_back(i); }
        void clear(void)
        {
            times.clear();
            mzs.clear();
            mzMins.clear();
            mzMaxs.clear();
            intensities.clear();
            mobMins.clear();
            mobMaxs.clear();
            mzsBP.clear();
            intensitiesBP.clear();
        }
        size_t size(void) const { return times.size(); }
        bool empty(void) const { return times.empty(); }
    };
    
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
    
    const unsigned EIMFlank = (sumEIMs > 1) ? (sumEIMs - 1) / 2 : 0;
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
    
    const auto updateSummedEIMs = [saveEIMs, smoothWindowMZ, smoothWindowMob, smoothExtMZ, smoothExtMob, &specMeta]
        (SpectrumRawTypes::Scan scanInd, const IMSFrameSummer &frameSummer, EIC &eic,
         SpectrumRawTypes::Mass mzStart, SpectrumRawTypes::Mass mzEnd,
         SpectrumRawTypes::Mobility mobStart, SpectrumRawTypes::Mobility mobEnd)
    {
        if (eic.empty())
            return;
        
        const auto scanTime = specMeta.first.times[scanInd];
        
        for (size_t i=eic.size()-1; ; --i)
        {
#ifdef LOG_EIM
            Rcpp::Rcout << " compare to: " << i << "/" << scanInd << "/" << scanTime << "/" << eic.times[i] << std::endl;
#endif
            if (scanTime > eic.times[i])
                break; // point not actually added
            if (compareTol(scanTime, eic.times[i]))
            {
                const auto mobsAndInts = frameSummer.getMobilities(smoothWindowMob, mobStart - smoothExtMob, mobEnd + smoothExtMob);
                SpectrumRawTypes::Intensity totInten = 0, maxInten = 0;
                for (size_t j=0; j<mobsAndInts.first.size(); ++j)
                {
                    // omit smoothing extension mobility range
                    // NOTE: mobilities are sorted
                    if (mobsAndInts.first[j] < mobStart)
                        continue;
                    if (mobEnd != 0.0 && mobsAndInts.first[j] > mobEnd)
                        break;
                    
                    eic.mobilities[i] += (mobsAndInts.first[j] * mobsAndInts.second[j]);
                    totInten += mobsAndInts.second[j];
                    if (mobsAndInts.second[j] > maxInten)
                    {
                        maxInten = mobsAndInts.second[j];
                        eic.mobilitiesBP[i] = mobsAndInts.first[j];
                    }
                }
                
                if (totInten > 0)
                    eic.mobilities[i] /= totInten;
                
                if (saveEIMs)
                    eic.EIMs[i] = std::make_pair(mobsAndInts.first, mobsAndInts.second);
#ifdef LOG_EIM
                Rcpp::Rcout << "  set: " << totInten << "/" << maxInten << "/" << eic.mobilities[i] << "/" << scanTime << std::endl;
#endif
                break;
            }
            if (i == 0)
                break;
        }
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
    std::vector<EIC> allEICs(EICCount);
    std::vector<SpectrumRawTypes::Intensity> allEICMaxIntensities(EICCount);
    
    #pragma omp parallel for
    for (size_t i=0; i<EICCount; ++i)
    {
        EIC eic;
        SpectrumRawTypes::Intensity maxInten = 0;
        SpectrumRawTypes::Time startTimeAboveThr = 0;
        unsigned adjPointsAboveThr = 0;
        bool enoughTimeAboveThr = false, enoughPointsAboveThr = false;
        
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
        if (itStart == allPeaksSorted.mzs.cend())
            continue;
        const auto itEnd = std::prev(std::upper_bound(itStart, allPeaksSorted.mzs.cend(), mzEnd));
        const auto startInd = std::distance(allPeaksSorted.mzs.cbegin(), itStart);
        auto endInd = std::distance(allPeaksSorted.mzs.cbegin(), itEnd);
        if (startInd > endInd)
            endInd = startInd; // only one peak at the end
        const auto sortedInds = getSortedInds(allPeaksSorted.indices.cbegin() + startInd,
                                              allPeaksSorted.indices.cbegin() + endInd);
        
        EICPoint curPoint;
        size_t curScanInd = allPeaksSorted.indices.size(), prvScanInd = allPeaksSorted.indices.size();
        // for frame summing & smoothing
        IMSFrameSummer frameSummer(sumEIMs);
        std::vector<SpectrumRawTypes::Mass> curPointMZs;
        std::vector<SpectrumRawTypes::Mobility> curPointMobs;
        std::vector<SpectrumRawTypes::Intensity> curPointInts;
        // for summing spectra in an IMS frame
        std::map<SpectrumRawTypes::Mass, SpectrumRawTypes::Intensity> curPointMS;
        
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
                if (curPoint.intensity > 0)
                {
                    curPoint.time = curTime;
                    if (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ)
                    {
                        if (anySpecHasMob)
                        {
                            for (const auto &p : curPointMS)
                            {
                                curPoint.mz += (p.first * p.second);
                                if (curPoint.intensityBP == 0.0 || p.second > curPoint.intensityBP)
                                {
                                    curPoint.intensityBP = p.second;
                                    curPoint.mzBP = p.first;
                                }
                            }
                            curPointMS.clear();
                        }
                        
                        curPoint.mz /= curPoint.intensity;
                        eic.addPoint(curPoint, anySpecHasMob, saveEIMs);
                        
                        if (anySpecHasMob)
                        {
                            SpectrumRawTypes::Scan sc = (prvScanInd == allPeaksSorted.indices.size()) ? curScanInd : (prvScanInd + 1);
                            
#ifdef LOG_EIM
                            Rcpp::Rcout << "EIC sum: " << i << "/" << j << "/" << sc << "/" << scanInd
                                        << "/" << frameSummer.size() << "/" <<  specMeta.first.times[sc] << "/" << curTime << std::endl;
#endif
                            
                            while (true)
                            {
                                if (sc != curScanInd)
                                {
                                    // add zero-intensity points for missing scans
                                    frameSummer.addZero();
#ifdef LOG_EIM
                                    Rcpp::Rcout << "added zero:" << i << "/" << j << "/" << sc << "/" << curScanInd << std::endl;
#endif
                                }
                                else
                                {
                                    if (eicMode == EICMode::FULL)
                                        frameSummer.add(std::move(curPointMZs), std::move(curPointMobs),
                                                        std::move(curPointInts));
                                    else
                                        frameSummer.add(std::move(curPointMZs), std::move(curPointInts));
#ifdef LOG_EIM
                                    Rcpp::Rcout << "added points:" << i << "/" << j << "/" << curPointMZs.size() << "/" << curScanInd << std::endl;
#endif
                                }
                                
                                // sumEIMs == 3 --> flank == 1 --> set position 1 back from current
                                // sumEIMs == 5 --> flank == 2 --> set position 2 back from current
                                // etc
                                
                                if (!eic.empty() && sc >= EIMFlank)
                                    updateSummedEIMs(sc - EIMFlank, frameSummer, eic, mobStart, mobEnd);
#ifdef LOG_EIM
                                Rcpp::Rcout << "EIC update: " << specMeta.first.times[sc] << "/" << j << "/" << sc << "/" << curScanInd
                                            << "/" << prvScanInd << "/" << frameSummer.size() << "/" << frameSummer.sizeNoZero() << std::endl;
#endif
                                
                                if (sc == curScanInd)
                                    break;
                                
                                if (sc >= (prvScanInd + sumEIMs))
                                    sc = curScanInd; // jump to current: doesn't make sense to add more zeroes
                                else
                                    ++sc;
                            }
                        }
                    }
                    else
                        eic.addPoint(curPoint.time, curPoint.intensity);
                    if (curPoint.intensity > maxInten)
                        maxInten = curPoint.intensity;
                }
                
                if ((!enoughTimeAboveThr && minEICAdjTime > 0.0) || (!enoughPointsAboveThr && minEICAdjPoints > 0))
                {
                    // if first point, gap in scans, or below intensity threshold, reset
                    if (prvScanInd == allPeaksSorted.indices.size() || prvScanInd != (curScanInd - 1) ||
                        !numberGTE(curPoint.intensity, minEICAdjIntensity))
                    {
                        // Rcpp::Rcout << "Stop! EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
                        startTimeAboveThr = 0.0;
                        adjPointsAboveThr = 0;   
                    }
                    else
                    {
                        if (minEICAdjTime > 0.0)
                        {
                            if (startTimeAboveThr == 0.0)
                                startTimeAboveThr = curTime;
                            else if (numberGTE(curTime - startTimeAboveThr, minEICAdjTime))
                                enoughTimeAboveThr = true;
                        }
                        if (minEICAdjPoints > 0)
                        {
                            ++adjPointsAboveThr;
                            // Rcpp::Rcout << "above: EIC: " << i << "/" << j << "/" << specMeta.first.scans[curScanInd] << "/" << curTime << "/" << curPoint.intensity << "/" << minEICAdjPoints << "/" << adjPointsAboveThr << "\n";
                            if (adjPointsAboveThr >= minEICAdjPoints)
                                enoughPointsAboveThr = true;
                        }
                    }
                }
                
                if (curPoint.intensity > 0)
                    prvScanInd = curScanInd;
                if (ended)
                    break;
                
                curScanInd = scanInd;
                curPoint = EICPoint();
                curPointMZs.clear();
                curPointMobs.clear();
                curPointInts.clear();
            }
            
            const auto time = specMeta.first.times[scanInd];
            if (!timeGTE(time, timeStart) || (timeEnd != 0.0 && !timeLTE(time, timeEnd)))
                continue;
            
            //Rcpp::Rcout << "EIC: " << allPeaksSortedInd << "/" << j << "/" << specMeta.first.times[allPeaksSorted.indices[allPeaksSortedInd]] << "/" << mz << "/" << mzStart << "/" << mzEnd << "\n";
            const auto mz = allPeaksSorted.mzs[allPeaksSortedInd];
            const auto mob = (anySpecHasMob) ? allPeaksSorted.mobilities[allPeaksSortedInd] : 0;
            const auto inten = allPeaksSorted.intensities[allPeaksSortedInd];
            if (anySpecHasMob && (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ))
            {
                // first check if outside smooth window.
                // NOTE: mobExtStart/mobExtEnd == mobStart/mobEnd if no smoothing is performed
                if (mz < mzExtStart || (mzExtEnd != 0.0 && mz > mzExtEnd) ||
                    (mob < mobExtStart || (mobExtEnd != 0.0 && mob > mobExtEnd)))
                    continue;
                
                // add data to frame summer. NOTE: mz and/or at this point may be outside target ranges, but should be
                // within smoothing extension
                curPointMZs.push_back(mz);
                if (eicMode == EICMode::FULL)
                    curPointMobs.push_back(mob);
                curPointInts.push_back(inten);
            }
            
            if (mz < mzStart || (mzEnd != 0.0 && mz > mzEnd) ||
                (anySpecHasMob && mob < mobStart || (mobEnd != 0.0 && mob > mobEnd)))
                continue; // outside target mz or mobility range
            
            if (init)
            {
                curScanInd = scanInd;
                init = false;
            }
            
            if (inten > 0)
            {
                curPoint.intensity += inten;
                if (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ)
                {
                    if (curPoint.mzMin == 0.0 || mz < curPoint.mzMin)
                        curPoint.mzMin = mz;
                    if (mz > curPoint.mzMax)
                        curPoint.mzMax = mz;
                    if (anySpecHasMob)
                        curPointMS[mz] += inten;
                    else
                    {
                        curPoint.mz += mz * inten;
                        if (inten > curPoint.intensityBP)
                        {
                            curPoint.intensityBP = inten;
                            curPoint.mzBP = mz;
                        }
                    }
                    
                    if (anySpecHasMob && eicMode != EICMode::FULL_MZ)
                    {
                        if (curPoint.mobMin == 0.0 || mob < curPoint.mobMin)
                            curPoint.mobMin = mob;
                        if (mob > curPoint.mobMax)
                            curPoint.mobMax = mob;
                    }
                    // Rcpp::Rcout << "EIC: " << curPoint.time << "/" << j << "\t" << std::fixed << mz << "\t" << inten << "\t" << curPoint.mzMin << "/" << curPoint.mzMax << "\t" << mob << std::endl;
                }
            }
        }
        
        if (minEICIntensity != 0.0 && maxInten < minEICIntensity)
            continue;
        if (!enoughTimeAboveThr && minEICAdjTime > 0.0)
            continue;
        if (!enoughPointsAboveThr && minEICAdjPoints > 0)
            continue;
        
        if (eicMode == EICMode::SIMPLE && pad && !eic.empty())
        {
            const auto padded = padEIC(specMeta.first.times, timeStart,
                                       (timeEnd == 0) ? specMeta.first.times.back() : timeEnd, eic.times,
                                       eic.intensities);
            eic.times = std::move(padded.first);
            eic.intensities = std::move(padded.second);
        }
        
        // fill in left-over trailing scans
        // NOTE: at this point, prvScanInd refers to the last filled in scan
        for (SpectrumRawTypes::Scan s=prvScanInd+1-std::min(EIMFlank, static_cast<unsigned>(frameSummer.size()));
             !frameSummer.empty() && s<specMeta.first.times.size(); ++s)
        {
#ifdef LOG_EIM
            Rcpp::Rcout << "EIM purge: " << s << "/" << curScanInd << "/" << specMeta.first.times[s] << "/"
                        << frameSummer.size() << std::endl;
#endif
            if (specMeta.first.times[s] > eic.times.back())
                break; // nothing more to set
            updateSummedEIMs(s, frameSummer, eic, mobStart, mobEnd);
            frameSummer.pop();
        }
        
        allEICs[i] = std::move(eic);
        allEICMaxIntensities[i] = maxInten;
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
    Rcpp::List ret(EICCount);
    for (size_t i=0; i<EICCount; ++i)
    {
        auto &eic = allEICs[i];
        
        if (allEICMaxIntensities[i] < minMaxIntens)
        {
            // just clear the EIC so an empty will be added below. Note that we would clear the EIC at the end anyway.
            eic.clear();
        }
        
        // NOTE: constructing the ouput list directly with create() is much faster than assigning items afterwards.
        // Hence, we do quite some if statements here...
        
        if (eicMode == EICMode::SIMPLE)
        {
            auto mat = Rcpp::NumericMatrix(eic.times.size(), 2);
            for (size_t i=0; i<eic.times.size(); ++i)
            {
                mat(i, 0) = eic.times[i];
                mat(i, 1) = eic.intensities[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity");
            ret[i] = mat;
        }
        else if (eicMode == EICMode::FULL && anySpecHasMob)
        {
            auto mat = Rcpp::NumericMatrix(eic.times.size(), 11);
            for (size_t i=0; i<eic.times.size(); ++i)
            {
                mat(i, 0) = eic.times[i];
                mat(i, 1) = eic.intensities[i];
                mat(i, 2) = eic.intensitiesBP[i];
                mat(i, 3) = eic.mzs[i];
                mat(i, 4) = eic.mzsBP[i];
                mat(i, 5) = eic.mzMins[i];
                mat(i, 6) = eic.mzMaxs[i];
                mat(i, 7) = eic.mobilities[i];
                mat(i, 8) = eic.mobMins[i];
                mat(i, 9) = eic.mobMaxs[i];
                mat(i, 10) = eic.mobilitiesBP[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity", "intensityBP", "mz", "mzBP",
                           "mzmin", "mzmax", "mobility", "mobmin", "mobmax",
                           "mobilityBP");
            ret[i] = mat;
        }
        else if (eicMode == EICMode::FULL || eicMode == EICMode::FULL_MZ)
        {
            auto mat = Rcpp::NumericMatrix(eic.times.size(), 7);
            for (size_t i=0; i<eic.times.size(); ++i)
            {
                mat(i, 0) = eic.times[i];
                mat(i, 1) = eic.intensities[i];
                mat(i, 2) = eic.intensitiesBP[i];
                mat(i, 3) = eic.mzs[i];
                mat(i, 4) = eic.mzsBP[i];
                mat(i, 5) = eic.mzMins[i];
                mat(i, 6) = eic.mzMaxs[i];
            }
            Rcpp::colnames(mat) = Rcpp::CharacterVector::create("time", "intensity", "intensityBP", "mz", "mzBP",
                           "mzmin", "mzmax");
            ret[i] = mat;
        }
        else // if (eicMode == EICMode::TEST)
        {
            ret[i] = !eic.empty();
        }
        
        eic.clear(); // free memory as EICs may consume a lot
    }
    
    if (eicMode != EICMode::TEST)
    {
        if (specMeta.first.times.size() > 1 && gapFactor > 0.0)
        {
            // NOTE: Bruker TIMS data (and maybe others?) seem to omit zero intensity scans, leading to time gaps. Since
            // this leads to incorrect EICs, we pad here. To detect gaps, we take the median difference between scans
            // and assume a gap is at least X higher than that. For padding we just have to add additional time points,
            // EIC filling and padding functions will assume these are zero intensity points.
            std::vector<SpectrumRawTypes::Time> diffs(specMeta.first.times.size() - 1, 0);
            for (size_t i=1; i<specMeta.first.times.size(); ++i)
                diffs[i - 1] = specMeta.first.times[i] - specMeta.first.times[i - 1];
            const double medianDiff = median(diffs);
            std::vector<SpectrumRawTypes::Time> allXValues;
            for (size_t i=0; i<specMeta.first.times.size(); ++i)
            {
                allXValues.push_back(specMeta.first.times[i]);
                if (i == (specMeta.first.times.size() - 1))
                    break;
                
                const auto diff = specMeta.first.times[i + 1] - specMeta.first.times[i];
                if (diff > medianDiff * gapFactor) // add dummy point after current
                {
                    allXValues.push_back(specMeta.first.times[i] + medianDiff);
                    if (diff > medianDiff * 2.0) // add dummy point before next
                        allXValues.push_back(specMeta.first.times[i + 1] - medianDiff);
                }
            }
            ret.attr("allXValues") = allXValues;
        }
        else
            ret.attr("allXValues") = specMeta.first.times;
    }
    
    if (eicMode == EICMode::FULL && anySpecHasMob && saveEIMs)
    {
        // add EIMs as attribute
        Rcpp::List allEIMs(EICCount);
        for (size_t i=0; i<EICCount; ++i)
        {
            const auto &eic = allEICs[i];
            
            if (eic.EIMs.empty())
                continue;
            
            Rcpp::List eimList(eic.EIMs.size());
            for (size_t j=0; j<eic.EIMs.size(); ++j)
            {
                Rcpp::NumericMatrix emat(eic.EIMs[j].first.size(), 2);
                for (size_t k=0; k<eic.EIMs[j].first.size(); ++k)
                {
                    emat(k, 0) = eic.EIMs[j].first[k];
                    emat(k, 1) = eic.EIMs[j].second[k];
                }
                Rcpp::colnames(emat) = Rcpp::CharacterVector::create("mobility", "intensity");
                eimList[j] = emat;
            }
            allEIMs[i] = eimList;
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
