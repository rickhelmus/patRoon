#include <Rcpp.h>

#include <map>
#include <queue>

#include "utils.h"

namespace{

struct FeatureDim
{
    double ret, mz, mob, inten;
    FeatureDim(void) = default;
    FeatureDim(double r, double m, double mb, double i) : ret(r), mz(m), mob(mb), inten(i) { }
};

struct Feature
{
    FeatureDim dims;
    int anaID, repID;
    size_t featID;
    Feature(void) = default;
    Feature(double r, double m, double mb, double i, int ana, int rep,
            size_t id) : dims(r, m, mb, i), anaID(ana), repID(rep), featID(id) { }
};

using ScoreWeights = FeatureDim;

double calcGroupScore(const std::vector<Feature> &tentativeGroup, double rtWindow, double mzWindow,
                      double mobWindow, const ScoreWeights &weights)
{
    if (tentativeGroup.size() < 2)
        return -1.0;
    
    std::map<int, double> maxRepIntensities;
    for (const auto &feat : tentativeGroup)
    {
        if (maxRepIntensities.find(feat.repID) == maxRepIntensities.end() ||
            feat.dims.inten > maxRepIntensities[feat.repID])
        {
            maxRepIntensities[feat.repID] = feat.dims.inten;
        }
    }
    
    double totRTDev = 0.0, totMZDev = 0.0, totMobDev = 0.0, totIntenScore = 0.0;
    const Feature &refFeat = tentativeGroup[0];
    for (const auto &feat : tentativeGroup)
    {
        totRTDev += std::abs(refFeat.dims.ret - feat.dims.ret) / rtWindow;
        totMZDev += std::abs(refFeat.dims.mz - feat.dims.mz) / mzWindow;
        totMobDev += std::abs(refFeat.dims.mob - feat.dims.mob) / mobWindow;
        totIntenScore += (feat.dims.inten / maxRepIntensities[feat.repID]);
    }
    
    const double dblSize = static_cast<double>(tentativeGroup.size());
    const double retScore = 1.0 - (totRTDev / dblSize);
    const double mzScore = 1.0 - (totMZDev / dblSize);
    const double mobScore = 1.0 - (totMobDev / dblSize);
    const double intenScore = (totIntenScore / dblSize);
    
#if 0
    if (true)
    {
        Rcpp::Rcout << "refFeat:" << refFeat.featID << "/" << refFeat.anaID << " scores: "
                    << "ret: " << retScore << " mz: " << mzScore << "mob: " << mobScore << "inten: " << intenScore << " total:"
                    << retScore * weights.ret + mzScore * weights.mz + mobScore * weights.mob + intenScore * weights.inten << "\n";
        Rcpp::Rcout << "group:\n";
        for (const auto &feat : tentativeGroup)
        {
            Rcpp::Rcout << "  featID: " << feat.featID << ", anaID: " << feat.anaID
                        << ", RT: " << feat.dims.ret << ", MZ: " << feat.dims.mz
                        << ", MOB: " << feat.dims.mob << "\n";
        }
    }
#endif
    
    return retScore * weights.ret + mzScore * weights.mz + mobScore * weights.mob + intenScore * weights.inten;
}


std::vector<Feature> getBestGroup(const std::vector<Feature> &allFeatures, double rtWindow,
                                  double mzWindow, double mobWindow, const ScoreWeights &weights)
{
    if (allFeatures.size() <= 1)
        return allFeatures;

    struct WorkItem
    {
        std::vector<Feature> tentativeGroup;
        size_t nextAnaIndex;
        WorkItem(const std::vector<Feature> &tg, size_t ni) : tentativeGroup(tg), nextAnaIndex(ni) { }
    };
    
    // NOTE: use a vector to keep the order the same
    std::vector<int> allAnaIDs;
    for (const auto &feat : allFeatures)
    {
        if (std::find(allAnaIDs.begin(), allAnaIDs.end(), feat.anaID) == allAnaIDs.end())
            allAnaIDs.push_back(feat.anaID);
    }
    
    std::queue<WorkItem> workQueue;
    std::vector<Feature> bestGroup;
    double bestScore = -1.0;
    
    // initialize queue
    for (const auto &feat : allFeatures)
    {
        if (feat.anaID == allAnaIDs[0])
            workQueue.emplace(std::vector<Feature>{feat}, 1);
    }
    
    while (!workQueue.empty())
    {
        const auto wi = workQueue.front();
        workQueue.pop();

        if (wi.nextAnaIndex >= allAnaIDs.size())
        {
            // reached the end of the analysis IDs, check the group
            const auto groupScore = calcGroupScore(wi.tentativeGroup, rtWindow, mzWindow, mobWindow, weights);
            if (groupScore > bestScore)
            {
                bestScore = groupScore;
                bestGroup = wi.tentativeGroup;
            }
            continue;
        }
        
        for (const auto &feat : allFeatures)
        {
            if (feat.anaID != allAnaIDs[wi.nextAnaIndex])
                continue;
            std::vector<Feature> tg = wi.tentativeGroup;
            tg.push_back(feat);
            workQueue.emplace(tg, wi.nextAnaIndex + 1);
        }
    }
    
    return bestGroup;
}

}

// [[Rcpp::export]]
Rcpp::IntegerVector getGroupIDs(const Rcpp::NumericVector &featRTs, const Rcpp::NumericVector  &featMZs,
                                const Rcpp::NumericVector &featMobs, const Rcpp::NumericVector &ints,
                                const Rcpp::IntegerVector &anaIDs, const Rcpp::IntegerVector &repIDs, double rtWindow,
                                double mzWindow, double mobWindow, const Rcpp::List &weightsList)
{
    const ScoreWeights weights(weightsList["retention"], weightsList["mz"], weightsList["mobility"],
                               weightsList["intensity"]);
    Rcpp::IntegerVector ret(featRTs.size(), -1);
    int curGroup = 0;
    const auto intSortedInds = getSortedInds(ints, true);
    const auto mzSortedInds = getSortedInds(featMZs);
    
    for (size_t i=0; i<intSortedInds.size(); ++i)
    {
        const auto refInd = intSortedInds[i];
        
        if (ret[refInd] != -1)
            continue; // already assigned to a group
        
        const auto refRT = featRTs[refInd], refMZ = featMZs[refInd], refMob = featMobs[refInd];
        
        const auto groupStartIndIt = std::lower_bound(mzSortedInds.begin(), mzSortedInds.end(), refInd,
            [&](size_t ind, size_t val) { return (featMZs[ind] < (refMZ - mzWindow)); });
        
        std::vector<Feature> tentativeGroup;
        for (auto it=groupStartIndIt; it!=mzSortedInds.end(); ++it)
        {
            const auto groupInd = *it;
            
            if (ret[groupInd] != -1)
                continue; // already assigned to a group
            
            if (featMZs[groupInd] > (refMZ + mzWindow))
                break; // no more features in this group
            
            if (std::abs(featRTs[groupInd] - refRT) > rtWindow)
                continue;
            if (std::abs(featMobs[groupInd] - refMob) > mobWindow)
                continue;
            
            tentativeGroup.emplace_back(featRTs[groupInd], featMZs[groupInd], featMobs[groupInd], ints[groupInd],
                                        anaIDs[groupInd], repIDs[groupInd], groupInd);
        }
        
        std::sort(tentativeGroup.begin(), tentativeGroup.end(),
            [&ints](const Feature &a, const Feature &b) { return ints[a.featID] > ints[b.featID]; });
        
        const auto bestGroup = getBestGroup(tentativeGroup, rtWindow, mzWindow, mobWindow, weights);
        if (bestGroup.empty())
            continue;
        
        for (const auto &feat : bestGroup)
            ret[feat.featID] = curGroup;
        
        ++curGroup;
    }
    
    // Fix unassigned groups
    for (int i=0; i<ret.size(); ++i)
    {
        if (ret[i] == -1)
        {
            ret[i] = curGroup;
            ++curGroup;
        }
    }
    
    return ret;
}
