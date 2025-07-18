#include <Rcpp.h>

#include <cmath>
#include <set>

#include "utils.h"

namespace{

struct FeatureDim
{
    double ret, mz, mob;
    FeatureDim(void) = default;
    FeatureDim(double r, double m, double mb) : ret(r), mz(m), mob(mb) { }
};

struct Feature
{
    FeatureDim dims;
    int anaID;
    size_t featID;
    Feature(void) = default;
    Feature(double r, double m, double mb, int ana, size_t id) : dims(r, m, mb), anaID(ana), featID(id) { }
};

using ScoreWeights = FeatureDim;

FeatureDim getFeatureDim(const std::vector<Feature> &group)
{
    double minRT = 0.0, maxRT = 0.0, minMZ = 0.0, maxMZ = 0.0, minMob = 0.0, maxMob = 0.0;
    bool init = true;
    for (const auto &feat : group)
    {
        if (init)
        {
            minRT = maxRT = feat.dims.ret;
            minMZ = maxMZ = feat.dims.mz;
            minMob = maxMob = feat.dims.mob;
            init = false;
        }
        else
        {
            minRT = std::min(minRT, feat.dims.ret); maxRT = std::max(maxRT, feat.dims.ret);
            minMZ = std::min(minMZ, feat.dims.mz); maxMZ = std::max(maxMZ, feat.dims.mz);
            minMob = std::min(minMob, feat.dims.mob); maxMob = std::max(maxMob, feat.dims.mob);
        }
    }
    return FeatureDim{(maxRT - minRT), (maxMZ - minMZ), (maxMob - minMob)};
}

double calcFeatureDist(const Feature &feat1, const Feature &feat2, double rtWindow, double mzWindow, double mobWindow)
{
    // calculate distance between two features
    const double retDiff = std::abs(feat1.dims.ret - feat2.dims.ret);
    const double mzDiff = std::abs(feat1.dims.mz - feat2.dims.mz);
    const double mobDiff = std::abs(feat1.dims.mob - feat2.dims.mob);
    
    return std::sqrt(std::pow(retDiff / rtWindow, 2.0) +
                     std::pow(mzDiff / mzWindow, 2.0) +
                     std::pow(mobDiff / mobWindow, 2.0));
}

double calcGroupScore(const std::vector<Feature> &tentativeGroup, const Feature &anaFeat, double rtWindow,
                      double mzWindow, double mobWindow, const ScoreWeights &weights)
{
    std::vector<Feature> otherFeatures; // all features from other analyses, and only closest if >1 in analysis
    otherFeatures.push_back(anaFeat); // add the feature itself to the group
    
    std::set<int> visitedAnas;
    visitedAnas.insert(anaFeat.anaID);
    
    double totRTDev = 0.0, totMZDev = 0.0, totMobDev = 0.0;
    
    for (const auto &feat : tentativeGroup)
    {
        if (visitedAnas.find(feat.anaID) != visitedAnas.end())
            continue;

        visitedAnas.insert(feat.anaID);
        
        // take the best if there are more features in this analysis
        double closestDist = -1.0;
        Feature closestFeat = feat;
        for (const auto &sameAnaFeat : tentativeGroup)
        {
            if (sameAnaFeat.anaID != feat.anaID || sameAnaFeat.featID == anaFeat.featID)
                continue;
            
            const double dist = calcFeatureDist(anaFeat, sameAnaFeat, rtWindow, mzWindow, mobWindow);
            if (closestDist == -1.0 || dist < closestDist)
            {
                closestDist = dist;
                closestFeat = sameAnaFeat;
            }
        }
        
        if (closestDist != -1.0)
        {
            otherFeatures.push_back(closestFeat);
            totRTDev += std::abs(anaFeat.dims.ret - closestFeat.dims.ret) / rtWindow;
            totMZDev += std::abs(anaFeat.dims.mz - closestFeat.dims.mz) / mzWindow;
            totMobDev += std::abs(anaFeat.dims.mob - closestFeat.dims.mob) / mobWindow;
        }
    }
    
    // don't consider size==1 groups: these are undesired and if the feature stays unassigned, it will be put in a
    // separate group afterwards in groupFeaturesGreedy()
    if (otherFeatures.size() <= 1)
        return -1.0;
    
    const double dblSize = static_cast<double>(otherFeatures.size());
    const double retScore = 1.0 - (totRTDev / dblSize);
    const double mzScore = 1.0 - (totMZDev / dblSize);
    const double mobScore = 1.0 - (totMobDev / dblSize);
    
#if 0
    if (false)
    {
        Rcpp::Rcout << "anaFeat:" << anaFeat.anaID << " scores: " << "ret: " << retScore << " mz: " << mzScore
                    << "mob: " << mobScore << ", " << "size: " << sizeScore << " total:"
                    << retScore * weights.ret + mzScore * weights.mz + mobScore * weights.mob << "\n";
        Rcpp::Rcout << "tentative group size: " << tentativeGroup.size() << ", other features size: " << otherFeatures.size() << "\n";
        Rcpp::Rcout << "group:\n";
        for (const auto &feat : otherFeatures)
        {
            Rcpp::Rcout << "  featID: " << feat.featID << ", anaID: " << feat.anaID
                        << ", RT: " << feat.dims.ret << ", MZ: " << feat.dims.mz
                        << ", MOB: " << feat.dims.mob << "\n";
        }
    }
#endif
    
    return retScore * weights.ret + mzScore * weights.mz + mobScore * weights.mob;
}

std::vector<int> checkTentativeGroup(const std::vector<Feature> &tentativeGroup, double rtWindow, double mzWindow,
                                     double mobWindow, const ScoreWeights &weights)
{
    std::vector<int> ret;
    if (tentativeGroup.empty())
        return ret;
    
    // NOTE: use a vector to keep the order the same
    std::vector<int> allAnaIDs;
    for (const auto &feat : tentativeGroup)
    {
        if (std::find(allAnaIDs.begin(), allAnaIDs.end(), feat.anaID) == allAnaIDs.end())
            allAnaIDs.push_back(feat.anaID);
    }
    
    for (int anaID : allAnaIDs)
    {
        double highestScore = -1.0;
        size_t bestFeatID;
        for (const auto &anaFeat : tentativeGroup)
        {
            if (anaFeat.anaID != anaID)
                continue;
            
            const auto score = calcGroupScore(tentativeGroup, anaFeat, rtWindow, mzWindow, mobWindow, weights);
            if (score > highestScore)
            {
                highestScore = score;
                bestFeatID = anaFeat.featID; // store the feature ID of the best feature
            }
        }
        if (highestScore > 0.0)
            ret.push_back(bestFeatID); // add this feature to the group
    }
    
    return ret;
}

}

// [[Rcpp::export]]
Rcpp::IntegerVector getGroupIDs(const Rcpp::NumericVector &featRTs, const Rcpp::NumericVector  &featMZs,
                                const Rcpp::NumericVector &featMobs, const Rcpp::NumericVector &ints,
                                const Rcpp::IntegerVector &anaIDs, double rtWindow, double mzWindow,
                                double mobWindow, const Rcpp::List &weightsList)
{
    const ScoreWeights weights(weightsList["retention"], weightsList["mz"], weightsList["mobility"]);
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
            
            tentativeGroup.emplace_back(featRTs[groupInd], featMZs[groupInd], featMobs[groupInd], anaIDs[groupInd],
                                        groupInd);
        }
        
        std::sort(tentativeGroup.begin(), tentativeGroup.end(),
            [&ints](const Feature &a, const Feature &b) { return ints[a.featID] > ints[b.featID]; });
        
        const auto groupIDs = checkTentativeGroup(tentativeGroup, rtWindow, mzWindow, mobWindow, weights);
        if (groupIDs.empty())
            continue;
        
        for (const auto &gID : groupIDs)
            ret[gID] = curGroup;
        
        ++curGroup;
    }
    
    // Fix unassigned groups
    for (size_t i=0; i<ret.size(); ++i)
    {
        if (ret[i] == -1)
        {
            ret[i] = curGroup;
            ++curGroup;
        }
    }
    
    return ret;
}
