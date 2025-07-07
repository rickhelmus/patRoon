#include <Rcpp.h>

#include <cmath>

#include "utils.h"

namespace{

struct Feature
{
    double ret, mz, mob;
    int anaID;
    size_t featID;
    Feature(void) = default;
    Feature(double rt, double mzVal, double mobVal, int ana, size_t id) : ret(rt), mz(mzVal), mob(mobVal), anaID(ana), featID(id) { }
};

struct FeatureDim
{
    double ret, mz, mob;
    FeatureDim(void) = default;
    FeatureDim(double r, double m, double mb) : ret(r), mz(m), mob(mb) { }
};

FeatureDim getFeatureDim(const std::vector<Feature> &group)
{
    double minRT, maxRT, minMZ, maxMZ, minMob, maxMob;
    bool init = true;
    for (const auto &feat : group)
    {
        if (init)
        {
            minRT = maxRT = feat.ret;
            minMZ = maxMZ = feat.mz;
            minMob = maxMob = feat.mob;
            init = false;
        }
        else
        {
            minRT = std::min(minRT, feat.ret); maxRT = std::max(maxRT, feat.ret);
            minMZ = std::min(minMZ, feat.mz); maxMZ = std::max(maxMZ, feat.mz);
            minMob = std::min(minMob, feat.mob); maxMob = std::max(maxMob, feat.mob);
        }
    }
    return FeatureDim{(maxRT - minRT), (maxMZ - minMZ), (maxMob - minMob)};
}

double calcFeatureDist(const Feature &feat1, const Feature &feat2, double rtWindow, double mzWindow, double mobWindow)
{
    // calculate distance between two features
    const double retDiff = std::abs(feat1.ret - feat2.ret);
    const double mzDiff = std::abs(feat1.mz - feat2.mz);
    const double mobDiff = std::abs(feat1.mob - feat2.mob);
    
    return std::sqrt(std::pow(retDiff / rtWindow, 2.0) +
                     std::pow(mzDiff / mzWindow, 2.0) +
                     std::pow(mobDiff / mobWindow, 2.0));
}

double getFeatGroupVolume(const std::vector<Feature> &group, double rtWindow, double mzWindow, double mobWindow)
{
    if (group.empty())
        return 0.0;

    FeatureDim dim = getFeatureDim(group);
    
    const double eps = 1E-8; // small epsilon to avoid division by zero
    const double dimRT = (dim.ret / rtWindow) + eps;
    const double dimMZ = (dim.mz / mzWindow) + eps;
    const double dimMob = (dim.mob / mobWindow) + eps;
    return dimRT * dimMZ * dimMob;
}
    
double calcGroupScore(const std::vector<Feature> &tentativeGroup, const Feature &anaFeat, double rtWindow,
                      double mzWindow, double mobWindow)
{
    std::vector<Feature> otherFeatures; // all features from other analyses, and only closest if >1 in analysis
    for (const auto &feat : tentativeGroup)
    {
        if (feat.anaID == anaFeat.anaID)
            continue;

        // take the best if there are more features in this analysis
        double closestDist = std::numeric_limits<double>::max();
        Feature closestFeat = feat;
        for (const auto &sameAnaFeat : tentativeGroup)
        {
            if (sameAnaFeat.anaID != feat.anaID || sameAnaFeat.featID == anaFeat.featID)
                continue;
            
            const double dist = calcFeatureDist(anaFeat, sameAnaFeat, rtWindow, mzWindow, mobWindow);
            if (dist < closestDist)
            {
                closestDist = dist;
                closestFeat = sameAnaFeat;
            }
        }
        
        otherFeatures.push_back(closestFeat);
    }
    otherFeatures.push_back(anaFeat); // add the feature itself to the group
    
    // check if this feature would lead to an invalid group
    FeatureDim dim = getFeatureDim(otherFeatures);
    if (dim.ret > (rtWindow * 2.0) || dim.mz > (mzWindow * 2.0) || dim.mob > (mobWindow * 2.0))
        return -1.0;
    
    const double groupVolume = getFeatGroupVolume(otherFeatures, rtWindow, mzWindow, mobWindow);
    return 1.0 / groupVolume + static_cast<double>(otherFeatures.size()) * 0.1; // UNDONE: normalize group size?
}

std::vector<int> checkTentativeGroup(const std::vector<Feature> &tentativeGroup, double rtWindow, double mzWindow,
                                     double mobWindow)
{
    std::vector<int> ret;
    if (tentativeGroup.empty())
        return ret;
    
    std::set<int> allAnaIDs;
    for (const auto &feat : tentativeGroup)
        allAnaIDs.insert(feat.anaID);

    for (int anaID : allAnaIDs)
    {
        double highestScore = -1.0;
        size_t bestFeatID;
        for (const auto &anaFeat : tentativeGroup)
        {
            if (anaFeat.anaID != anaID)
                continue;
            
            const auto score = calcGroupScore(tentativeGroup, anaFeat, rtWindow, mzWindow, mobWindow);
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
                                double mobWindow, bool verbose)
{
    Rcpp::IntegerVector ret(featRTs.size(), -1);
    int curGroup = 0;
    const auto intSortedInds = getSortedInds(ints, true);
    const auto mzSortedInds = getSortedInds(featMZs);
    
    if (verbose)
        Rcpp::Rcout << "Grouping " << featRTs.size() << " features...\n";
    
    for (size_t i=0; i<intSortedInds.size(); ++i)
    {
        if (verbose && i > 0 && (i % 5000) == 0)
            Rcpp::Rcout << "Processed " << i << " features...\n";
        
        const auto refInd = intSortedInds[i];
        
        if (ret[refInd] != -1)
            continue; // already assigned to a group
        
        const auto refRT = featRTs[refInd], refMZ = featMZs[refInd], refMob = featMobs[refInd];
        
        auto groupStartIndIt = std::lower_bound(mzSortedInds.begin(), mzSortedInds.end(), refInd,
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
        
        const auto groupIDs = checkTentativeGroup(tentativeGroup, rtWindow, mzWindow, mobWindow);
        
        for (const auto &gID : groupIDs)
            ret[gID] = curGroup;
        
        ++curGroup;
    }
    
    if (verbose)
        Rcpp::Rcout << "Done!\n";
    
    return ret;
}
