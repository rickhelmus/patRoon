/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <string>
#include <vector>

#include <Rcpp.h>

#include "utils-xml.h"

namespace {

struct featureMRM
{
    std::string ID, chromID;
    numType ret, area, intensity;
    numType retMin, retMax;
    numType tailingFactor, asymmetryFactor, slopeOfBaseline;
    int pointsAcrossBaseline, pointsAcrossHalfHeight;
};

featureMRM parseFeatureXMLBlock(pugi::xml_document &doc)
{
    pugi::xml_node featNode = doc.child("feature");
    
    featureMRM ret;
    ret.ID = featNode.attribute("id").value();
    ret.area = getNumericFromXML(featNode.child("intensity").text()); // intensity reported by OpenMS is actually the area
    ret.retMin = getNumericFromXML(featNode.find_child_by_attribute("UserParam", "name", "leftWidth").attribute("value"));
    ret.retMax = getNumericFromXML(featNode.find_child_by_attribute("UserParam", "name", "rightWidth").attribute("value"));
    
    for (auto sub : featNode.child("subordinate").children("feature"))
    {
        ret.ret = getNumericFromXML(sub.find_child_by_attribute("UserParam", "name", "peak_apex_position").attribute("value"));
        ret.intensity = getNumericFromXML(sub.find_child_by_attribute("UserParam", "name", "peak_apex_int").attribute("value"));
        ret.chromID = sub.find_child_by_attribute("UserParam", "name", "native_id").attribute("value").as_string();
        ret.tailingFactor = getNumericFromXML(sub.find_child_by_attribute("UserParam", "name", "tailing_factor").attribute("value"));
        ret.asymmetryFactor = getNumericFromXML(sub.find_child_by_attribute("UserParam", "name", "asymmetry_factor").attribute("value"));
        ret.slopeOfBaseline = getNumericFromXML(sub.find_child_by_attribute("UserParam", "name", "slope_of_baseline").attribute("value"));
        ret.pointsAcrossBaseline = sub.find_child_by_attribute("UserParam", "name", "points_across_baseline").attribute("value").as_int();
        ret.pointsAcrossHalfHeight = sub.find_child_by_attribute("UserParam", "name", "points_across_half_height").attribute("value").as_int();
    }
    
    return ret;
}

}

// [[Rcpp::export]]
Rcpp::DataFrame parseFeatureMRMXMLFile(Rcpp::CharacterVector file)
{
    std::vector<std::string> IDs, chromIDs;
    std::vector<numType> areas, ints, rets, retMins, retMaxs, tFactors, aFactors, slopesBase;
    std::vector<int> pointsBase, pointsHH;
    
    parseXMLFile(Rcpp::as<const char *>(file), "<feature id=", "</feature>",
                 [&](pugi::xml_document &doc)
                 {
                     featureMRM f(parseFeatureXMLBlock(doc));
                     IDs.push_back(f.ID);
                     chromIDs.push_back(f.chromID);
                     areas.push_back(f.area);
                     ints.push_back(f.intensity);
                     rets.push_back(f.ret);
                     retMins.push_back(f.retMin);
                     retMaxs.push_back(f.retMax);
                     tFactors.push_back(f.tailingFactor);
                     aFactors.push_back(f.asymmetryFactor);
                     slopesBase.push_back(f.slopeOfBaseline);
                     pointsBase.push_back(f.pointsAcrossBaseline);
                     pointsHH.push_back(f.pointsAcrossHalfHeight);
                 });
    
    return Rcpp::DataFrame::create(Rcpp::Named("ID") = IDs,
                                   Rcpp::Named("chromID") = chromIDs,
                                   Rcpp::Named("ret") = rets,
                                   Rcpp::Named("area") = areas,
                                   Rcpp::Named("intensity") = ints,
                                   Rcpp::Named("retmin") = retMins,
                                   Rcpp::Named("retmax") = retMaxs,
                                   Rcpp::Named("tailingFactor") = tFactors,
                                   Rcpp::Named("asymmetryFactor") = aFactors,
                                   Rcpp::Named("slopeOfBaseline") = slopesBase,
                                   Rcpp::Named("pointsAcrossBaseline") = pointsBase,
                                   Rcpp::Named("pointsAcrossHalfHeight") = pointsHH,
                                   Rcpp::Named("stringsAsFactors") = false);
}
