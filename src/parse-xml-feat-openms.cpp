#include <string>
#include <vector>

#include <Rcpp.h>

#include "parse-xml.h"

namespace {

struct feature
{
    std::string id;
    numType ret, mz, intensity;
    numType retMin, retMax;
    numType mzMin, mzMax;
    int isoCount;
};

feature parseFeatureXMLBlock(pugi::xml_document &doc)
{
    pugi::xml_node featNode = doc.child("feature");
    
    feature ret;
    ret.id = featNode.attribute("id").value();
    ret.intensity = getNumericFromXML(featNode.child("intensity").text());
    
    for (auto pos: featNode.children("position"))
    {
        if (pos.attribute("dim").as_int() == 0)
            ret.ret = getNumericFromXML(pos.text());
        else
            ret.mz = getNumericFromXML(pos.text());
    }
    
    ret.isoCount = 0;
    for (auto hull: featNode.children("convexhull"))
    {
        ++ret.isoCount;
        
        // only take first for now
        if (hull.attribute("nr").as_int() == 0)
        {
            bool init = true;
            for (auto pt: hull.children("pt"))
            {
                const numType rt = getNumericFromXML(pt.attribute("x"));
                const numType mz = getNumericFromXML(pt.attribute("y"));
                if (init)
                {
                    init = false;
                    ret.retMin = ret.retMax = rt;
                    ret.mzMin = ret.mzMax = mz;
                }
                else
                {
                    ret.retMin = std::min(ret.retMin, rt);
                    ret.retMax = std::max(ret.retMax, rt);
                    ret.mzMin = std::min(ret.mzMin, mz);
                    ret.mzMax = std::max(ret.mzMax, mz);
                }
            }
        }
    }
    
    return ret;
}

}

// [[Rcpp::export]]
Rcpp::DataFrame parseFeatureXMLFile(Rcpp::CharacterVector file)
{
    std::vector<std::string> ids;
    std::vector<numType> ints, rets, mzs, retMins, retMaxs, mzMins, mzMaxs;
    std::vector<int> isoCounts;
    
    parseXMLFile(Rcpp::as<const char *>(file), "<feature id=", "</feature>",
                 [&](pugi::xml_document &doc)
                 {
                     feature f(parseFeatureXMLBlock(doc));
                     ids.push_back(f.id);
                     ints.push_back(f.intensity);
                     rets.push_back(f.ret);
                     mzs.push_back(f.mz);
                     retMins.push_back(f.retMin);
                     retMaxs.push_back(f.retMax);
                     mzMins.push_back(f.mzMin);
                     mzMaxs.push_back(f.mzMax);
                     isoCounts.push_back(f.isoCount);
                 });
    
    return Rcpp::DataFrame::create(Rcpp::Named("ID") = ids,
                                   Rcpp::Named("ret") = rets,
                                   Rcpp::Named("mz") = mzs,
                                   Rcpp::Named("area") = ints, // intensity reported by OpenMS is actually the area
                                   Rcpp::Named("retmin") = retMins,
                                   Rcpp::Named("retmax") = retMaxs,
                                   Rcpp::Named("mzmin") = mzMins,
                                   Rcpp::Named("mzmax") = mzMaxs,
                                   Rcpp::Named("isocount") = isoCounts,
                                   Rcpp::Named("stringsAsFactors") = false);
}
