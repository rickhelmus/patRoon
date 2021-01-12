#include <iomanip>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Rcpp.h>

#include "parse-xml.h"

namespace {

struct adductCons
{
    numType ret, mz;
    std::vector<int> charges;
    std::vector<std::string> featIDs, adducts;
};
    
adductCons parseAddConsXMLBlock(pugi::xml_document &doc)
{
    pugi::xml_node consNode = doc.child("consensusElement");
    
    adductCons ret;
    
    auto centr = consNode.child("centroid");
    ret.ret = getNumericFromXML(centr.attribute("rt"));
    ret.mz = getNumericFromXML(centr.attribute("mz"));
    
    for (auto el : consNode.child("groupedElementList").children("element"))
    {
        ret.featIDs.push_back(el.attribute("id").value());
        ret.charges.push_back(el.attribute("charge").as_int());
    }
    
    ret.adducts.resize(ret.featIDs.size());
    for (auto up : consNode.children("UserParam"))
    {
        const std::string type = up.attribute("type").value();
        if (type != "string")
            continue;
        const std::string id = up.attribute("name").value();
        for (std::size_t i=0; i<ret.featIDs.size(); ++i)
        {
            if (ret.featIDs[i] == id)
            {
                ret.adducts[i] = up.attribute("value").value();
                break;
            }
        }
    }
 
    return ret;
}
    
}


// [[Rcpp::export]]
Rcpp::List parseAdductConsXMLFile(Rcpp::CharacterVector file)
{
    std::vector<adductCons> consElements;
    
    parseXMLFile(Rcpp::as<const char *>(file), "<consensusElement id=", "</consensusElement>",
                 [&](pugi::xml_document &doc)
                 {
                     adductCons fc(parseAddConsXMLBlock(doc));
                     consElements.push_back(fc);
                 });
    
    Rcpp::List ret;
    for (auto ce : consElements)
    {
        ret.push_back(Rcpp::DataFrame::create(Rcpp::Named("avgRet") = ce.ret,
                                              Rcpp::Named("neutralMass") = ce.mz,
                                              Rcpp::Named("ID") = ce.featIDs,
                                              Rcpp::Named("charge") = ce.charges,
                                              Rcpp::Named("adduct") = ce.adducts));
    }
    
    return ret;    
}
