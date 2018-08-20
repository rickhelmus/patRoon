#include <string>
#include <vector>

#include <Rcpp.h>

#include "parse-xml.h"

namespace {

struct featCons
{
    struct item
    {
        int fileID, featID;
    };
    
    std::vector<item> items;
    numType ret, mz;
};

featCons parseFeatConsXMLBlock(pugi::xml_document &doc)
{
    pugi::xml_node consNode = doc.child("consensusElement");
    
    featCons ret;
    auto centr = consNode.child("centroid");
    ret.ret = getNumericFromXML(centr.attribute("rt"));
    ret.mz = getNumericFromXML(centr.attribute("mz"));
    
    for (auto el : consNode.child("groupedElementList").children("element"))
    {
        const int map = el.attribute("map").as_int();
        const int id = el.attribute("id").as_int();
        ret.items.push_back(featCons::item{map, id});
    }
    
    return ret;
}

}


// [[Rcpp::export]]
Rcpp::List parseFeatConsXMLFile(Rcpp::CharacterVector file, Rcpp::IntegerVector anaCount)
{
    std::vector<featCons> consElements;
    
    parseXMLFile(Rcpp::as<const char *>(file), "<consensusElement id=", "</consensusElement>",
                 [&](pugi::xml_document &doc)
                 {
                     featCons fc(parseFeatConsXMLBlock(doc));
                     consElements.push_back(fc);
                 });
    
    const auto gCount = consElements.size();
    const int aCount = Rcpp::as<int>(anaCount);
    std::vector<numType> rets(gCount), mzs(gCount);
    Rcpp::IntegerMatrix ftindex(aCount, gCount);
    
    for (std::size_t gi=0; gi<gCount; ++gi)
    {
        rets[gi] = consElements[gi].ret;
        mzs[gi] = consElements[gi].mz;
        for (const auto &el : consElements[gi].items)
            ftindex(el.fileID, gi) = el.featID;
    }
    
    Rcpp::DataFrame gInfo = Rcpp::DataFrame::create(Rcpp::Named("rts") = rets,
                                                    Rcpp::Named("mzs") = mzs);
    
    return Rcpp::List::create(Rcpp::Named("gInfo") = gInfo,
                              Rcpp::Named("ftindex") = ftindex);
}
