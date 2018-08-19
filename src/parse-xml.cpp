#include <fstream>
#include <string>
#include <vector>

#include <Rcpp.h>

#define PUGIXML_HEADER_ONLY
#define PUGIXML_NO_XPATH
#include "pugixml/pugixml.hpp"

#include "utils.h"

typedef double numType;
template<typename T> numType getNumericFromXML(const T &x) { return x.as_double(); }

// startTag: uses partial match (prefix) as tag may contain attributes etc
void parseXMLFile(const char *file, const std::string &startTag,
                  const std::string &endTag,
                  std::function<void(pugi::xml_document &)> func)
{
    std::ifstream ifs(file);
    std::string line, block;
    bool inBlock = false;
    
    while (std::getline(ifs, line))
    {
        trim(line);
        if (strStartsWith(line, startTag))
            inBlock = true;
        
        if (inBlock)
            block += '\n' + line;
        
        if (line == endTag)
        {
            pugi::xml_document doc;
            doc.load_string(block.c_str());
            func(doc);
            inBlock = false;
            block.clear();
        }
    }
}


struct feature
{
    std::string id;
    numType ret, mz, intensity;
    numType retMin, retMax;
    numType mzMin, mzMax;
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
    
    for (auto hull: featNode.children("convexhull"))
    {
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
            break;
        }
    }
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::DataFrame parseFeatureXMLFile(Rcpp::CharacterVector file)
{
    std::vector<std::string> ids;
    std::vector<numType> ints, rets, mzs, retMins, retMaxs, mzMins, mzMaxs;
    
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
    });
    
    return Rcpp::DataFrame::create(Rcpp::Named("ID") = ids,
                                   Rcpp::Named("ret") = rets,
                                   Rcpp::Named("mz") = mzs,
                                   Rcpp::Named("area") = ints, // intensity reported by OpenMS is actually the area
                                   Rcpp::Named("retmin") = retMins,
                                   Rcpp::Named("retmax") = retMaxs,
                                   Rcpp::Named("mzmin") = mzMins,
                                   Rcpp::Named("mzmax") = mzMaxs,
                                   Rcpp::Named("stringsAsFactors") = false);
}

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
