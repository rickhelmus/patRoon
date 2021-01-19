#include <iomanip>
#include <fstream>
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

// generating XML via package is too slow...http://r.789695.n4.nabble.com/Creating-XML-document-extremely-slow-td4376088.html
// generate by simply writing text to file instead
// [[Rcpp::export]]
void writeFeatureXML(Rcpp::DataFrame featList, Rcpp::CharacterVector out, Rcpp::LogicalVector hulls)
{
    const char *outStr = Rcpp::as<const char *>(out);
    const Rcpp::NumericVector rets = featList["ret"];
    const Rcpp::NumericVector mzs = featList["mz"];
    const Rcpp::NumericVector areas = featList["area"];
    const Rcpp::NumericVector retmins = featList["retmin"];
    const Rcpp::NumericVector retmaxs = featList["retmax"];
    const Rcpp::NumericVector mzmins = featList["mzmin"];
    const Rcpp::NumericVector mzmaxs = featList["mzmax"];
    const int ftCount = rets.length();
    const bool doHulls = Rcpp::as<bool>(hulls);
    
    std::ofstream ofile(outStr);
    ofile << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<featureMap version=\"1.9\" id=\"fm\" "
          << "xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeatureXML_1_9.xsd\" "
          << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    
    auto indent = [](int level) { return(std::string(level * 4, ' ')); };
    
    ofile << indent(1) << "<featureList count=\"" << ftCount << "\">\n";
    ofile << std::fixed << std::setprecision(6);
    for (int fti=0; fti<ftCount; ++fti)
    {
        ofile << indent(2) << "<feature id=\"f_" << fti + 1 << "\">\n";
        
        ofile << indent(3) << "<position dim=\"0\">" << rets[fti] << "</position>\n";
        ofile << indent(3) << "<position dim=\"1\">" << mzs[fti] << "</position>\n";
        ofile << indent(3) << "<intensity>" << areas[fti] << "</intensity>\n";
        ofile << indent(3) << "<quality dim=\"0\">0</quality>\n";
        ofile << indent(3) << "<quality dim=\"1\">0</quality>\n";
        ofile << indent(3) << "<overallquality>0</overallquality>\n";
        ofile << indent(3) << "<charge>0</charge>\n";
        
        if (doHulls)
        {
            // NOTE: for now just include min/max RT and m/z, as this is sufficient for
            // MetaboliteAdductDecharger...
            ofile << indent(3) << "<convexhull nr=\"0\">\n";
            ofile << indent(4) << "<pt x=\"" << retmins[fti] << "\" y=\"" << mzmins[fti] << "\" />\n";
            ofile << indent(4) << "<pt x=\"" << retmaxs[fti] << "\" y=\"" << mzmaxs[fti] << "\" />\n";
            ofile << indent(3) << "</convexhull>\n";
        }
        
        ofile << indent(2) << "</feature>\n";
    }
    
    ofile << indent(1) << "</featureList>\n";
    ofile << "</featureMap>\n";
    
    ofile.close();
}
