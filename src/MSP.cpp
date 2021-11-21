#include <fstream>
#include <set>
#include <string>
#include <unordered_map>

#include <Rcpp.h>

#include "utils.h"

struct MSPSpectrum
{
    std::vector<double> mzs, intensities;
};

struct MSPRecord
{
    std::unordered_map<std::string, std::string> values;
    MSPSpectrum spectrum;
};

// [[Rcpp::export]]
Rcpp::List readMSP(Rcpp::CharacterVector file)
{
    std::ifstream fs;
    std::vector<MSPRecord> records;
    std::set<std::string> keys;
    
    fs.open(Rcpp::as<const char *>(file));
    if (fs.is_open())
    {
        std::string line;
        MSPRecord curRec;
        
        Rcpp::Rcout << "Parsing file...";
        
        while (std::getline(fs, line))
        {
            auto cPos = line.find(":");
            if (cPos != std::string::npos)
            {
                std::string key = line.substr(0, cPos);
                std::string val = line.substr(cPos + 1);
                trim(key); trim(val);
                
                if (key == "Num Peaks")
                {
                    for (int n = std::stoi(val); n; --n)
                    {
                        // UNDONE: MSP also allows other formats than one space separated pair per line
                        
                        double m, i;
                        fs >> m >> i;
                        curRec.spectrum.mzs.push_back(m);
                        curRec.spectrum.intensities.push_back(i);
                    }
                    records.push_back(curRec);
                    curRec = MSPRecord();
                }
                else
                {
                    if (key == "DB#")
                        key = "DB_ID"; // rename for R name compat
                    
                    keys.insert(key);
                    
                    auto p = curRec.values.insert({key, val});
                    if (!p.second) // NOT inserted, i.e. already present?
                        curRec.values[key] = p.first->second + ";" + val;
                }
            }
        }
    }
    Rcpp::Rcout << " Done!" << std::endl << "Converting to R data...";
    
    Rcpp::List recordsList(keys.size());
    recordsList.names() = Rcpp::wrap(std::vector<std::string>(keys.begin(), keys.end()));
    for (const std::string &k : keys)
    {
        std::vector<std::string> vals;
        for (const auto &r : records)
        {
            const auto it = r.values.find(k);
            vals.push_back((it == r.values.end()) ? "NA" : it->second);
        }
        recordsList[k] = vals;
    }

    Rcpp::List specList(records.size());
    specList.names() = recordsList["DB_ID"];
    for (size_t i=0; i<specList.size(); ++i)
    {
        Rcpp::NumericMatrix nm(records[i].spectrum.mzs.size(), 2);
        nm(Rcpp::_, 0) = Rcpp::NumericVector(records[i].spectrum.mzs.begin(), records[i].spectrum.mzs.end());
        nm(Rcpp::_, 1) = Rcpp::NumericVector(records[i].spectrum.intensities.begin(),
                                             records[i].spectrum.intensities.end());
        Rcpp::colnames(nm) = Rcpp::CharacterVector({"mz", "intensity"});
        specList[i] = nm;
    }
    
    Rcpp::Rcout << " Done!" << std::endl;
    
    return Rcpp::List::create(Rcpp::Named("records") = Rcpp::DataFrame(recordsList),
                              Rcpp::Named("spectra") = specList);
}
