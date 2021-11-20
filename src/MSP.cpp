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
    
    Rcpp::List recordsList(keys.size());
    recordsList.names() = Rcpp::wrap(std::vector<std::string>(keys.begin(), keys.end()));
    for (std::string k : keys)
    {
        std::vector<std::string> vals;
        for (const auto r : records)
        {
            const auto it = r.values.find(k);
            vals.push_back((it == r.values.end()) ? "NA" : it->second);
        }
        recordsList[k] = vals;
    }

    Rcpp::List specList(records.size());
    specList.names() = recordsList["Name"];
    for (size_t i=0; i<specList.size(); ++i)
        specList[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = records[i].spectrum.mzs,
                                              Rcpp::Named("intensity") = records[i].spectrum.intensities);
    
    return Rcpp::List::create(Rcpp::Named("records") = Rcpp::DataFrame(recordsList),
                              Rcpp::Named("spectra") = specList);
}
