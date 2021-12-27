// [[Rcpp::depends(rapidjsonr)]]

#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <unordered_map>

#include <Rcpp.h>

#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "utils.h"

namespace {

struct MSLibSpectrum
{
    std::vector<double> mzs, intensities;
};

struct MSLibAnnotations
{
    std::vector<double> mzs;
    std::vector<std::string> formulas;
};
    
struct MSLibRecord
{
    std::unordered_map<std::string, std::string> values;
    MSLibSpectrum spectrum;
    MSLibAnnotations annotations;
};

bool parseMSPComments(const std::string &comments, const std::string &field, std::string &out)
{
    const std::string toMatch = '"' + field + '=';
    auto start = comments.find(toMatch);
    if (start != std::string::npos)
    {
        start += toMatch.length();
        const auto end = comments.find('"', start);
        if (end != std::string::npos)
        {
            out = comments.substr(start, end-start);
            return true;
        }
    }
    return false;
}

Rcpp::List convertRecordsToRData(const std::vector<MSLibRecord> &records, const std::vector<std::string> &keys)
{
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
    for (int i=0; i<specList.size(); ++i)
    {
        Rcpp::NumericMatrix nm(records[i].spectrum.mzs.size(), 2);
        nm(Rcpp::_, 0) = Rcpp::NumericVector(records[i].spectrum.mzs.begin(), records[i].spectrum.mzs.end());
        nm(Rcpp::_, 1) = Rcpp::NumericVector(records[i].spectrum.intensities.begin(),
           records[i].spectrum.intensities.end());
        Rcpp::colnames(nm) = Rcpp::CharacterVector({"mz", "intensity"});
        specList[i] = nm;
    }
    
    Rcpp::List annList(records.size());
    annList.names() = recordsList["DB_ID"];
    for (int i=0; i<annList.size(); ++i)
    {
        annList[i] = Rcpp::DataFrame::create(Rcpp::Named("mz") = records[i].annotations.mzs,
                                             Rcpp::Named("formula") = records[i].annotations.formulas);
    }
    
    return Rcpp::List::create(Rcpp::Named("records") = Rcpp::DataFrame(recordsList),
                              Rcpp::Named("spectra") = specList,
                              Rcpp::Named("annotations") = annList);
}

}


// [[Rcpp::export]]
Rcpp::List readMSP(Rcpp::CharacterVector file, Rcpp::LogicalVector pc)
{
    const bool pComments = Rcpp::as<bool>(pc);
    std::ifstream fs;
    std::vector<MSLibRecord> records;
    std::vector<std::string> keys;
    
    auto addKey = [&keys](const std::string &k)
    {
        if (std::find(keys.begin(), keys.end(), k) == keys.end())
            keys.push_back(k);
    };
    
    fs.open(Rcpp::as<const char *>(file));
    if (fs.is_open())
    {
        std::string line;
        MSLibRecord curRec;
        
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
                        std::string m, i;
                        fs >> m >> i;
                        
                        curRec.spectrum.mzs.push_back(stod(m));
                        curRec.spectrum.intensities.push_back(stod(i));
                    }
                    
                    // NOTE: Num Peaks is always the last entry, finish up record
                    
                    // Parse comments?
                    if (pComments && curRec.values.find("Comments") != curRec.values.end())
                    {
                        const std::string com = curRec.values["Comments"];
                        std::string cv;
                        if (curRec.values.find("SMILES") == curRec.values.end() &&
                            (parseMSPComments(com, "SMILES", cv) || parseMSPComments(com, "computed SMILES", cv)))
                        {
                            addKey("SMILES");
                            curRec.values["SMILES"] = cv;
                        }
                        if (curRec.values.find("InChI") == curRec.values.end() &&
                            (parseMSPComments(com, "InChI", cv) || parseMSPComments(com, "computed InChI", cv)))
                        {
                            addKey("InChI");
                            curRec.values["InChI"] = cv;
                        }
                        // NOTE: MoNA saves Splash as uppercase in comments
                        if ((curRec.values.find("Splash") == curRec.values.end() ||
                             curRec.values.find("SPLASH") == curRec.values.end()) && parseMSPComments(com, "SPLASH", cv))
                        {
                            addKey("SPLASH");
                            curRec.values["SPLASH"] = cv;
                        }
                    }
                    
                    records.push_back(curRec);
                    curRec = MSLibRecord();
                    
                    if ((records.size() % 25000) == 0)
                        Rcpp::Rcout << "Read " << records.size() << " records\n";
                }
                else
                {
                    if (key == "DB#")
                        key = "DB_ID"; // rename for R name compat
                    
                    addKey(key);
                    
                    auto p = curRec.values.insert({key, val});
                    if (!p.second) // NOT inserted, i.e. already present?
                        curRec.values[key] = p.first->second + ";" + val;
                }
            }
        }
        fs.close();
        Rcpp::Rcout << "Read " << records.size() << " records\n";
    }
    
    return convertRecordsToRData(records, keys);
}

// [[Rcpp::export]]
void writeMSPLibrary(Rcpp::CharacterMatrix recordsM, Rcpp::List spectraList, Rcpp::CharacterVector outCV)
{
    // UNDONE: numeric precision seems to reduce
    const char *out = Rcpp::as<const char *>(outCV);
    std::ofstream outf(out);
    outf << std::fixed << std::setprecision(6);
    const Rcpp::CharacterVector fields = Rcpp::colnames(recordsM);
    if (outf.is_open())
    {
        for (int row=0; row<recordsM.nrow(); ++row)
        {
            for (int col=0; col<recordsM.ncol(); ++col)
            {
                const char *f = fields[col], *v = recordsM(row, col);
                outf << ((!strcmp(f, "DB_ID")) ? "DB#" : f) << ": " << v << "\n";
            }
            const Rcpp::NumericMatrix spec = spectraList[row];
            outf << "Num Peaks: " << spec.nrow() << "\n";
            for (int srow=0; srow<spec.nrow(); ++srow)
                outf << spec(srow, 0) << " " << spec(srow, 1) << "\n";
            outf << "\n";
            
            if (row > 0 && (row == (recordsM.nrow()-1) || (row % 25000) == 0))
                Rcpp::Rcout << "Wrote " << row << " records\n";
        }
        outf.close();
    }
}

// [[Rcpp::export]]
Rcpp::List readMoNAJSON(Rcpp::CharacterVector file)
{
    const std::unordered_map<std::string, std::string> JSONCompMDMapping =  {
        { "molecular formula", "Formula" },
        { "SMILES", "SMILES" },
        { "InChI", "InChI" },
        { "InChIKey", "InChIKey" },
        { "pubchem cid", "PubChemCID" },
        { "chemspider", "ChemSpiderID" },
        { "total exact mass", "ExactMass" }
    };
    const std::unordered_map<std::string, std::string> JSONRecMDMapping =  {
        { "instrument", "Instrument" },
        { "instrument type", "Instrument_type" },
        { "ms level", "Spectrum_type" },
        { "ionization", "Ionization" },
        { "Fragmentation", "Fragmentation" },
        { "collision energy", "Collision_energy" },
        { "resolution", "Resolution" },
        { "ionization mode", "Ion_mode" },
        { "precursor m/z", "PrecursorMZ" },
        { "precursor type", "Precursor_type" }
    };
    
    std::ifstream fs;
    std::vector<MSLibRecord> records;
    std::vector<std::string> keys;
    
    auto addKey = [&keys](const std::string &k)
    {
        if (std::find(keys.begin(), keys.end(), k) == keys.end())
            keys.push_back(k);
    };
    
    const auto getString = [&](const rapidjson::Value &val, const char *var, MSLibRecord &record,
                               const char *outVar = NULL)
    {
        rapidjson::Value::ConstMemberIterator it = val.FindMember(var);
        if (it != val.MemberEnd() && it->value.IsString())
        {
            if (outVar == NULL)
                outVar = var;
            record.values[outVar] = it->value.GetString();
            addKey(outVar);
            return true;
        }
        return false;
    };
    const auto getStringMD = [&](const rapidjson::Value &val, MSLibRecord &record,
                                 const std::unordered_map<std::string, std::string> &mapping)
    {
        rapidjson::Value::ConstMemberIterator itn = val.FindMember("name");
        if (itn != val.MemberEnd() && itn->value.IsString())
        {
            const auto mappedVal = mapping.find(itn->value.GetString());
            if (mappedVal != mapping.end())
            {
                rapidjson::Value::ConstMemberIterator itv = val.FindMember("value");
                if (itv != val.MemberEnd() && itv->value.IsString())
                {
                    record.values[mappedVal->second] = itv->value.GetString();
                    addKey(mappedVal->second);
                }
            }
            return true;
        }
        return false;
    };
    
        
    fs.open(Rcpp::as<const char *>(file));
    if (fs.is_open())
    {
        std::string line;
        while (std::getline(fs, line))
        {
            if (line.empty() || line[0] == '[')
                continue;
            if (line[0] == ']')
                break;
            const auto len = line.length();
            if (line[len-1] == ',')
                line.pop_back();
            rapidjson::Document d;
            //d.ParseInsitu<rapidjson::kParseNumbersAsStringsFlag>(&line[0]); // see https://stackoverflow.com/a/4152881
            // BUG: kParseNumbersAsStringsFlag doesn't seem to work with Insitu parsing
            d.Parse<rapidjson::kParseNumbersAsStringsFlag>(line.c_str());
            if (d.HasParseError())
            {
                // UNDONE
                /*fprintf(stderr, "\nError(offset %u): %s\n",
                        (unsigned)d.GetErrorOffset(),
                        GetParseError_En(d.GetParseError()));*/
            }
            
            MSLibRecord curRec;

            // UNDONE: something better than 'continue'? or throw error in getString?
            
            if (!getString(d, "id", curRec, "DB_ID")) continue;
            
            rapidjson::Value::ConstMemberIterator cit = d.FindMember("compound");
            if (cit != d.MemberEnd() && cit->value.IsArray() && !cit->value.Empty())
            {
                // NOTE: these will also be taken from the metaData if present
                getString(cit->value[0], "inchi", curRec, "InChI");
                getString(cit->value[0], "inchiKey", curRec, "InChIKey");
                
                rapidjson::Value::ConstMemberIterator mit = cit->value[0].FindMember("metaData");
                if (mit != cit->value[0].MemberEnd() && mit->value.IsArray())
                {
                    for (auto &item : mit->value.GetArray())
                        getStringMD(item, curRec, JSONCompMDMapping);
                }
            }
            
            cit = d.FindMember("metaData");
            if (cit != d.MemberEnd() && cit->value.IsArray())
            {
                // NOTE: Some metaData might be present in duplicate (e.g. given/calculated SMILES)
                // For now the last value is always kept
                for (auto &item : cit->value.GetArray())
                    getStringMD(item, curRec, JSONRecMDMapping);
            }
            
            cit = d.FindMember("score");
            if (cit != d.MemberEnd() && cit->value.IsObject())
                getString(cit->value, "score", curRec, "Score");
            
            cit = d.FindMember("spectrum");
            if (cit != d.MemberEnd() && cit->value.IsString())
            {
                // m/z / intensity pair is separated by colon, pairs are separated by space
                std::stringstream strm(cit->value.GetString());
                double mz, inten;
                char colon; // dummy
                while(strm >> mz >> colon >> inten)
                {
                    curRec.spectrum.mzs.push_back(mz);
                    curRec.spectrum.intensities.push_back(inten);
                }
            }

            cit = d.FindMember("annotations");
            if (cit != d.MemberEnd() && cit->value.IsArray())
            {
                for (auto &item : cit->value.GetArray())
                {
                    rapidjson::Value::ConstMemberIterator itf = item.FindMember("name"), itm = item.FindMember("value");
                    if (itf != item.MemberEnd() && itf->value.IsString() &&
                        itm != item.MemberEnd() && itm->value.IsString())
                    {
                        curRec.annotations.mzs.push_back(std::stod(itm->value.GetString()));
                        curRec.annotations.formulas.push_back(itf->value.GetString());
                    }
                }
            }
            
            records.push_back(curRec);
            
            if ((records.size() % 25000) == 0)
                Rcpp::Rcout << "Read " << records.size() << " records\n";
        }
        
        fs.close();
        Rcpp::Rcout << "Read " << records.size() << " records\n";
    }
    
    return convertRecordsToRData(records, keys);
}
