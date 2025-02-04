/*
 * SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

// [[Rcpp::depends(rapidjsonr)]]

#include <cctype>
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
    std::vector<std::string> annotations;
};

struct MSLibRecord
{
    std::unordered_map<std::string, std::string> values;
    MSLibSpectrum spectrum;
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

void fixAnnFormula(std::string &form)
{
    trim(form);
    
    // get rid of trailing charges
    const auto len = form.length();
    if (len > 0 && (form[len-1] == '+' || form[len-1] == '-'))
        form.pop_back();
    if (hasWS(form))
        form.clear(); // formulas with whitespace are invalid and cleared out
}

Rcpp::List convertRecordsToRData(const std::vector<MSLibRecord> &records, const std::vector<std::string> &keys)
{
    if (records.empty())
        return Rcpp::List::create(Rcpp::Named("records") = Rcpp::DataFrame(),
                                  Rcpp::Named("spectraMZs") = Rcpp::List(),
                                  Rcpp::Named("spectraInts") = Rcpp::List(),
                                  Rcpp::Named("annotations") = Rcpp::List());
    
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
    
    Rcpp::List specListMZs(records.size()), specListInts(records.size()), annList(records.size());
    specListMZs.names() = specListInts.names() = annList.names() = recordsList["DB_ID"];
    for (auto i=0; i<records.size(); ++i)
    {
        specListMZs[i] = records[i].spectrum.mzs;
        specListInts[i] = records[i].spectrum.intensities;
        annList[i] = records[i].spectrum.annotations;
    }
    
    return Rcpp::List::create(Rcpp::Named("records") = Rcpp::DataFrame(recordsList),
                              Rcpp::Named("spectraMZs") = specListMZs,
                              Rcpp::Named("spectraInts") = specListInts,
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
    auto parseComment = [&](const std::string &comment, const std::string &outVar, MSLibRecord &record,
                            const std::string &var, const std::string &altVar = "")
    {
        std::string cv;
        if (record.values.find(outVar) == record.values.end() &&
            (parseMSPComments(comment, var, cv) || (!altVar.empty() && parseMSPComments(comment, altVar, cv))))
        {
            addKey(outVar);
            record.values[outVar] = cv;
        }
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
                
                auto keyLower = key;
                std::transform(keyLower.begin(), keyLower.end(), keyLower.begin(), ::tolower);
                
                if (keyLower == "num peaks")
                {
                    bool getMZ = true;
                    const char *separators = " \t,;:()[]{}"; // separators from NIST manual
                    int numPeaks = std::stoi(val);
                    while (numPeaks && std::getline(fs, line))
                    {
                        const auto lineLen = line.length();
                        std::string::size_type startPos = 0, endPos = 0;
                        do
                        {
                            endPos = line.find_first_of(separators, startPos);
                            const std::string s = line.substr(startPos, endPos - startPos);
                            startPos = endPos + 1;
                            
                            if (s[0] == '"')
                                break; // annotation string --> not supported yet, they seem to be at the end of the line so just skip everything

                            const double val = stod(s);
                            if (getMZ)
                                curRec.spectrum.mzs.push_back(val);
                            else
                            {
                                curRec.spectrum.intensities.push_back(val);
                                --numPeaks;
                            }
                            getMZ = !getMZ;
                        } while (endPos != std::string::npos && (startPos + 1) <= lineLen);
                        // UNDONE: sort peaks?
                    }

                    // NOTE: Num Peaks is always the last entry, finish up record
                    
                    // Parse comments?
                    if (pComments && curRec.values.find("Comments") != curRec.values.end())
                    {
                        const std::string com = curRec.values["Comments"];
                        parseComment(com, "SMILES", curRec, "SMILES", "computed SMILES");
                        parseComment(com, "InChI", curRec, "InChI", "computed InChI");
                        // NOTE: MoNA saves Splash as uppercase in comments
                        parseComment(com, "SPLASH", curRec, "SPLASH", "Splash");
                        parseComment(com, "CAS", curRec, "cas");
                        parseComment(com, "PubChemCID", curRec, "pubchem cid");
                        parseComment(com, "ChemSpiderID", curRec, "chemspider");
                        parseComment(com, "Ionization", curRec, "ionization");
                        parseComment(com, "Resolution", curRec, "resolution");
                    }
                    
                    // ensure there is a DB_ID
                    if (curRec.values.count("DB_ID") == 0)
                    {
                        addKey("DB_ID");
                        curRec.values.insert({std::string("DB_ID"), std::string("ID") + std::to_string(records.size() + 1)});
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
                    else if (key == "Synon" && val == "$:00in-source")
                        continue; // Skip these weird markers
                    
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
            const Rcpp::DataFrame spec = Rcpp::as<Rcpp::DataFrame>(spectraList[row]);
            outf << "Num Peaks: " << spec.nrow() << "\n";
            const std::vector<double> mzs = spec[0], ints = spec[1];
            for (int srow=0; srow<spec.nrow(); ++srow)
                outf << mzs[srow] << " " << ints[srow] << "\n";
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
        { "cas", "CAS" },
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
                    // use insert: some values like SMILES etc may be specified double. The MSP files seem to contain
                    // the first value, so stick to that for consistency.
                    const auto p = record.values.insert({mappedVal->second, itv->value.GetString()});
                    if (p.second) // TRUE if inserted
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
                Rcpp::stop("JSON parse error (offset: %u): %s", (unsigned)d.GetErrorOffset(),
                           rapidjson::GetParseError_En(d.GetParseError()));
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
                
                rapidjson::Value::ConstMemberIterator cit2 = cit->value[0].FindMember("metaData");
                if (cit2 != cit->value[0].MemberEnd() && cit2->value.IsArray())
                {
                    for (auto &item : cit2->value.GetArray())
                        getStringMD(item, curRec, JSONCompMDMapping);
                }
                
                cit2 = cit->value[0].FindMember("names");
                if (cit2 != cit->value[0].MemberEnd() && cit2->value.IsArray())
                {
                    for (auto &item : cit2->value.GetArray())
                    {
                        rapidjson::Value::ConstMemberIterator itn = item.FindMember("name");
                        if (itn != item.MemberEnd() && itn->value.IsString())
                        {
                            const std::string n = itn->value.GetString();
                            auto p = curRec.values.insert({"Name", n});
                            if (!p.second) // NOT inserted, i.e. already present?
                            {
                                // Try the same for Synon, paste if already present
                                p = curRec.values.insert({"Synon", n});
                                if (!p.second)
                                    curRec.values["Synon"] = p.first->second + ";" + n;
                                else
                                    addKey("Synon");
                            }
                            else
                                addKey("Name");
                        }
                    }
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
                curRec.spectrum.annotations.resize(curRec.spectrum.mzs.size());
                for (auto &item : cit->value.GetArray())
                {
                    rapidjson::Value::ConstMemberIterator itf = item.FindMember("name"), itm = item.FindMember("value");
                    if (itf != item.MemberEnd() && itf->value.IsString() &&
                        itm != item.MemberEnd() && itm->value.IsString())
                    {
                        const double mz = std::stod(itm->value.GetString());
                        for (size_t i=0; i<curRec.spectrum.mzs.size(); ++i)
                        {
                            if (compareTol(mz, curRec.spectrum.mzs[i], 1E-6))
                            {
                                std::string form = itf->value.GetString();
                                fixAnnFormula(form);
                                curRec.spectrum.annotations[i] = form;
                                break;
                            }
                        }
                    }
                }
            }
            
            cit = d.FindMember("splash");
            if (cit != d.MemberEnd() && cit->value.IsObject())
            {
                curRec.values["SPLASH"] = cit->value["splash"].GetString();
                addKey("SPLASH");
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
