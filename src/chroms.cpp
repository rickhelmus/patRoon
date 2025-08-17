#include <fstream>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "msdata.h"

#include "cpp-base64/base64.cpp"

namespace {

std::string indentStr(unsigned level)
{
    return std::string(level * 4, ' ');
}

}

// [[Rcpp::export]]
void writeTraML(const std::vector<std::string> &IDs, const std::string &out)
{
    // generate a minimal acceptable TraML file for MRMTransitionGroupPicker
    
    std::ofstream ofile(out);
    
    ofile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
          << R"(<TraML version="1.0.0">)" << "\n";
    
    ofile << indentStr(1) << "<CompoundList>\n";
    for (size_t i=0; i<IDs.size(); ++i)
    {
        ofile << indentStr(2) << "<Peptide id=\"pep_" << i << "\" sequence=\"A\">\n"
              << indentStr(3) << "<RetentionTimeList>\n"
              << indentStr(4) << "<RetentionTime>\n"
              << indentStr(5) << R"(<cvParam cvRef="MS" accession="MS:1002005" name="iRT retention time normalization standard" value="0"/>)" << "\n"
              << indentStr(4) << "</RetentionTime>\n"
              << indentStr(3) << "</RetentionTimeList>\n"
              << indentStr(2) << "</Peptide>\n";
    }
    ofile << indentStr(1) << "</CompoundList>\n";
    
    ofile << indentStr(1) << "<TransitionList>\n";
    for (size_t i=0; i<IDs.size(); ++i)
    {
        ofile << indentStr(2) << "<Transition id=\"" << IDs[i] << "\" peptideRef=\"pep_" << i << "\">\n"
              << indentStr(3) << "<Precursor>\n"
              << indentStr(3) << "</Precursor>\n"
              << indentStr(2) << "</Transition>\n";
    }
    ofile << indentStr(1) << "</TransitionList>\n";

    ofile << "</TraML>\n";    
}

// [[Rcpp::export]]
void writeChromsToMzML(Rcpp::List EICs, bool fillEICs, const std::vector<std::string> &IDs, const std::string &out)
{
    // generate a minimal acceptable mzML file for MRMTransitionGroupPicker
    
    std::ofstream ofile(out);
    //ofile << std::fixed << std::setprecision(6);

    ofile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
          << "<mzML>\n"
          << indentStr(1) << "<run id=\"runID\" defaultInstrumentConfigurationRef=\"IC1\">\n"
          << indentStr(2) << "<chromatogramList count = \"" << EICs.size() << "\" defaultDataProcessingRef=\"patRoon\">\n";

    // BUG: seems files with only double arrays are supported by OpenMS...
    
    const auto dataSizeTime = sizeof(double) * 8;
    const auto dataSizeIntens = sizeof(double) * 8;
    const auto allTimes = (fillEICs) ? Rcpp::as<std::vector<double>>(EICs.attr("allXValues")) : std::vector<double>();
    
    for (int i = 0; i < EICs.size(); ++i)
    {
        const Rcpp::NumericMatrix item = EICs[i];
        const auto times = Rcpp::as<std::vector<double>>(static_cast<Rcpp::NumericVector>(item(Rcpp::_, 0)));
        const auto ints = Rcpp::as<std::vector<double>>(static_cast<Rcpp::NumericVector>(item(Rcpp::_, 1)));
        std::string timesEnc, intsEnc;
        if (fillEICs)
        {
            const auto allInts = fillEIXIntensities(allTimes, times, ints);
            timesEnc = base64_encode(reinterpret_cast<unsigned const char*>(allTimes.data()), allTimes.size() * sizeof(double));
            intsEnc = base64_encode(reinterpret_cast<unsigned const char*>(allInts.data()), allInts.size() * sizeof(double));
        }
        else
        {
            timesEnc = base64_encode(reinterpret_cast<unsigned const char*>(times.data()), times.size() * sizeof(double));
            intsEnc = base64_encode(reinterpret_cast<unsigned const char*>(ints.data()), ints.size() * sizeof(double));
        }
        const auto arrayLen = (fillEICs) ? allTimes.size() : times.size();

        ofile << indentStr(3) << "<chromatogram index=\"" << i << "\" id=\"" << IDs[i] << "\" defaultArrayLength=\""
              << arrayLen << "\">\n";

        ofile << indentStr(4) << R"(<cvParam cvRef="MS" accession="MS:1001473" name="selected reaction monitoring chromatogram" value=""/>)" << "\n";
        ofile << indentStr(4) << R"(<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>)" << "\n";

        ofile << indentStr(4) << "<product>\n"
              << indentStr(5) << "<isolationWindow>\n"
              << indentStr(6) << R"(<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value=")" << i+1 << R"(" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>)" << "\n"
              << indentStr(5) << "</isolationWindow>\n"
              << indentStr(4) << "</product>\n";

        ofile << indentStr(4) << R"(<binaryDataArrayList count="2">)" << "\n"
              << indentStr(5) << R"(<binaryDataArray encodedLength="1400">)" << "\n"
              << indentStr(6) << R"(<cvParam cvRef="MS" accession="MS:1000523" name=")" << dataSizeTime << R"(-bit float" value=""/>)" << "\n"
              << indentStr(6) << R"(<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>)" << "\n"
              << indentStr(6) << R"(<cvParam cvRef="MS" accession="MS:1000595" name="time array" value="" unitCvRef="UO" unitAccession="UO:0000010" unitName="seconds"/>)" << "\n"
              << indentStr(6) << "<binary>" << timesEnc << "</binary>\n"
              << indentStr(5) << "</binaryDataArray>\n"
              << indentStr(5) << R"(<binaryDataArray encodedLength="1400">)" << "\n"
              << indentStr(6) << R"(<cvParam cvRef="MS" accession="MS:1000523" name=")" << dataSizeIntens << R"(-bit float" value=""/>)" << "\n"
              << indentStr(6) << R"(<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>)" << "\n"
              << indentStr(6) << R"(<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>)" << "\n"
              << indentStr(6) << "<binary>" << intsEnc << "</binary>\n"
              << indentStr(5) << "</binaryDataArray>\n"
              << indentStr(4) << "</binaryDataArrayList>\n";

        ofile << indentStr(3) << "</chromatogram>\n";
    }

    ofile << indentStr(2) << "</chromatogramList>\n"
          << indentStr(1) << "</run>\n"
          << "</mzML>\n";
}
