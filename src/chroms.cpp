#include <fstream>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "cpp-base64/base64.cpp"

// [[Rcpp::export]]
void writeChromsToMzML(Rcpp::List EICs, const std::string &out)
{
    auto indent = [](int level) { return std::string(level * 4, ' '); };

    std::ofstream ofile(out);
    //ofile << std::fixed << std::setprecision(6);

    // generate a minimal acceptable mzML file for FeatureFinderMRM
    ofile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
          << "<mzML>\n"
          << indent(1) << "<run id=\"runID\" defaultInstrumentConfigurationRef=\"IC1\">\n"
          << indent(2) << "<chromatogramList count = \"" << EICs.size() << "\" defaultDataProcessingRef=\"patRoon\">\n";

    for (int i = 0; i < EICs.size(); ++i)
    {
        const Rcpp::List item = EICs[i];
        const auto times = Rcpp::as<std::vector<double>>(item["time"]);
        const auto ints = Rcpp::as<std::vector<double>>(item["intensity"]);
        const auto timesEnc = base64_encode(reinterpret_cast<unsigned const char*>(times.data()), times.size() * sizeof(double));
        const auto intsEnc = base64_encode(reinterpret_cast<unsigned const char*>(ints.data()), ints.size() * sizeof(double));
        const auto dataSize = sizeof(double) * 8;

        ofile << indent(3) << "<chromatogram index=\"" << i << "\" id=\"chrom_" << i << "\" defaultArrayLength=\""
              << times.size() << "\">\n";

        ofile << indent(4) << R"(<cvParam cvRef="MS" accession="MS:1001473" name="selected reaction monitoring chromatogram" value=""/>)" << "\n";
        ofile << indent(4) << R"(<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>)" << "\n";

        ofile << indent(4) << "<product>\n"
              << indent(5) << "<isolationWindow>\n"
              << indent(6) << R"(<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value=")" << i+1 << R"(" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>)" << "\n"
              << indent(5) << "</isolationWindow>\n"
              << indent(4) << "</product>\n";

        ofile << indent(4) << R"(<binaryDataArrayList count="2">)" << "\n"
              << indent(5) << R"(<binaryDataArray encodedLength="1400">)" << "\n"
              << indent(6) << R"(<cvParam cvRef="MS" accession="MS:1000523" name=")" << dataSize << R"(-bit float" value=""/>)" << "\n"
              << indent(6) << R"(<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>)" << "\n"
              << indent(6) << R"(<cvParam cvRef="MS" accession="MS:1000595" name="time array" value="" unitCvRef="UO" unitAccession="UO:0000010" unitName="seconds"/>)" << "\n"
              << indent(6) << "<binary>" << timesEnc << "</binary>\n"
              << indent(5) << "</binaryDataArray>\n"
              << indent(5) << R"(<binaryDataArray encodedLength="1400">)" << "\n"
              << indent(6) << R"(<cvParam cvRef="MS" accession="MS:1000523" name=")" << dataSize << R"(-bit float" value=""/>)" << "\n"
              << indent(6) << R"(<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>)" << "\n"
              << indent(6) << R"(<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>)" << "\n"
              << indent(6) << "<binary>" << intsEnc << "</binary>\n"
              << indent(5) << "</binaryDataArray>\n"
              << indent(4) << "</binaryDataArrayList>\n";

        ofile << indent(3) << "</chromatogram>\n";
    }

    ofile << indent(2) << "</chromatogramList>\n"
          << indent(1) << "</run>\n"
          << "</mzML>\n";
}
