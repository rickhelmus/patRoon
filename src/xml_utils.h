#ifndef XML_UTILS_H
#define XML_UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include "external_libraries.hpp"
#include <Rcpp.h>

namespace xml_utils {

  // // encode tools // //

  std::string encode_little_endian(const Rcpp::NumericVector& input);

  std::string encode_big_endian(const Rcpp::NumericVector& input);

  std::string compress_zlib(const std::string& str);

  std::string encode_base64(const std::string& str);

  // // decode tools // //

  std::string decode_base64(const std::string& encoded_string);

  std::string decompress_zlib(const std::string& compressed_string);

  Rcpp::NumericVector decode_little_endian(const std::string& str, const int& precision);

  Rcpp::NumericVector decode_big_endian(const std::string& str, const int& precision);

  // // structures // //

  struct spectraHeaders {
    std::vector<int> index;
    std::vector<int> scan;
    std::vector<int> traces;
    std::vector<std::string> polarity;
    std::vector<double> bpcmz;
    std::vector<double> bpcint;
    std::vector<double> ticint;
    std::vector<int> level;
    std::vector<std::string> mode;
    std::vector<double> mzlow;
    std::vector<double> mzhigh;
    std::vector<double> rt;
    std::vector<double> drift;
    std::vector<int> pre_scan;
    std::vector<double> pre_mz;
    std::vector<double> pre_loweroffset;
    std::vector<double> pre_upperoffset;
    std::vector<double> pre_ce;
  }; // spectraHeaders

  struct chromatogramsHeaders {
    std::vector<int> index;
    std::vector<std::string> id;
    std::vector<int> traces;
    std::vector<std::string> polarity;
    std::vector<double> pre_mz;
    std::vector<double> pre_ce;
    std::vector<double> pro_mz;
  }; // spectraHeaders

  struct runSummary {
    int spectra_number;
    Rcpp::CharacterVector polarity;
    Rcpp::CharacterVector mode;
    std::vector<int> levels;
    double mz_low;
    double mz_high;
    double rt_start;
    double rt_end;
    bool has_ion_mobility;
    Rcpp::List spectra;
    int chromatograms_number;
    Rcpp::List chromatograms;
  }; // runSummary

  struct spectraHeadersOriginal {
    std::string fileFormat;
    std::string run_id;
    std::string run_defaultInstrumentConfigurationRef;
    std::string run_startTimeStamp;
    std::string run_defaultSourceFileRef;
    int specList_count;
    std::string specList_defaultDataProcessingRef;
    std::vector<int> spec_index;
    std::vector<int> spec_id;
    std::vector<int> scan;
    std::vector<int> spec_defaultArrayLength;
    std::vector<std::string> spec_polarity;
    std::vector<double> spec_bpcmz;
    std::vector<double> spec_bpcint;
    std::vector<double> spec_ticint;
    std::vector<int> spec_level;
    std::vector<std::string> spec_mode;
    std::vector<double> spec_mzlow;
    std::vector<double> spec_mzhigh;
    std::string spec_title;
    std::vector<double> scan_rt;
    std::vector<double> scan_drift;
    std::string scan_filter_string;
    std::vector<double> scan_injectionIonTime;
    std::vector<int> pre_scan;
    std::vector<double> pre_mz;
    std::vector<double> pre_loweroffset;
    std::vector<double> pre_upperoffset;
    std::vector<double> ion_mz;
    std::vector<double> ion_charge;
    std::vector<double> ion_intensity;
    std::vector<double> activation_type;
    std::vector<double> activation_ce;
  }; // spectraHeadersOriginal


  // // parsing tools // //

  std::vector<int> mzml_get_precision(pugi::xml_node& node);
  int mzxml_get_precision(pugi::xml_node& node);

  std::vector<std::string> mzml_get_compression(pugi::xml_node& node);
  std::string mzxml_get_compression(pugi::xml_node& node);

  Rcpp::CharacterVector mzml_get_binary_type(pugi::xml_node& node);

  Rcpp::CharacterVector mzml_get_binary_names(pugi::xml_node& node);

  Rcpp::NumericMatrix mzml_parse_binary_data_from_spectrum_node(
      const pugi::xml_node& node,
      const std::vector<int>& precision,
      const std::vector<std::string>& compression,
      const Rcpp::CharacterVector& cols
  );

  Rcpp::NumericMatrix mzxml_parse_binary_data_from_spectrum_node(
      const pugi::xml_node& node,
      const int& precision,
      const std::string& compression,
      const Rcpp::CharacterVector& cols
  );

  Rcpp::List parse_spectra(const pugi::xml_node& root);

  Rcpp::List parse_partial_spectra(
      const pugi::xml_node& root,
      Rcpp::IntegerVector& index
  );

  Rcpp::List parse_chromatograms(const pugi::xml_node& root);

  Rcpp::List parse_partial_chromatograms(
      const pugi::xml_node& root,
      Rcpp::IntegerVector& index
  );

  std::list<std::vector<std::string>> mzml_instrument_parser(const pugi::xml_node& node_mzml);
  std::list<std::vector<std::string>> mzxml_instrument_parser(const pugi::xml_node& node_mzxml);
  Rcpp::List parse_instrument(const pugi::xml_node& root);

  std::list<std::vector<std::string>> mzml_software_parser(const pugi::xml_node& node_mzml);
  std::list<std::vector<std::string>> mzxml_software_parser(const pugi::xml_node& node_mzxml);
  Rcpp::List parse_software(const pugi::xml_node& root);

  void mzml_spectra_headers_parser(const pugi::xml_node& node_mzml, spectraHeaders& output);
  void mzxml_spectra_headers_parser(const pugi::xml_node& node_mzxml, spectraHeaders& output);
  spectraHeaders parse_spectra_headers(const pugi::xml_node& root);
  Rcpp::List spectraHeaders_to_list(const spectraHeaders& headers_cpp);
  spectraHeaders list_to_spectraHeaders(const Rcpp::List& run);

  void mzml_chromatograms_headers_parser(const pugi::xml_node& node_mzml, chromatogramsHeaders& output);
  chromatogramsHeaders parse_chromatograms_headers(const pugi::xml_node& root);
  Rcpp::List chromatogramsHeaders_to_list(const chromatogramsHeaders& headers_cpp);

  runSummary run_summary(spectraHeaders& spec_headers, chromatogramsHeaders& chrom_headers);

  // // other functions // //

  // in file xml_encoding_decoding_test_function.cpp
  Rcpp::List encoding_decoding_test_function(Rcpp::NumericVector input);

} // xml_utils

#endif
