#include <iostream>
#include "external_libraries.hpp"
#include "xml_utils.h"
#include <string>
#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::List rcpp_parse_chromatograms(std::string file_path, Rcpp::IntegerVector index = NA_INTEGER) {

  Rcpp::List list_output;

  const char * path = file_path.c_str();

  pugi::xml_document doc;

  pugi::xml_parse_result result = doc.load_file(path);

  if (result) {
    pugi::xml_node root = doc.document_element();

    if (Rcpp::IntegerVector::is_na(index[0])) {
      return xml_utils::parse_chromatograms(root);

    } else {
      return xml_utils::parse_partial_chromatograms(root, index);
    }

  } else {
    std::cout << "\u2717 XML file could not be opened! " << result.description() << std::endl;
  }

  return list_output;
}
