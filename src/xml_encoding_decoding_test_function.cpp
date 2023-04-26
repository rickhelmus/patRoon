
#include "xml_utils.h"
#include <Rcpp.h>

Rcpp::List xml_utils::encoding_decoding_test_function(Rcpp::NumericVector input) {

  Rcpp::List list_out;

  list_out["input"] = input;


  // encode binary

  std::string little_endian =  xml_utils::encode_little_endian(input);

  list_out["little_endian"] = little_endian;

  std::string big_endian =  xml_utils::encode_big_endian(input);

  list_out["big_endian"] = big_endian;


  // compress with zlib

  std::string little_endian_compressed = xml_utils::compress_zlib(little_endian);

  list_out["little_endian_compressed"] = little_endian_compressed;

  std::string big_endian_compressed = xml_utils::compress_zlib(big_endian);

  list_out["big_endian_compressed"] = big_endian_compressed;


  // encode to base 64

  std::string little_endian_base64 = xml_utils::encode_base64(little_endian_compressed);

  list_out["little_endian_base64"] = little_endian_base64;

  std::string big_endian_base64 = xml_utils::encode_base64(big_endian_compressed);

  list_out["big_endian_base64"] = big_endian_base64;


  // decode to base 64

  std::string little_endian_base64_back = xml_utils::decode_base64(little_endian_base64);

  list_out["little_endian_base64_back"] = little_endian_base64_back;

  std::string big_endian_base64_back = xml_utils::decode_base64(big_endian_base64);

  list_out["big_endian_base64_back"] = big_endian_base64_back;


 // decompress with zlib

  std::string little_endian_compressed_back =  xml_utils::decompress_zlib(little_endian_base64_back);

  list_out["little_endian_compressed_back"] = little_endian_compressed_back;

  std::string big_endian_compressed_back =  xml_utils::decompress_zlib(big_endian_base64_back);

  list_out["big_endian_compressed_back"] = big_endian_compressed_back;


  // decode binary

  list_out["little_endian_back"] = xml_utils::decode_little_endian(little_endian_compressed_back, 8);

  list_out["big_endian_back"] = xml_utils::decode_big_endian(big_endian_compressed_back, 8);

  Rcpp::NumericVector l_back = list_out["little_endian_back"];

  Rcpp::NumericVector b_back = list_out["big_endian_back"];

  Rcpp::LogicalVector l_check = l_back == input;

  Rcpp::LogicalVector b_check = b_back == input;

  list_out["little_endian_check"] = l_check;

  list_out["big_endian_check"] = b_check;


  return list_out;
}
