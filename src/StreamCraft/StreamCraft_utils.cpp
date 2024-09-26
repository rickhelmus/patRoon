#include "StreamCraft_utils.hpp"
#include <vector>
#include <string>
#include <cstring>
#include <cstdint>
#include <stdexcept>
#include <zlib.h>

std::string sc::utils::encode_little_endian(const std::vector<double>& input, const int& precision) {

  if (precision == 8) {
    std::vector<uint8_t> bytes(sizeof(double) * input.size());
    std::memcpy(bytes.data(), input.data(), bytes.size());
    std::string result(bytes.begin(), bytes.end());
    return result;

  } else if (precision == 4) {
    std::vector<uint8_t> bytes(sizeof(float) * input.size());
    for (size_t i = 0; i < input.size(); ++i) {
      float floatValue = static_cast<float>(input[i]);
      std::memcpy(bytes.data() + i * sizeof(float), &floatValue, sizeof(float));
    }
    std::string result(bytes.begin(), bytes.end());
    return result;

  } else {
    throw std::runtime_error("Precision must be 4 (32-bit) or 8 (64-bit)!");
  }
};

std::string sc::utils::encode_big_endian(const std::vector<double>& input, const int& precision) {

  if (precision == 8) {
    std::vector<uint8_t> bytes(sizeof(double) * input.size());
    
    for (size_t i = 0; i < input.size(); ++i) {
      double doubleValue = input[i];
      uint64_t value;
      std::memcpy(&value, &doubleValue, sizeof(double));
      for (size_t j = 0; j < sizeof(double); ++j) {
        bytes[i * sizeof(double) + j] = (value >> (8 * (sizeof(double) - 1 - j))) & 0xFF;
      }
    }

    std::string result(bytes.begin(), bytes.end());

    return result;

  } else if (precision == 4) {
    std::vector<uint8_t> bytes(sizeof(float) * input.size());
    
    for (size_t i = 0; i < input.size(); ++i) {
      float floatValue = static_cast<float>(input[i]);
      uint32_t value;
      std::memcpy(&value, &floatValue, sizeof(float));
      for (size_t j = 0; j < sizeof(float); ++j) {
        bytes[i * sizeof(float) + j] = (value >> (8 * (sizeof(float) - 1 - j))) & 0xFF;
      }
    }

    std::string result(bytes.begin(), bytes.end());

    return result;

  } else {
    throw std::runtime_error("Precision must be 4 (32-bit) or 8 (64-bit)!");
  }
}

std::vector<double> sc::utils::decode_little_endian(const std::string& str, const int& precision) {

  std::vector<unsigned char> bytes(str.begin(), str.end());

  if (precision != sizeof(double) && precision != sizeof(float)) {
    throw std::invalid_argument("Precision must be sizeof(double) or sizeof(float)!");
  }

  size_t bytes_size = bytes.size() / precision;
  std::vector<double> result(bytes_size);

  for (size_t i = 0; i < bytes_size; ++i) {
    if (precision == sizeof(double)) {
      double doubleValue;
      std::memcpy(&doubleValue, &bytes[i * precision], sizeof(double));
      result[i] = doubleValue;
    } else if (precision == sizeof(float)) {
      float floatValue;
      std::memcpy(&floatValue, &bytes[i * precision], sizeof(float));
      result[i] = static_cast<double>(floatValue);
    } else {
      throw std::runtime_error("Precision must be 4 (32-bit) or 8 (64-bit)!");
    }
  }

  return result;
};

std::vector<double> sc::utils::decode_big_endian(const std::string& str, const int& precision) {

  std::vector<unsigned char> bytes(str.begin(), str.end());

  if (precision != sizeof(double) && precision != sizeof(float)) {
    throw std::invalid_argument("Precision must be sizeof(double) or sizeof(float)!");
  }

  size_t bytes_size = bytes.size() / precision;
  std::vector<double> result(bytes_size);

  for (size_t i = 0; i < bytes_size; ++i) {

    if (precision == sizeof(double)) {
      uint64_t value = 0;
      for (int j = 0; j < precision; ++j) {
        value = (value << 8) | bytes[i * precision + j];
      }
      double doubleValue;
      std::memcpy(&doubleValue, &value, sizeof(double));
      result[i] = doubleValue;
    } else if (precision == sizeof(float)) {
      uint32_t value = 0;
      for (int j = 0; j < precision; ++j) {
        value = (value << 8) | bytes[i * precision + j];
      }
      float floatValue;
      std::memcpy(&floatValue, &value, sizeof(float));
      result[i] = static_cast<double>(floatValue);
    } else {
      throw std::runtime_error("Precision must be 4 (32-bit) or 8 (64-bit)!");
    }
  }

  return result;
};

std::string sc::utils::compress_zlib(const std::string& str) {

  std::vector<char> compressed_data;

  z_stream zs;

  memset(&zs, 0, sizeof(zs));

  if (deflateInit(&zs, Z_DEFAULT_COMPRESSION) != Z_OK) {
    throw std::runtime_error("deflateInit failed while initializing zlib for compression");
  }

  zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(str.data()));
  zs.avail_in = str.size();

  int ret;
  char outbuffer[32768];

  do {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
    zs.avail_out = sizeof(outbuffer);

    ret = deflate(&zs, Z_FINISH);

    if (compressed_data.size() < zs.total_out) {
      compressed_data.insert(compressed_data.end(), outbuffer, outbuffer + (zs.total_out - compressed_data.size()));
    }

  } while (ret == Z_OK);

  deflateEnd(&zs);

  return std::string(compressed_data.begin(), compressed_data.end());
};

std::string sc::utils::decompress_zlib(const std::string& compressed_string) {

  z_stream zs;

  memset(&zs, 0, sizeof(zs));

  inflateInit(&zs);

  zs.next_in = (Bytef*)compressed_string.data();

  zs.avail_in = compressed_string.size();

  int ret;

  char outbuffer[32768];

  std::string outstring;

  do {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);

    zs.avail_out = sizeof(outbuffer);

    ret = inflate(&zs, 0);

    if (outstring.size() < zs.total_out) {
      outstring.append(outbuffer, zs.total_out - outstring.size());
    }
  }

  while (ret == Z_OK);

  inflateEnd(&zs);

  return outstring;
};

std::string sc::utils::encode_base64(const std::string& str) {

  static const char* base64_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

  std::string encoded_data;

  encoded_data.reserve(((str.size() + 2) / 3) * 4);

  for (size_t i = 0; i < str.size(); i += 3) {

    int b = (str[i] & 0xFC) >> 2;
    encoded_data.push_back(base64_chars[b]);

    if (i + 1 < str.size()) {
      b = ((str[i] & 0x03) << 4) | ((str[i + 1] & 0xF0) >> 4);
      encoded_data.push_back(base64_chars[b]);

      if (i + 2 < str.size()) {
        b = ((str[i + 1] & 0x0F) << 2) | ((str[i + 2] & 0xC0) >> 6);
        encoded_data.push_back(base64_chars[b]);
        b = str[i + 2] & 0x3F;
        encoded_data.push_back(base64_chars[b]);

      } else {
        b = (str[i + 1] & 0x0F) << 2;
        encoded_data.push_back(base64_chars[b]);
        encoded_data.push_back('=');
      }

    } else {
      b = (str[i] & 0x03) << 4;
      encoded_data.push_back(base64_chars[b]);
      encoded_data.push_back('=');
      encoded_data.push_back('=');
    }
  }

  return encoded_data;
};

std::string sc::utils::decode_base64(const std::string& encoded_string) {

  std::string decoded_string;

  decoded_string.reserve((encoded_string.size() * 3) / 4);

  int val = 0;
  int valb = -8;
  for (char c : encoded_string) {
    if (c == '=') {
      valb -= 6;
      continue;
    }
    if (c >= 'A' && c <= 'Z') {
      c -= 'A';
    } else if (c >= 'a' && c <= 'z') {
      c -= 'a' - 26;
    } else if (c >= '0' && c <= '9') {
      c -= '0' - 52;
    } else if (c == '+') {
      c = 62;
    } else if (c == '/') {
      c = 63;
    } else {
      continue;
    }
    val = (val << 6) + c;
    valb += 6;
    if (valb >= 0) {
      decoded_string.push_back(char((val >> valb) & 0xFF));
      valb -= 8;
    }
  }

  return decoded_string;
};
