#ifndef STREAMCRAFT_MZXML_HPP
#define STREAMCRAFT_MZXML_HPP

#include <iostream>
#include <vector>
#include <string>
#include <numeric>

#define PUGIXML_HEADER_ONLY

#ifndef PUGIXML_PATH
#define PUGIXML_PATH "../../pugixml-1.14/src/pugixml.hpp"
#endif

#include PUGIXML_PATH

#define STREAMCRAFT_HEADER_ONLY
#include "StreamCraft_utils.hpp"

namespace sc {

  inline namespace mzxml {

    class MZXML_BINARY_METADATA {
      public:
        int precision;
        std::string compression;
        bool compressed;
        std::string byte_order;
    };

    class MZXML_SPECTRUM {
      public:
        MZXML_SPECTRUM(const pugi::xml_node& node) : spec(node) {};
        int extract_spec_index() const;
        std::string extract_spec_id() const;
        int extract_spec_scan() const;
        int extract_spec_array_length() const;
        int extract_spec_level() const;
        int extract_spec_mode() const;
        int extract_spec_polarity() const;
        double extract_spec_lowmz() const;
        double extract_spec_highmz() const;
        double extract_spec_bpmz() const;
        double extract_spec_bpint() const;
        double extract_spec_tic() const;
        double extract_scan_rt() const;
        double extract_ion_mz() const;
        double extract_activation_ce() const;
        bool has_precursor() const { return spec.child("precursorMz"); }
        MZXML_BINARY_METADATA extract_binary_metadata() const;
        std::vector<std::vector<double>> extract_binary_data(const MZXML_BINARY_METADATA&  mtd) const;
      
      private:
        const pugi::xml_node& spec;
    };

    class MZXML : public MS_READER {
      private:
        std::vector<pugi::xml_node> link_vector_spectra_nodes() const;

      public:
        std::string file_path;
        std::string file_dir;
        std::string file_name;
        std::string file_extension;
        pugi::xml_document doc;
        pugi::xml_parse_result loading_result;
        pugi::xml_node root;
        std::string format;
        std::string name;
        std::vector<pugi::xml_node> spectra_nodes;

        MZXML(const std::string& file);

        std::vector<std::string> get_spectra_binary_short_names();
        MZXML_BINARY_METADATA get_spectra_binary_metadata();

        std::string get_format() override { return format; };
        int get_number_spectra() override;
        int get_number_chromatograms() override { return 0; };
        int get_number_spectra_binary_arrays() override;
        std::string get_time_stamp() override { return ""; }
        std::string get_type() override;
        std::vector<int> get_spectra_index(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_scan_number(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_array_length(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_level(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_configuration(std::vector<int> indices = {}) override {
          std::vector<int> configuration;
          return configuration;
        };
        std::vector<int> get_spectra_mode(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_polarity(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_lowmz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_highmz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_bpmz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_bpint(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_tic(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_rt(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_mobility(std::vector<int> indices = {}) override {
          std::vector<double> drift;
          return drift;
        };
        std::vector<int> get_spectra_precursor_scan(std::vector<int> indices = {}) override {
          std::vector<int> precursor_scan;
          return precursor_scan;
        };
        std::vector<double> get_spectra_precursor_mz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_precursor_window_mz(std::vector<int> indices = {}) override {
          std::vector<double> precursor_window_mz;
          return precursor_window_mz;
        };
        std::vector<double> get_spectra_precursor_window_mzlow(std::vector<int> indices = {}) override {
          std::vector<double> precursor_window_mzlow;
          return precursor_window_mzlow;
        };
        std::vector<double> get_spectra_precursor_window_mzhigh(std::vector<int> indices = {}) override {
          std::vector<double> precursor_window_mzhigh;
          return precursor_window_mzhigh;
        };
        std::vector<double> get_spectra_collision_energy(std::vector<int> indices = {}) override;
        std::vector<int> get_polarity() override;
        std::vector<int> get_mode() override;
        std::vector<int> get_level() override;
        std::vector<int> get_configuration() override {
          std::vector<int> configuration;
          return configuration;
        };
        double get_min_mz() override;
        double get_max_mz() override;
        double get_start_rt() override;
        double get_end_rt() override;
        bool has_ion_mobility() override { return false; };
        MS_SUMMARY get_summary() override;
        MS_SPECTRA_HEADERS get_spectra_headers(std::vector<int> indices = {}) override;
        MS_CHROMATOGRAMS_HEADERS get_chromatograms_headers(std::vector<int> indices = {}) override {
          MS_CHROMATOGRAMS_HEADERS chromatograms_headers;
          return chromatograms_headers;
        };
        std::vector<std::vector<std::vector<double>>> get_spectra(std::vector<int> indices = {}) override;
        std::vector<std::vector<std::vector<double>>> get_chromatograms(std::vector<int> indices = {}) override {
          std::vector<std::vector<std::vector<double>>> chromatograms;
          return chromatograms;
        };
        MS_SPECTRUM get_spectrum(const int& idx) override;
        std::vector<std::vector<std::string>> get_software() override;
        std::vector<std::vector<std::string>> get_hardware() override;
    }; // class MZXML
  }; // namespace mzxml
}; // namespace sc

#endif // STREAMCRAFT_MZXML_HPP

#if defined(STREAMCRAFT_HEADER_ONLY) && !defined(STREAMCRAFT_MZXML_SOURCE)
#	define STREAMCRAFT_MZXML_SOURCE "StreamCraft_mzxml.cpp"
#	include STREAMCRAFT_MZXML_SOURCE
#endif
