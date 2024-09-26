#ifndef STREAMCRAFT_MZML_HPP
#define STREAMCRAFT_MZML_HPP

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <memory>

#define PUGIXML_HEADER_ONLY

#ifndef PUGIXML_PATH
#define PUGIXML_PATH "../../pugixml-1.14/src/pugixml.hpp"
#endif

#include PUGIXML_PATH

#include "StreamCraft_utils.hpp"

namespace sc {

  inline namespace mzml {

    const std::vector<std::string> mzml_possible_accessions_binary_data = {
      "MS:1000514", "MS:1000515", "MS:1000516", "MS:1000517",
      "MS:1000595", "MS:1000617", "MS:1000786", "MS:1000820",
      "MS:1000821", "MS:1000822", "MS:1002478", "MS:1002529",
      "MS:1002530", "MS:1002742", "MS:1002743", "MS:1002744",
      "MS:1002745", "MS:1002893", "MS:1003143", "MS:1003157",
      "MS:1003158", "MS:1003006"
    };

    const std::vector<std::string> mzml_possible_short_name_binary_data = {
      "mz", "intensity", "charge", "sn",
      "time", "wavelength", "other", "flowrate",
      "pressure", "temperature", "mean_charge", "resolution",
      "baseline", "noise", "sampled_noise_mz", "sampled_noise_intensity",
      "sampled_noise_baseline", "ion_mobility", "mass", "quadrupole_position_lower_bound_mz",
      "quadrupole_position_upper_bound_mz", "ion_mobility"
    };

    class MZML_BINARY_METADATA {
      public:
        int index;
        std::string precision_name;
        std::string precision_accession;
        int precision_int;
        std::string precision_type;
        std::string compression;
        bool compressed;
        std::string data_name;
        std::string data_accession;
        std::string data_value;
        std::string data_unit;
        std::string data_name_short;
    };

    class MZML_SPECTRUM {
      public:
        MZML_SPECTRUM(const pugi::xml_node& node) : spec(node) {};
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
        std::string extract_spec_title() const;
        double extract_scan_rt() const;
        int extract_scan_configuration_number() const;
        double extract_scan_mobility() const;
        std::string extract_scan_filter_string() const;
        int extract_scan_config() const;
        double extract_scan_injection_ion_time() const;
        int extract_precursor_scan() const;
        double extract_window_mz() const;
        double extract_window_mzlow() const;
        double extract_window_mzhigh() const;
        double extract_ion_mz() const;
        double extract_ion_intensity() const;
        int extract_ion_charge() const;
        std::string extract_activation_type() const;
        double extract_activation_ce() const;
        bool has_precursor() const { return spec.child("precursorList").child("precursor"); }
        bool has_selected_ion() const { return spec.child("precursorList").child("precursor").child("selectedIonList").child("selectedIon"); }
        bool has_activation() const { return spec.child("precursorList").child("precursor").child("activation"); }
        std::vector<MZML_BINARY_METADATA> extract_binary_metadata() const;
        std::vector<std::vector<double>> extract_binary_data(const std::vector<MZML_BINARY_METADATA>&  mtd) const;
      
      private:
        const pugi::xml_node& spec;
    };

    class MZML_CHROMATOGRAM {
      public:
        MZML_CHROMATOGRAM(const pugi::xml_node& node) : chrom(node) {};
        int extract_index() const;
        std::string extract_id() const;
        int extract_array_length() const;
        int extract_polarity() const;
        double extract_precursor_mz() const;
        std::string extract_activation_type() const;
        double extract_activation_ce() const;
        double extract_product_mz() const;
        bool has_precursor() const { return chrom.child("precursor"); }
        bool has_activation() const { return chrom.child("precursor").child("activation"); }
        bool has_product() const { return chrom.child("product"); }
        std::vector<std::vector<double>> extract_binary_data() const;
      
      private:
        const pugi::xml_node& chrom;
    };

    class MZML : public MS_READER {
      private:
        std::vector<pugi::xml_node> link_vector_spectra_nodes() const;
        std::vector<pugi::xml_node> link_vector_chrom_nodes() const;

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
        std::vector<pugi::xml_node> chrom_nodes;

        MZML(const std::string& file);

        std::vector<std::string> get_spectra_binary_short_names();
        std::vector<MZML_BINARY_METADATA> get_spectra_binary_metadata();
        void write_spectra(const std::vector<std::vector<std::vector<double>>>& spectra, const std::vector<std::string>& names, MS_SPECTRA_MODE mode, bool compress, bool save, std::string save_suffix);

        std::string get_format() override { return format; };
        int get_number_spectra() override;
        int get_number_chromatograms() override;
        int get_number_spectra_binary_arrays() override;
        std::string get_time_stamp() override;
        std::string get_type() override;
        std::vector<int> get_spectra_index(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_scan_number(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_array_length(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_level(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_configuration(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_mode(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_polarity(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_lowmz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_highmz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_bpmz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_bpint(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_tic(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_rt(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_mobility(std::vector<int> indices = {}) override;
        std::vector<int> get_spectra_precursor_scan(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_precursor_mz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_precursor_window_mz(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_precursor_window_mzlow(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_precursor_window_mzhigh(std::vector<int> indices = {}) override;
        std::vector<double> get_spectra_collision_energy(std::vector<int> indices = {}) override;
        std::vector<int> get_polarity() override;
        std::vector<int> get_mode() override;
        std::vector<int> get_level() override;
        std::vector<int> get_configuration() override;
        double get_min_mz() override;
        double get_max_mz() override;
        double get_start_rt() override;
        double get_end_rt() override;
        bool has_ion_mobility() override;
        MS_SUMMARY get_summary() override;
        MS_SPECTRA_HEADERS get_spectra_headers(std::vector<int> indices = {}) override;
        MS_CHROMATOGRAMS_HEADERS get_chromatograms_headers(std::vector<int> indices = {}) override;
        std::vector<std::vector<std::vector<double>>> get_spectra(std::vector<int> indices = {}) override;
        std::vector<std::vector<std::vector<double>>> get_chromatograms(std::vector<int> indices = {}) override;
        std::vector<std::vector<std::string>> get_software() override;
        std::vector<std::vector<std::string>> get_hardware() override;
        MS_SPECTRUM get_spectrum(const int& idx) override;
    }; // class MZML
  }; // namespace mzml
}; // namespace sc

#endif // STREAMCRAFT_MZML_HPP

#if defined(STREAMCRAFT_HEADER_ONLY) && !defined(STREAMCRAFT_MZML_SOURCE)
#	define STREAMCRAFT_MZML_SOURCE "StreamCraft_mzml.cpp"
#	include STREAMCRAFT_MZML_SOURCE
#endif
