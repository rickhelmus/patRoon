#ifndef STREAMCRAFT_LIB_H
#define STREAMCRAFT_LIB_H

#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <unordered_set>
#include <cmath>
#include <numeric>
#include <memory>

namespace sc
{

  enum MS_SPECTRA_MODE
  {
    UNDEFINED,
    PROFILE,
    CENTROID
  };

  enum MS_SPECTRA_POLARITY
  {
    NEUTRAL,
    POSITIVE,
    NEGATIVE
  };

  struct MS_SPECTRUM
  {
    int index;
    int scan;
    int array_length;
    int level;
    int mode;
    int polarity;
    float lowmz;
    float highmz;
    float bpmz;
    float bpint;
    float tic;
    int configuration;
    float rt;
    float mobility;
    float window_mz;
    float window_mzlow;
    float window_mzhigh;
    float precursor_mz;
    float precursor_intensity;
    int precursor_charge;
    float activation_ce;
    int binary_arrays_count;
    std::vector<std::string> binary_names;
    std::vector<std::vector<float>> binary_data;
  };

  struct MS_SPECTRA_HEADERS
  {
    std::vector<int> index;
    std::vector<int> scan;
    std::vector<int> array_length;
    std::vector<int> level;
    std::vector<int> mode;
    std::vector<int> polarity;
    std::vector<float> lowmz;
    std::vector<float> highmz;
    std::vector<float> bpmz;
    std::vector<float> bpint;
    std::vector<float> tic;
    std::vector<int> configuration;
    std::vector<float> rt;
    std::vector<float> mobility;
    std::vector<float> window_mz;
    std::vector<float> window_mzlow;
    std::vector<float> window_mzhigh;
    std::vector<float> precursor_mz;
    std::vector<float> precursor_intensity;
    std::vector<int> precursor_charge;
    std::vector<float> activation_ce;

    void resize_all(int n)
    {
      index.resize(n);
      scan.resize(n);
      array_length.resize(n);
      level.resize(n);
      mode.resize(n);
      polarity.resize(n);
      lowmz.resize(n);
      highmz.resize(n);
      bpmz.resize(n);
      bpint.resize(n);
      tic.resize(n);
      configuration.resize(n);
      rt.resize(n);
      mobility.resize(n);
      window_mz.resize(n);
      window_mzlow.resize(n);
      window_mzhigh.resize(n);
      precursor_mz.resize(n);
      precursor_intensity.resize(n);
      precursor_charge.resize(n);
      activation_ce.resize(n);
    }

    size_t size() const
    {
      return index.size();
    }
  };

  struct MS_SUMMARY
  {
    std::string file_name;
    std::string file_path;
    std::string file_dir;
    std::string file_extension;
    int number_spectra;
    int number_chromatograms;
    int number_spectra_binary_arrays;
    std::string format;
    std::string time_stamp;
    std::vector<int> polarity;
    std::vector<int> mode;
    std::vector<int> level;
    std::vector<int> configuration;
    std::string type;
    float min_mz;
    float max_mz;
    float start_rt;
    float end_rt;
    bool has_ion_mobility;
  };

  struct MS_CHROMATOGRAMS_HEADERS
  {
    std::vector<int> index;
    std::vector<std::string> id;
    std::vector<int> array_length;
    std::vector<int> polarity;
    std::vector<float> precursor_mz;
    std::vector<float> activation_ce;
    std::vector<float> product_mz;

    void resize_all(int n)
    {
      index.resize(n);
      id.resize(n);
      array_length.resize(n);
      polarity.resize(n);
      precursor_mz.resize(n);
      activation_ce.resize(n);
      product_mz.resize(n);
    }

    size_t size() const
    {
      return index.size();
    }
  };

  struct MS_TARGETS
  {
    std::vector<int> index;
    std::vector<std::string> id;
    std::vector<int> level;
    std::vector<int> polarity;
    std::vector<bool> precursor;
    std::vector<float> mzmin;
    std::vector<float> mzmax;
    std::vector<float> rtmin;
    std::vector<float> rtmax;
    std::vector<float> mobilitymin;
    std::vector<float> mobilitymax;

    void resize_all(int n)
    {
      index.resize(n);
      id.resize(n);
      level.resize(n);
      polarity.resize(n);
      precursor.resize(n);
      mzmin.resize(n);
      mzmax.resize(n);
      rtmin.resize(n);
      rtmax.resize(n);
      mobilitymin.resize(n);
      mobilitymax.resize(n);
    }

    MS_TARGETS operator[](int i)
    {
      MS_TARGETS target;
      target.index.push_back(index[i]);
      target.id.push_back(id[i]);
      target.level.push_back(level[i]);
      target.polarity.push_back(polarity[i]);
      target.precursor.push_back(precursor[i]);
      target.mzmin.push_back(mzmin[i]);
      target.mzmax.push_back(mzmax[i]);
      target.rtmin.push_back(rtmin[i]);
      target.rtmax.push_back(rtmax[i]);
      target.mobilitymin.push_back(mobilitymin[i]);
      target.mobilitymax.push_back(mobilitymax[i]);
      return target;
    }

    size_t size() const
    {
      return index.size();
    }
  };

  struct MS_TARGETS_SPECTRA
  {
    std::vector<std::string> id;
    std::vector<int> polarity;
    std::vector<int> level;
    std::vector<float> pre_mz;
    std::vector<float> pre_mzlow;
    std::vector<float> pre_mzhigh;
    std::vector<float> pre_ce;
    std::vector<float> rt;
    std::vector<float> mobility;
    std::vector<float> mz;
    std::vector<float> intensity;

    void resize_all(int n)
    {
      id.resize(n);
      polarity.resize(n);
      level.resize(n);
      pre_mz.resize(n);
      pre_mzlow.resize(n);
      pre_mzhigh.resize(n);
      pre_ce.resize(n);
      rt.resize(n);
      mobility.resize(n);
      mz.resize(n);
      intensity.resize(n);
    }

    size_t size() const
    {
      return id.size();
    }

    int number_ids() const
    {
      std::unordered_set<std::string> unique_ids(id.begin(), id.end());
      return unique_ids.size();
    }

    MS_TARGETS_SPECTRA operator[](const std::string &unique_id) const
    {
      MS_TARGETS_SPECTRA target;
      for (size_t i = 0; i < id.size(); ++i)
      {
        if (id[i] == unique_id)
        {
          target.id.push_back(id[i]);
          target.polarity.push_back(polarity[i]);
          target.level.push_back(level[i]);
          target.pre_mz.push_back(pre_mz[i]);
          target.pre_mzlow.push_back(pre_mzlow[i]);
          target.pre_mzhigh.push_back(pre_mzhigh[i]);
          target.pre_ce.push_back(pre_ce[i]);
          target.rt.push_back(rt[i]);
          target.mobility.push_back(mobility[i]);
          target.mz.push_back(mz[i]);
          target.intensity.push_back(intensity[i]);
        }
      }
      return target;
    }
  };

  class MS_READER
  {
  public:
    MS_READER(const std::string &file) : file_(file) {}
    virtual ~MS_READER() = default;
    virtual int get_number_spectra() = 0;
    virtual int get_number_chromatograms() = 0;
    virtual int get_number_spectra_binary_arrays() = 0;
    virtual std::string get_format() = 0;
    virtual std::string get_type() = 0;
    virtual std::string get_time_stamp() = 0;
    virtual std::vector<int> get_polarity() = 0;
    virtual std::vector<int> get_mode() = 0;
    virtual std::vector<int> get_level() = 0;
    virtual std::vector<int> get_configuration() = 0;
    virtual float get_min_mz() = 0;
    virtual float get_max_mz() = 0;
    virtual float get_start_rt() = 0;
    virtual float get_end_rt() = 0;
    virtual bool has_ion_mobility() = 0;
    virtual MS_SUMMARY get_summary() = 0;
    virtual std::vector<int> get_spectra_index(std::vector<int> indices = {}) = 0;
    virtual std::vector<int> get_spectra_scan_number(std::vector<int> indices = {}) = 0;
    virtual std::vector<int> get_spectra_array_length(std::vector<int> indices = {}) = 0;
    virtual std::vector<int> get_spectra_level(std::vector<int> indices = {}) = 0;
    virtual std::vector<int> get_spectra_configuration(std::vector<int> indices = {}) = 0;
    virtual std::vector<int> get_spectra_mode(std::vector<int> indices = {}) = 0;
    virtual std::vector<int> get_spectra_polarity(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_lowmz(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_highmz(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_bpmz(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_bpint(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_tic(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_rt(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_mobility(std::vector<int> indices = {}) = 0;
    virtual std::vector<int> get_spectra_precursor_scan(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_precursor_mz(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_precursor_window_mz(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_precursor_window_mzlow(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_precursor_window_mzhigh(std::vector<int> indices = {}) = 0;
    virtual std::vector<float> get_spectra_collision_energy(std::vector<int> indices = {}) = 0;
    virtual MS_SPECTRA_HEADERS get_spectra_headers(std::vector<int> indices = {}) = 0;
    virtual MS_CHROMATOGRAMS_HEADERS get_chromatograms_headers(std::vector<int> indices = {}) = 0;
    virtual std::vector<std::vector<std::vector<float>>> get_spectra(std::vector<int> indices = {}) = 0;
    virtual std::vector<std::vector<std::vector<float>>> get_chromatograms(std::vector<int> indices = {}) = 0;
    virtual std::vector<std::vector<std::string>> get_software() = 0;
    virtual std::vector<std::vector<std::string>> get_hardware() = 0;
    virtual MS_SPECTRUM get_spectrum(const int &idx) = 0;

  protected:
    std::string file_;
  };

  // MARK: FUNCTIONS

  std::string encode_little_endian_from_float(const std::vector<float> &input, const int &precision);

  std::string encode_little_endian_from_double(const std::vector<double> &input, const int &precision);

  std::vector<float> decode_little_endian_to_float(const std::string &str, const int &precision);

  std::vector<double> decode_little_endian_to_double(const std::string &str, const int &precision);

  std::string encode_big_endian_from_float(const std::vector<float> &input, const int &precision);

  std::string encode_big_endian_from_double(const std::vector<double> &input, const int &precision);

  std::vector<float> decode_big_endian_to_float(const std::string &str, const int &precision);

  std::vector<double> decode_big_endian_to_double(const std::string &str, const int &precision);

  std::string compress_zlib(const std::string &str);

  std::string decompress_zlib(const std::string &compressed_string);

  std::string encode_base64(const std::string &str);

  std::string decode_base64(const std::string &encoded_string);

  std::unique_ptr<MS_READER> create_ms_reader(const std::string &file);

  // MARK: MZML
  inline namespace mzml
  {

    const std::vector<std::string> mzml_possible_accessions_binary_data = {
        "MS:1000514", "MS:1000515", "MS:1000516", "MS:1000517",
        "MS:1000595", "MS:1000617", "MS:1000786", "MS:1000820",
        "MS:1000821", "MS:1000822", "MS:1002478", "MS:1002529",
        "MS:1002530", "MS:1002742", "MS:1002743", "MS:1002744",
        "MS:1002745", "MS:1002893", "MS:1003143", "MS:1003157",
        "MS:1003158", "MS:1003006"};

    const std::vector<std::string> mzml_possible_short_name_binary_data = {
        "mz", "intensity", "charge", "sn",
        "time", "wavelength", "other", "flowrate",
        "pressure", "temperature", "mean_charge", "resolution",
        "baseline", "noise", "sampled_noise_mz", "sampled_noise_intensity",
        "sampled_noise_baseline", "ion_mobility", "mass", "quadrupole_position_lower_bound_mz",
        "quadrupole_position_upper_bound_mz", "ion_mobility"};

    class MZML_BINARY_METADATA
    {
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

    class MZML_SPECTRUM
    {
    public:
      MZML_SPECTRUM(const pugi::xml_node &node) : spec(node) {};
      int extract_spec_index() const;
      std::string extract_spec_id() const;
      int extract_spec_scan() const;
      int extract_spec_array_length() const;
      int extract_spec_level() const;
      int extract_spec_mode() const;
      int extract_spec_polarity() const;
      float extract_spec_lowmz() const;
      float extract_spec_highmz() const;
      float extract_spec_bpmz() const;
      float extract_spec_bpint() const;
      float extract_spec_tic() const;
      std::string extract_spec_title() const;
      float extract_scan_rt() const;
      int extract_scan_configuration_number() const;
      float extract_scan_mobility() const;
      std::string extract_scan_filter_string() const;
      int extract_scan_config() const;
      float extract_scan_injection_ion_time() const;
      int extract_precursor_scan() const;
      float extract_window_mz() const;
      float extract_window_mzlow() const;
      float extract_window_mzhigh() const;
      float extract_ion_mz() const;
      float extract_ion_intensity() const;
      int extract_ion_charge() const;
      std::string extract_activation_type() const;
      float extract_activation_ce() const;
      bool has_precursor() const { return spec.child("precursorList").child("precursor"); }
      bool has_selected_ion() const { return spec.child("precursorList").child("precursor").child("selectedIonList").child("selectedIon"); }
      bool has_activation() const { return spec.child("precursorList").child("precursor").child("activation"); }
      std::vector<MZML_BINARY_METADATA> extract_binary_metadata() const;
      std::vector<std::vector<float>> extract_binary_data(const std::vector<MZML_BINARY_METADATA> &mtd) const;

    private:
      const pugi::xml_node &spec;
    };

    class MZML_CHROMATOGRAM
    {
    public:
      MZML_CHROMATOGRAM(const pugi::xml_node &node) : chrom(node) {};
      int extract_index() const;
      std::string extract_id() const;
      int extract_array_length() const;
      int extract_polarity() const;
      float extract_precursor_mz() const;
      std::string extract_activation_type() const;
      float extract_activation_ce() const;
      float extract_product_mz() const;
      bool has_precursor() const { return chrom.child("precursor"); }
      bool has_activation() const { return chrom.child("precursor").child("activation"); }
      bool has_product() const { return chrom.child("product"); }
      std::vector<std::vector<float>> extract_binary_data() const;

    private:
      const pugi::xml_node &chrom;
    };

    class MZML : public sc::MS_READER
    {
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

      MZML(const std::string &file);

      std::vector<std::string> get_spectra_binary_short_names();
      std::vector<MZML_BINARY_METADATA> get_spectra_binary_metadata();
      void write_spectra(const std::vector<std::vector<std::vector<float>>> &spectra, const std::vector<std::string> &names, MS_SPECTRA_MODE mode, bool compress, bool save, std::string save_suffix);

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
      std::vector<float> get_spectra_lowmz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_highmz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_bpmz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_bpint(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_tic(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_rt(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_mobility(std::vector<int> indices = {}) override;
      std::vector<int> get_spectra_precursor_scan(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_precursor_mz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_precursor_window_mz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_precursor_window_mzlow(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_precursor_window_mzhigh(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_collision_energy(std::vector<int> indices = {}) override;
      std::vector<int> get_polarity() override;
      std::vector<int> get_mode() override;
      std::vector<int> get_level() override;
      std::vector<int> get_configuration() override;
      float get_min_mz() override;
      float get_max_mz() override;
      float get_start_rt() override;
      float get_end_rt() override;
      bool has_ion_mobility() override;
      MS_SUMMARY get_summary() override;
      MS_SPECTRA_HEADERS get_spectra_headers(std::vector<int> indices = {}) override;
      MS_CHROMATOGRAMS_HEADERS get_chromatograms_headers(std::vector<int> indices = {}) override;
      std::vector<std::vector<std::vector<float>>> get_spectra(std::vector<int> indices = {}) override;
      std::vector<std::vector<std::vector<float>>> get_chromatograms(std::vector<int> indices = {}) override;
      std::vector<std::vector<std::string>> get_software() override;
      std::vector<std::vector<std::string>> get_hardware() override;
      MS_SPECTRUM get_spectrum(const int &idx) override;
    }; // class MZML
  }; // namespace mzml

  // MARK: MZXML
  inline namespace mzxml
  {

    class MZXML_BINARY_METADATA
    {
    public:
      int precision;
      std::string compression;
      bool compressed;
      std::string byte_order;
    };

    class MZXML_SPECTRUM
    {
    public:
      MZXML_SPECTRUM(const pugi::xml_node &node) : spec(node) {};
      int extract_spec_index() const;
      std::string extract_spec_id() const;
      int extract_spec_scan() const;
      int extract_spec_array_length() const;
      int extract_spec_level() const;
      int extract_spec_mode() const;
      int extract_spec_polarity() const;
      float extract_spec_lowmz() const;
      float extract_spec_highmz() const;
      float extract_spec_bpmz() const;
      float extract_spec_bpint() const;
      float extract_spec_tic() const;
      float extract_scan_rt() const;
      float extract_ion_mz() const;
      float extract_activation_ce() const;
      bool has_precursor() const { return spec.child("precursorMz"); }
      MZXML_BINARY_METADATA extract_binary_metadata() const;
      std::vector<std::vector<float>> extract_binary_data(const MZXML_BINARY_METADATA &mtd) const;

    private:
      const pugi::xml_node &spec;
    };

    class MZXML : public MS_READER
    {
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

      MZXML(const std::string &file);

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
      std::vector<int> get_spectra_configuration(std::vector<int> indices = {}) override
      {
        std::vector<int> configuration;
        return configuration;
      };
      std::vector<int> get_spectra_mode(std::vector<int> indices = {}) override;
      std::vector<int> get_spectra_polarity(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_lowmz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_highmz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_bpmz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_bpint(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_tic(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_rt(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_mobility(std::vector<int> indices = {}) override
      {
        std::vector<float> drift;
        return drift;
      };
      std::vector<int> get_spectra_precursor_scan(std::vector<int> indices = {}) override
      {
        std::vector<int> precursor_scan;
        return precursor_scan;
      };
      std::vector<float> get_spectra_precursor_mz(std::vector<int> indices = {}) override;
      std::vector<float> get_spectra_precursor_window_mz(std::vector<int> indices = {}) override
      {
        std::vector<float> precursor_window_mz;
        return precursor_window_mz;
      };
      std::vector<float> get_spectra_precursor_window_mzlow(std::vector<int> indices = {}) override
      {
        std::vector<float> precursor_window_mzlow;
        return precursor_window_mzlow;
      };
      std::vector<float> get_spectra_precursor_window_mzhigh(std::vector<int> indices = {}) override
      {
        std::vector<float> precursor_window_mzhigh;
        return precursor_window_mzhigh;
      };
      std::vector<float> get_spectra_collision_energy(std::vector<int> indices = {}) override;
      std::vector<int> get_polarity() override;
      std::vector<int> get_mode() override;
      std::vector<int> get_level() override;
      std::vector<int> get_configuration() override
      {
        std::vector<int> configuration;
        return configuration;
      };
      float get_min_mz() override;
      float get_max_mz() override;
      float get_start_rt() override;
      float get_end_rt() override;
      bool has_ion_mobility() override { return false; };
      MS_SUMMARY get_summary() override;
      MS_SPECTRA_HEADERS get_spectra_headers(std::vector<int> indices = {}) override;
      MS_CHROMATOGRAMS_HEADERS get_chromatograms_headers(std::vector<int> indices = {}) override
      {
        MS_CHROMATOGRAMS_HEADERS chromatograms_headers;
        return chromatograms_headers;
      };
      std::vector<std::vector<std::vector<float>>> get_spectra(std::vector<int> indices = {}) override;
      std::vector<std::vector<std::vector<float>>> get_chromatograms(std::vector<int> indices = {}) override
      {
        std::vector<std::vector<std::vector<float>>> chromatograms;
        return chromatograms;
      };
      MS_SPECTRUM get_spectrum(const int &idx) override;
      std::vector<std::vector<std::string>> get_software() override;
      std::vector<std::vector<std::string>> get_hardware() override;
    }; // class MZXML
  }; // namespace mzxml

  // MARK: MS_FILE
  class MS_FILE
  {
  private:
    const std::vector<std::string> possible_formats = {"mzML", "mzXML"};

  public:
    std::string file_path;
    std::string file_dir;
    std::string file_name;
    std::string file_extension;
    std::string format;
    int format_case;
    std::unique_ptr<MS_READER> ms;

    MS_FILE(const std::string &file);

    int get_number_spectra() { return ms->get_number_spectra(); }
    int get_number_chromatograms() { return ms->get_number_chromatograms(); }
    int get_number_spectra_binary_arrays() { return ms->get_number_spectra_binary_arrays(); }
    std::string get_format() { return ms->get_format(); }
    std::string get_type() { return ms->get_type(); }
    std::string get_time_stamp() { return ms->get_time_stamp(); }
    std::vector<int> get_polarity() { return ms->get_polarity(); }
    std::vector<int> get_mode() { return ms->get_mode(); }
    std::vector<int> get_level() { return ms->get_level(); }
    std::vector<int> get_configuration() { return ms->get_configuration(); }
    float get_min_mz() { return ms->get_min_mz(); }
    float get_max_mz() { return ms->get_max_mz(); }
    float get_start_rt() { return ms->get_start_rt(); }
    float get_end_rt() { return ms->get_end_rt(); }
    bool has_ion_mobility() { return ms->has_ion_mobility(); }
    MS_SUMMARY get_summary() { return ms->get_summary(); }
    std::vector<int> get_spectra_index(std::vector<int> indices = {}) { return ms->get_spectra_index(indices); }
    std::vector<int> get_spectra_scan_number(std::vector<int> indices = {}) { return ms->get_spectra_scan_number(indices); }
    std::vector<int> get_spectra_array_length(std::vector<int> indices = {}) { return ms->get_spectra_array_length(indices); }
    std::vector<int> get_spectra_level(std::vector<int> indices = {}) { return ms->get_spectra_level(indices); }
    std::vector<int> get_spectra_mode(std::vector<int> indices = {}) { return ms->get_spectra_mode(indices); }
    std::vector<int> get_spectra_polarity(std::vector<int> indices = {}) { return ms->get_spectra_polarity(indices); }
    std::vector<float> get_spectra_lowmz(std::vector<int> indices = {}) { return ms->get_spectra_lowmz(indices); }
    std::vector<float> get_spectra_highmz(std::vector<int> indices = {}) { return ms->get_spectra_highmz(indices); }
    std::vector<float> get_spectra_bpmz(std::vector<int> indices = {}) { return ms->get_spectra_bpmz(indices); }
    std::vector<float> get_spectra_bpint(std::vector<int> indices = {}) { return ms->get_spectra_bpint(indices); }
    std::vector<float> get_spectra_tic(std::vector<int> indices = {}) { return ms->get_spectra_tic(indices); }
    std::vector<int> get_spectra_configuration(std::vector<int> indices = {}) { return ms->get_spectra_configuration(indices); }
    std::vector<float> get_spectra_rt(std::vector<int> indices = {}) { return ms->get_spectra_rt(indices); }
    std::vector<float> get_spectra_drift(std::vector<int> indices = {}) { return ms->get_spectra_mobility(indices); }
    std::vector<int> get_spectra_precursor_scan(std::vector<int> indices = {}) { return ms->get_spectra_precursor_scan(indices); }
    std::vector<float> get_spectra_precursor_mz(std::vector<int> indices = {}) { return ms->get_spectra_precursor_mz(indices); }
    std::vector<float> get_spectra_precursor_window_mz(std::vector<int> indices = {}) { return ms->get_spectra_precursor_window_mz(indices); }
    std::vector<float> get_spectra_precursor_window_mzlow(std::vector<int> indices = {}) { return ms->get_spectra_precursor_window_mzlow(indices); }
    std::vector<float> get_spectra_precursor_window_mzhigh(std::vector<int> indices = {}) { return ms->get_spectra_precursor_window_mzhigh(indices); }
    std::vector<float> get_spectra_collision_energy(std::vector<int> indices = {}) { return ms->get_spectra_collision_energy(indices); }
    MS_SPECTRA_HEADERS get_spectra_headers(std::vector<int> indices = {}) { return ms->get_spectra_headers(indices); }
    MS_CHROMATOGRAMS_HEADERS get_chromatograms_headers(std::vector<int> indices = {}) { return ms->get_chromatograms_headers(indices); }
    std::vector<std::vector<std::vector<float>>> get_spectra(std::vector<int> indices = {}) { return ms->get_spectra(indices); }
    std::vector<std::vector<std::vector<float>>> get_chromatograms(std::vector<int> indices = {}) { return ms->get_chromatograms(indices); }
    std::vector<std::vector<std::string>> get_software() { return ms->get_software(); }
    std::vector<std::vector<std::string>> get_hardware() { return ms->get_hardware(); }
    MS_SPECTRUM get_spectrum(const int &index) { return ms->get_spectrum(index); }
    MS_TARGETS_SPECTRA get_spectra_targets(const MS_TARGETS &targets, const sc::MS_SPECTRA_HEADERS &hd, const float &minIntLv1, const float &minIntLv2);
  };
}; // namespace sc
#endif // STREAMCRAFT_LIB_H
