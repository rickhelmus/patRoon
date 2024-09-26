#ifndef STREAMCRAFT_ANIML_HPP
#define STREAMCRAFT_ANIML_HPP

#include <iostream>
#include <vector>
#include <numeric>

// #include <string>
// #include <regex>
// #include <tuple> 
// #include <cstring>
// #include <algorithm>
// #include <omp.h>

#define PUGIXML_HEADER_ONLY

#ifndef PUGIXML_PATH
#define PUGIXML_PATH "../../pugixml-1.14/src/pugixml.hpp"
#endif

#include PUGIXML_PATH

#include "StreamCraft_utils.hpp"

namespace animl {

    // Forward declaration of EXPSTEP
    class EXPSTEP;
    class RESULT;
    class SERIES;

    std::tuple<std::string, int> extract_encoding_parameters(const std::string& str);

    int extract_set_size(const pugi::xml_node& node, const std::string& name);

    std::vector<EXPSTEP> extract_experiment_step_by_indices(const pugi::xml_node& node, const std::vector<int>& idxs);

    std::vector<EXPSTEP> extract_experiment_step_by_attr(const pugi::xml_node& node, const std::string& attr, const std::vector<std::string>& names);

    std::vector<EXPSTEP> extract_experiment_step_by_child_attr(const pugi::xml_node& node, const std::string& child_name, const std::string& attr, const std::vector<std::string>& names);

    class TAG {
      public:
        std::string name;
        std::string value;
    };



    class UNIT {
      public:
        std::string label;
        std::string quantity;
        double factor;
        double exponent;
    };



    class SAMPLE {
      public:
        std::string name;
        std::string sampleID;
        std::vector<TAG> TagSet;
    };



    class PARAMETER {
      public:
        std::string name;
        std::string id;
        std::string type;
        unsigned value;
        UNIT value_unit;
    };



    class SERIES {
      public:
        std::string name;
        std::string dependency;
        std::string seriesID;
        std::string seriesType;
        std::string plotScale;
        std::string ValueSetName;
        std::vector<double> EncodedValueSet;
        std::vector<double> AutoIncrementedValueSet;
        std::vector<double> IndividualValueSet;
        UNIT valueUnit;
        void extract(const pugi::xml_node& node);
    };



    class CATEGORY {
      public:
        std::vector<PARAMETER> ParameterSet;
        std::vector<SERIES> SeriesSet;
        std::vector<CATEGORY> CategorySet;
    };



    class RESULT {
      public:
        std::string name;
        std::vector<SERIES> SeriesSet;
        int SeriesSetLength;
        std::string SeriesSetName;
        std::vector<EXPSTEP> ExpStepSet;
        void extract(const pugi::xml_node& node);
    };



    class EXPSTEP {
      public:
        std::string name;
        std::string expID;
        std::string technique_name;
        std::string technique_uri;
        std::vector<TAG> TagSet;
        std::vector<RESULT> ResultSet;
        void extract(const pugi::xml_node& node);
    };



    /// @brief AnIMl file representation.
    class ANIML {

      private:

        std::string file_path;

        pugi::xml_document doc;

        pugi::xml_parse_result loading_result;

        pugi::xml_node root;

        std::string format;

        std::string name;

        std::vector<SAMPLE> extract_samples_by_indices(const std::vector<int>& idxs);

        std::vector<SAMPLE> extract_samples_by_attr(const std::string& attr, const std::vector<std::string>& names);

      public:

        ANIML(const std::string& file) {

          file_path = file;

          const char* path = file.c_str();

          loading_result = doc.load_file(path, pugi::parse_default | pugi::parse_declaration | pugi::parse_pi);

          if (loading_result) {
            root = doc.document_element();

            if (!root) {
              std::cerr << "Warning: Root element is empty!" << std::endl;

            } else {
              format = root.name();

              if ("AnIML" != format) {
                std::cerr << "Warning: Root element must be AnIML!" << std::endl;
              }
            }

          } else {
            std::cerr << "Warning: AnIML file could not be opened!" << std::endl << loading_result.description() << std::endl;
          }

          name = root.name();

        };

        const std::string& get_name() {
          return name;
        };

        const pugi::xml_node& get_root() {
          return root;
        };

        int get_number_samples() {
          return extract_set_size(root, "SampleSet");
        };

        int get_number_experiment_steps() {
          return extract_set_size(root, "ExperimentStepSet");
        };

        std::vector<SAMPLE> get_samples(const std::vector<int>& indices = {},
                                              const std::vector<std::string>& names = {},
                                              const std::vector<std::string>& sampleID = {}) {

          std::vector<SAMPLE> samples;

          int number_samples = get_number_samples();

          if (number_samples == 0) {
            std::cerr << "Warning: No samples to return!" << std::endl;
            return samples;
          }

          if (indices.size() != 0) {
            samples = extract_samples_by_indices(indices);

          } else if (names.size() != 0) {
            samples = extract_samples_by_attr("name", names);

          } else if (sampleID.size() != 0) {
            samples = extract_samples_by_attr("sampleID", sampleID);

          } else {
            std::vector<int> idxs(number_samples);
            std::iota(idxs.begin(), idxs.end(), 0);
            samples = extract_samples_by_indices(idxs);
          }

          return samples;
        };

        std::vector<EXPSTEP> get_experiment_steps(const std::vector<int>& indices = {},
                                                  const std::vector<std::string>& names = {},
                                                  const std::vector<std::string>& techniques = {},
                                                  const std::vector<std::string>& experimentStepIDs = {}) {

          std::vector<EXPSTEP> exps;

          int number_exps = get_number_experiment_steps();

          if (number_exps == 0) {
            std::cerr << "Warning: No samples to return!" << std::endl;
            return exps;
          }

          if (indices.size() != 0) {
            exps = extract_experiment_step_by_indices(root, indices);

          } else if (names.size() != 0) {
            exps = extract_experiment_step_by_attr(root, "name", names);

          } else if (techniques.size() != 0) {
            exps = extract_experiment_step_by_child_attr(root, "Technique", "name", techniques);


          } else if (experimentStepIDs.size() != 0) {
            exps = extract_experiment_step_by_attr(root, "experimentStepID", experimentStepIDs);

          } else {
            std::vector<int> idxs(number_exps);
            std::iota(idxs.begin(), idxs.end(), 0);
            exps = extract_experiment_step_by_indices(root, idxs);
          }

          return exps;
        };
    };

  }; // animl

#endif // ANIML_HPP

#if defined(STREAMCRAFT_HEADER_ONLY) && !defined(STREAMCRAFT_ANIML_SOURCE)
#	define STREAMCRAFT_ANIML_SOURCE "StreamCraft_animl.cpp"
#	include STREAMCRAFT_ANIML_SOURCE
#endif
