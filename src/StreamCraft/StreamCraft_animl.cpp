#include "StreamCraft_animl.hpp"
#include <string>
#include <regex>

int animl::extract_set_size(const pugi::xml_node& node, const std::string& name) {

  int set_size = 0;

  pugi::xml_node set = node.child(name.c_str());

  if (set) {

    std::vector<pugi::xml_node> set_nodes;

    for (pugi::xml_node child = set.first_child(); child; child = child.next_sibling()) {
      set_nodes.push_back(child);
    }

    set_size = set_nodes.size();
  }

  return set_size;
};



std::vector<animl::SAMPLE> animl::ANIML::extract_samples_by_indices(const std::vector<int>& idxs) {

  std::vector<animl::SAMPLE> samples;

  pugi::xml_node sample_set = root.child("SampleSet");

  if (sample_set) {

    std::vector<pugi::xml_node> sample_nodes;

    for (pugi::xml_node child = sample_set.first_child(); child; child = child.next_sibling()) {
      sample_nodes.push_back(child);
    }

    int number_samples = idxs.size();

    if (number_samples == 0) {
      std::cerr << "Warning: No samples to return!" << std::endl;
      return samples;
    }

    samples.resize(number_samples);

    for (int i = 0; i < number_samples; i++) {

      const int& index = idxs[i];
      const pugi::xml_node& node = sample_nodes[index];

      if (node) {
        animl::SAMPLE sample;
        sample.name = node.attribute("name").as_string();
        sample.sampleID = node.attribute("sampleID").as_string();
        samples[i] = sample;
      }
    }
  }

  // DUMMY FOR NOW

  return samples;
};



std::vector<animl::SAMPLE> animl::ANIML::extract_samples_by_attr(const std::string& attr, const std::vector<std::string>& names) {

  std::vector<animl::SAMPLE> samples;

  pugi::xml_node sample_set = root.child("SampleSet");

  if (sample_set) {

    int number_samples = names.size();

    if (number_samples == 0) {
      std::cerr << "Warning: No samples to return!" << std::endl;
      return samples;
    }

    samples.resize(number_samples);

    for (int i = 0; i < number_samples; i++) {

      const std::string& nm = names[i];
      const pugi::xml_node node = sample_set.find_child_by_attribute(attr.c_str(), nm.c_str());

      if (node) {
        animl::SAMPLE sample;
        sample.name = node.attribute("name").as_string();
        sample.sampleID = node.attribute("sampleID").as_string();
        samples[i] = sample;
      }
    }
  }

  // DUMMY FOR NOW

  return samples;
};



std::tuple<std::string, int> animl::extract_encoding_parameters(const std::string& str) {

  std::tuple<std::string, int> out;

  std::regex pattern(R"((\D+)(\d+))");

  std::smatch match;

  if (std::regex_match(str, match, pattern)) {
    out = {match[1].str(), std::stoi(match[2])};

  } else {
     throw("Enconding type could not by extracted!");
  }

  return out;
};



void animl::SERIES::extract(const pugi::xml_node& node) {
  name = node.attribute("name").as_string();
  dependency = node.attribute("dependency").as_string();
  seriesID = node.attribute("seriesID").as_string();
  seriesType = node.attribute("seriesType").as_string();
  plotScale = node.attribute("plotScale").as_string();

  pugi::xml_node value_set_node = node.first_child();
  ValueSetName = value_set_node.name();

  if (ValueSetName == "EncodedValueSet") {
    std::string encoded_string = value_set_node.child_value();
    std::string decoded_string = sc::decode_base64(encoded_string);
    std::tuple<std::string, int> param = extract_encoding_parameters(seriesType);
    // always outputs double vector
    EncodedValueSet = sc::decode_little_endian(decoded_string, std::get<1>(param) / 8);
  }
};



void animl::RESULT::extract(const pugi::xml_node& node) {
  name = node.attribute("name").as_string();

  std::string search_series_nodes = "./SeriesSet/Series";
  pugi::xpath_node_set ser_nodes = node.select_nodes(search_series_nodes.c_str());
  int n_ser_nodes = ser_nodes.size();

  if (n_ser_nodes > 0) {
    SeriesSet.resize(n_ser_nodes);
    for (int i = 0; i < n_ser_nodes; ++i) {
      pugi::xpath_node n = ser_nodes[i];
      SeriesSet[i].extract(n.node());
    }
  }

  std::string search_exp_nodes = "./ExperimentStepSet/ExperimentStep";
  pugi::xpath_node_set exp_nodes = node.select_nodes(search_exp_nodes.c_str());
  int n_exp_nodes = exp_nodes.size();

  if (n_exp_nodes > 0) {
    ExpStepSet.resize(n_exp_nodes);
    for (int i = 0; i < n_exp_nodes; ++i) {
      ExpStepSet[i].extract(exp_nodes[i].node());
    }
  }
};



void animl::EXPSTEP::extract(const pugi::xml_node& node) {
  name = node.attribute("name").as_string();
  expID = node.attribute("experimentStepID").as_string();

  pugi::xml_node technique_node = node.child("Technique");
  if (technique_node) {
    technique_name = technique_node.attribute("name").as_string();
    technique_uri = technique_node.attribute("uri").as_string();
  }

  std::string search_result_nodes = "./Result";
  pugi::xpath_node_set res_nodes = node.select_nodes(search_result_nodes.c_str());
  int n_res_nodes = res_nodes.size();

  if (n_res_nodes > 0) {
    ResultSet.resize(n_res_nodes);
    for (int i = 0; i < n_res_nodes; ++i) {
      ResultSet[i].extract(res_nodes[i].node());
    }
  }
};



std::vector<animl::EXPSTEP> animl::extract_experiment_step_by_indices(const pugi::xml_node& node, const std::vector<int>& idxs) {

  std::vector<EXPSTEP> exps;

  pugi::xml_node exp_set = node.child("ExperimentStepSet");

  if (exp_set) {

    std::vector<pugi::xml_node> exp_nodes;

    for (pugi::xml_node child = exp_set.first_child(); child; child = child.next_sibling()) {
      exp_nodes.push_back(child);
    }

    int number_exps = idxs.size();

    if (number_exps == 0) {
      std::cerr << "Warning: No experiment steps to return!" << std::endl;
      return exps;
    }

    exps.resize(number_exps);

    for (int i = 0; i < number_exps; i++) {

      const int& index = idxs[i];
      const pugi::xml_node& node = exp_nodes[index];

      if (node) {
        exps[i].extract(node);
      }
    }
  }

  // DUMMY FOR NOW

  return exps;
};



std::vector<animl::EXPSTEP> animl::extract_experiment_step_by_attr(const pugi::xml_node& node, const std::string& attr, const std::vector<std::string>& names) {

  std::vector<EXPSTEP> exps;

  pugi::xml_node exp_set = node.child("ExperimentStepSet");


  if (exp_set) {

    int number_exps = names.size();

    if (number_exps == 0) {
      std::cerr << "Warning: No samples to return!" << std::endl;
      return exps;
    }

    exps.resize(number_exps);

    for (int i = 0; i < number_exps; i++) {

      const std::string& nm = names[i];
      const pugi::xml_node node = exp_set.find_child_by_attribute(attr.c_str(), nm.c_str());

      if (node) {
        exps[i].extract(node);
      }
    }
  }

  // DUMMY FOR NOW

  return exps;
};



std::vector<animl::EXPSTEP> animl::extract_experiment_step_by_child_attr(const pugi::xml_node& node, const std::string& child_name, const std::string& attr, const std::vector<std::string>& names) {

  std::vector<EXPSTEP> exps;

  pugi::xml_node exp_set = node.child("ExperimentStepSet");

  if (exp_set) {

    int number_exps = names.size();

    if (number_exps == 0) {
      std::cerr << "Warning: No samples to return!" << std::endl;
      return exps;
    }

    exps.resize(number_exps);

    for (int i = 0; i < number_exps; i++) {

      const std::string& nm = names[i];

      std::string search_technique_node = "./ExperimentStep/"+ child_name + "[@name='" + nm +"']/..";

      pugi::xml_node node = exp_set.select_node(search_technique_node.c_str()).node();

      if (node) {
        exps[i].extract(node);
      }
    }
  }

  // DUMMY FOR NOW

  return exps;
};
