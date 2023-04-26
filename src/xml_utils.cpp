#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <zlib.h>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include "external_libraries.hpp"
#include "xml_utils.h"
#include <Rcpp.h>




std::list<std::vector<std::string>> xml_utils::mzml_instrument_parser(const pugi::xml_node& node_mzml) {

  std::list<std::vector<std::string>> output;
  std::vector<std::string> names;
  std::vector<std::string> vals;

  std::string search_ref = "//referenceableParamGroup";
  pugi::xpath_node xp_ref = node_mzml.select_node(search_ref.c_str());

  if (xp_ref.node() != NULL) {
    for (pugi::xml_node temp: xp_ref.node().children()) {
      if (temp.attribute("name") != NULL) {
        names.push_back(temp.attribute("name").as_string());
        std::string val = temp.attribute("value").as_string();
        if (val != "") {
          vals.push_back(temp.attribute("value").as_string());
        } else {
          vals.push_back(temp.attribute("name").as_string());
        }
      }
    }
  }

  std::string search_inst = "//instrumentConfiguration";
  pugi::xpath_node xp_inst = node_mzml.select_node(search_inst.c_str());

  if (xp_inst.node() != NULL) {
    for (pugi::xml_node temp: xp_inst.node().children()) {
      if (temp.attribute("name") != NULL) {
        names.push_back(temp.attribute("name").as_string());
        std::string val = temp.attribute("value").as_string();
        if (val != "") {
          vals.push_back(temp.attribute("value").as_string());
        } else {
          vals.push_back(temp.attribute("name").as_string());
        }
      }
    }
  }

  std::string search_config = "//componentList/child::node()";
  pugi::xpath_node_set xps_config = node_mzml.select_nodes(search_config.c_str());

  if (xps_config.size() > 0) {
    for (pugi::xpath_node_set::const_iterator it = xps_config.begin(); it != xps_config.end(); ++it) {
      pugi::xpath_node node = *it;
      for (pugi::xml_node temp: node.node().children())
      {
        names.push_back(node.node().name());
        vals.push_back(temp.attribute("name").as_string());
      }
    }
  }

  if (names.size() > 0) {
    output.push_back(names);
    output.push_back(vals);
  }

  return output;
}




std::list<std::vector<std::string>> xml_utils::mzxml_instrument_parser(const pugi::xml_node& node_mzxml) {

  std::list<std::vector<std::string>> output;
  std::vector<std::string> names;
  std::vector<std::string> vals;

  std::string search_inst = "//msInstrument/child::node()[starts-with(name(), 'ms')]";
  pugi::xpath_node_set xps_inst = node_mzxml.select_nodes(search_inst.c_str());

  if (xps_inst.size() > 0) {
    for (pugi::xpath_node_set::const_iterator it = xps_inst.begin(); it != xps_inst.end(); ++it) {
      pugi::xpath_node node = *it;
      names.push_back(node.node().attribute("category").as_string());
      vals.push_back(node.node().attribute("value").as_string());
    }
  }

  if (names.size() > 0) {
    output.push_back(names);
    output.push_back(vals);
  }

  return output;
}




Rcpp::List xml_utils::parse_instrument(const pugi::xml_node& root) {

  Rcpp::List instrument_list;

  std::list<std::vector<std::string>> output;

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {
    output = xml_utils::mzml_instrument_parser(root.child("mzML"));

  } else if (strcmp("mzXML", root_name) == 0) {
    output = xml_utils::mzxml_instrument_parser(root);

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
    return instrument_list;
  }

  auto it = output.begin();
  instrument_list["category"] = Rcpp::wrap(*it);
  std::advance(it, 1);
  instrument_list["value"] = Rcpp::wrap(*it);
  instrument_list.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");

  return instrument_list;
}




std::list<std::vector<std::string>> xml_utils::mzml_software_parser(const pugi::xml_node& node_mzml) {

  std::list<std::vector<std::string>> output;
  std::vector<std::string> names;
  std::vector<std::string> ids;
  std::vector<std::string> version;

  std::string search_software = "//softwareList/child::node()";
  pugi::xpath_node_set xps_software = node_mzml.select_nodes(search_software.c_str());

  if (xps_software.size() > 0) {
    for (pugi::xpath_node_set::const_iterator it = xps_software.begin(); it != xps_software.end(); ++it) {
      pugi::xpath_node node = *it;
      for (pugi::xml_node temp: node.node().children()) {
        std::string name = temp.attribute("name").as_string();
        if (name != "") {
          names.push_back(name);
          ids.push_back(node.node().attribute("id").as_string());
          version.push_back(node.node().attribute("version").as_string());
        }
      }
    }
  }

  output.push_back(names);
  output.push_back(ids);
  output.push_back(version);

  return output;
}




std::list<std::vector<std::string>> xml_utils::mzxml_software_parser(const pugi::xml_node& node_mzxml) {

  std::list<std::vector<std::string>> output;
  std::vector<std::string> names;
  std::vector<std::string> ids;
  std::vector<std::string> version;

  std::string search_software = "//msInstrument/child::node()[starts-with(name(), 'soft')]";
  pugi::xpath_node_set xps_software = node_mzxml.select_nodes(search_software.c_str());

  if (xps_software.size() > 0) {
    for (pugi::xpath_node_set::const_iterator it = xps_software.begin(); it != xps_software.end(); ++it) {
      pugi::xpath_node node = *it;
      names.push_back(node.node().attribute("name").as_string());
      ids.push_back(node.node().attribute("type").as_string());
      version.push_back(node.node().attribute("version").as_string());
    }
  }

  output.push_back(names);
  output.push_back(ids);
  output.push_back(version);

  return output;
}




Rcpp::List xml_utils::parse_software(const pugi::xml_node& root) {

  Rcpp::List software_list;

  std::list<std::vector<std::string>> output;

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {
    output = xml_utils::mzml_software_parser(root.child("mzML"));

  } else if (strcmp("mzXML", root_name) == 0) {
    output = xml_utils::mzxml_software_parser(root);

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
    return software_list;
  }

  auto it2 = output.begin();
  software_list["name"] = Rcpp::wrap(*it2);
  std::advance(it2, 1);
  software_list["id"] = Rcpp::wrap(*it2);
  std::advance(it2, 1);
  software_list["version"] = Rcpp::wrap(*it2);
  software_list.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");

  return software_list;
}




void xml_utils::mzml_spectra_headers_parser(const pugi::xml_node& node_mzml, xml_utils::spectraHeaders& output) {

  std::string search_run = "//run";

  pugi::xpath_node xps_run = node_mzml.select_node(search_run.c_str());

  pugi::xml_node spec_list = xps_run.node().child("spectrumList");

  if (spec_list != NULL) {

    std::vector<pugi::xml_node> spectra;

    for (pugi::xml_node child = spec_list.first_child(); child; child = child.next_sibling()) {
      spectra.push_back(child);
    }

    int number_spectra = spectra.size();

    output.index.resize(number_spectra);
    output.scan.resize(number_spectra);
    output.traces.resize(number_spectra);
    output.polarity.resize(number_spectra);
    output.bpcmz.resize(number_spectra);
    output.bpcint.resize(number_spectra);
    output.ticint.resize(number_spectra);
    output.level.resize(number_spectra);
    output.mode.resize(number_spectra);
    output.mzlow.resize(number_spectra);
    output.mzhigh.resize(number_spectra);
    output.rt.resize(number_spectra);
    output.drift.resize(number_spectra);
    output.pre_scan.resize(number_spectra);
    output.pre_mz.resize(number_spectra);
    output.pre_loweroffset.resize(number_spectra);
    output.pre_upperoffset.resize(number_spectra);
    output.pre_ce.resize(number_spectra);

    for (int i = 0; i < number_spectra; i++) {

      const pugi::xml_node spec = spectra[i];

      output.index[i] = spec.attribute("index").as_int();

      std::string scan = spec.attribute("id").as_string();
      std::size_t poslastEqual = scan.rfind('=');
      int scan_n = std::stoi(scan.substr(poslastEqual + 1));
      // std::sscanf(scan.c_str(), "%*[^=]=%d", &scan_n);
      output.scan[i] = scan_n;

      output.traces[i] = spec.attribute("defaultArrayLength").as_int();
      if (std::isnan(output.traces[i])) output.traces[i] = NA_INTEGER;

      pugi::xml_node node_scan = spec.child("scanList").child("scan");

      pugi::xml_node rt_node = node_scan.find_child_by_attribute("cvParam", "name", "scan start time");
      std::string rt_unit = rt_node.attribute("unitName").as_string();
      double rt_val = rt_node.attribute("value").as_double();
      if (rt_unit == "minute") rt_val = rt_val * 60;
      if (std::isnan(rt_val)) rt_val = NA_INTEGER;
      output.rt[i] = rt_val;

      pugi::xml_node drift_node = node_scan.find_child_by_attribute("cvParam", "name", "ion mobility drift time");
      if (drift_node != NULL) {
        // std::string drift_unit = node_drift.attribute("unitName").as_string();
        double drift_val = drift_node.attribute("value").as_double();
        // if (rt_unit == "millisecond") drift_t = drift_t / 0.001;
        output.drift[i] = drift_val;
      } else {
        output.drift[i] = NA_REAL;
      }

      pugi::xml_node level_node = spec.find_child_by_attribute("cvParam", "name", "ms level");
      output.level[i] = level_node.attribute("value").as_int();
      if (std::isnan(output.level[i])) output.level[i] = NA_INTEGER;

      pugi::xml_node pol_pos = spec.find_child_by_attribute("cvParam", "accession", "MS:1000130");
      pugi::xml_node pol_neg = spec.find_child_by_attribute("cvParam", "accession", "MS:1000129");

      if (pol_pos != NULL) {
        output.polarity[i] = "positive";
      } else if (pol_neg != NULL) {
        output.polarity[i] = "negative";
      } else {
        output.polarity[i] = "NA";
      }

      pugi::xml_node centroid = spec.find_child_by_attribute("cvParam", "accession", "MS:1000127");
      pugi::xml_node profile = spec.find_child_by_attribute("cvParam", "accession", "MS:1000128");

      if (centroid != NULL) {
        output.mode[i] = "centroid";
      } else if (profile != NULL) {
        output.mode[i] = "profile";
      } else {
        output.mode[i] = "NA";
      }

      pugi::xml_node mzlow_node = spec.find_child_by_attribute("cvParam", "name", "lowest observed m/z");
      output.mzlow[i] = mzlow_node.attribute("value").as_double();
      if (std::isnan(output.mzlow[i])) output.mzlow[i] = NA_REAL;

      pugi::xml_node mzhigh_node = spec.find_child_by_attribute("cvParam", "name", "highest observed m/z");
      output.mzhigh[i] = mzhigh_node.attribute("value").as_double();
      if (std::isnan(output.mzhigh[i])) output.mzhigh[i] = NA_REAL;

      pugi::xml_node bpcmz_node = spec.find_child_by_attribute("cvParam", "name", "base peak m/z");
      output.bpcmz[i] = bpcmz_node.attribute("value").as_double();
      if (std::isnan(output.bpcmz[i])) output.bpcmz[i] = NA_REAL;

      pugi::xml_node bpcint_node = spec.find_child_by_attribute("cvParam", "name", "base peak intensity");
      output.bpcint[i] = bpcint_node.attribute("value").as_double();
      if (std::isnan(output.bpcint[i])) output.bpcint[i] = NA_REAL;

      pugi::xml_node ticint_node = spec.find_child_by_attribute("cvParam", "name", "total ion current");
      output.ticint[i] = ticint_node.attribute("value").as_double();
      if (std::isnan(output.ticint[i])) output.ticint[i] = NA_REAL;

      pugi::xml_node precursor = spec.child("precursorList").child("precursor");

      if (precursor != NULL) {

        std::string pre_scan_str = precursor.attribute("spectrumRef").as_string();
        if (pre_scan_str == "") {
          output.pre_scan[i] = NA_INTEGER;
        } else {
          std::size_t pre_poslastEqual = pre_scan_str.rfind('=');
          int pre_scan_n = std::stoi(pre_scan_str.substr(pre_poslastEqual + 1));
          output.pre_scan[i] = pre_scan_n;
        }

        pugi::xml_node isolation = precursor.child("isolationWindow");

        pugi::xml_node pre_mz_node = isolation.find_child_by_attribute("cvParam", "name", "isolation window target m/z");
        output.pre_mz[i] = pre_mz_node.attribute("value").as_double();
        if (std::isnan(output.pre_mz[i])) output.pre_mz[i] = NA_REAL;

        pugi::xml_node pre_loweroffset_node = isolation.find_child_by_attribute("cvParam", "name", "isolation window lower offset");
        output.pre_loweroffset[i] = pre_loweroffset_node.attribute("value").as_double();
        if (std::isnan(output.pre_loweroffset[i])) output.pre_loweroffset[i] = NA_REAL;


        pugi::xml_node pre_upperoffset_node = isolation.find_child_by_attribute("cvParam", "name", "isolation window upper offset");
        output.pre_upperoffset[i] = pre_upperoffset_node.attribute("value").as_double();
        if (std::isnan(output.pre_upperoffset[i])) output.pre_upperoffset[i] = NA_REAL;

        pugi::xml_node activation = precursor.child("activation");

        pugi::xml_node pre_ce_node = activation.find_child_by_attribute("cvParam", "name", "collision energy");
        output.pre_ce[i] = pre_ce_node.attribute("value").as_double();
        if (std::isnan(output.pre_ce[i])) output.pre_ce[i] = NA_REAL;

      } else {
        output.pre_scan[i] = NA_INTEGER;
        output.pre_mz[i] = NA_REAL;
        output.pre_loweroffset[i] = NA_REAL;
        output.pre_upperoffset[i] = NA_REAL;
        output.pre_ce[i] = NA_REAL;
      }
    }
  }
}




void xml_utils::mzxml_spectra_headers_parser(const pugi::xml_node& node_mzxml, xml_utils::spectraHeaders& output) {

  std::string search_run = "//msRun";

  pugi::xpath_node xps_run = node_mzxml.select_node(search_run.c_str());

  std::vector<pugi::xml_node> spectra;

  for (pugi::xml_node child = xps_run.node().child("scan"); child; child = child.next_sibling()) {
    spectra.push_back(child);
  }

  int number_spectra = spectra.size();

  output.index.resize(number_spectra);
  output.scan.resize(number_spectra);
  output.traces.resize(number_spectra);
  output.polarity.resize(number_spectra);
  output.bpcmz.resize(number_spectra);
  output.bpcint.resize(number_spectra);
  output.ticint.resize(number_spectra);
  output.level.resize(number_spectra);
  output.mode.resize(number_spectra);
  output.mzlow.resize(number_spectra);
  output.mzhigh.resize(number_spectra);
  output.rt.resize(number_spectra);
  output.drift.resize(number_spectra);
  output.pre_scan.resize(number_spectra);
  output.pre_mz.resize(number_spectra);
  output.pre_loweroffset.resize(number_spectra);
  output.pre_upperoffset.resize(number_spectra);
  output.pre_ce.resize(number_spectra);

  for (int i = 0; i < number_spectra; i++) {

    const pugi::xml_node spec = spectra[i];

    output.index[i] = i;

    output.scan[i] = spec.attribute("num").as_int();

    output.traces[i] = spec.attribute("peaksCount").as_int();
    if (std::isnan(output.traces[i])) output.traces[i] = NA_REAL;

    std::string rt = spec.attribute("retentionTime").as_string();
    double rt_n;
    std::sscanf(rt.c_str(), "%*[^0123456789]%lf", &rt_n);
    char last_char = '\0';
    std::sscanf(rt.c_str() + rt.size() - 1, "%c", &last_char);
    if (last_char != 'S') rt_n = rt_n * 60;
    output.rt[i] = rt_n;

    output.drift[i] = NA_REAL;

    output.level[i] = spec.attribute("msLevel").as_int();
    if (std::isnan(output.level[i])) output.level[i] = NA_INTEGER;

    std::string pol_sign = spec.attribute("polarity").as_string();

    if (pol_sign == "+") {
      output.polarity[i] = "positive";
    } else if (pol_sign == "-") {
      output.polarity[i] = "negative";
    } else {
      output.polarity[i] = "NA";
    }

    int centroided = spec.attribute("centroided").as_int();

    if (centroided == 1) {
      output.mode[i] = "centroid";
    } else if (centroided == 0) {
      output.mode[i] = "profile";
    } else {
      output.mode[i] = "NA";
    }

    output.mzlow[i] = spec.attribute("lowMz").as_double();
    if (std::isnan(output.mzlow[i])) output.mzlow[i] = NA_REAL;

    output.mzhigh[i] = spec.attribute("highMz").as_double();
    if (std::isnan(output.mzhigh[i])) output.mzhigh[i] = NA_REAL;

    output.bpcmz[i] = spec.attribute("basePeakMz").as_double();
    if (std::isnan(output.bpcmz[i])) output.bpcmz[i] = NA_REAL;

    output.bpcint[i] = spec.attribute("basePeakIntensity").as_double();
    if (std::isnan(output.bpcint[i])) output.bpcint[i] = NA_REAL;

    output.ticint[i] = spec.attribute("totIonCurrent").as_double();
    if (std::isnan(output.ticint[i])) output.ticint[i] = NA_REAL;

    pugi::xml_node precursor = spec.child("precursorMz");

    if (precursor != NULL) {
      output.pre_mz[i] = precursor.text().as_double();
      if (std::isnan(output.pre_mz[i])) output.pre_mz[i] = NA_REAL;

      output.pre_ce[i] = spec.attribute("collisionEnergy").as_double();
      if (std::isnan(output.pre_ce[i])) output.pre_ce[i] = NA_REAL;

    } else {
      output.pre_mz[i] = NA_REAL;
      output.pre_ce[i] = NA_REAL;
    }

    output.pre_scan[i] = NA_INTEGER;
    output.pre_loweroffset[i] = NA_REAL;
    output.pre_upperoffset[i] = NA_REAL;
  }
}




xml_utils::spectraHeaders xml_utils::parse_spectra_headers(const pugi::xml_node& root) {

  xml_utils::spectraHeaders output;

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {
    xml_utils::mzml_spectra_headers_parser(root.child("mzML"), output);

  } else if (strcmp("mzXML", root_name) == 0) {
    xml_utils::mzxml_spectra_headers_parser(root, output);

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
  }

  return output;
}




Rcpp::List xml_utils::spectraHeaders_to_list(const xml_utils::spectraHeaders& headers_cpp) {

  Rcpp::StringVector polarity = Rcpp::wrap(headers_cpp.polarity);
  Rcpp::StringVector mode = Rcpp::wrap(headers_cpp.mode);

  for (int i = 0; i < polarity.size(); i++) {
    if (polarity(i) == "NA") {
      polarity(i) = NA_STRING;
    }
    if (mode(i) == "NA") {
      mode(i) = NA_STRING;
    }
  }

  Rcpp::List headers;

  headers["index"] = headers_cpp.index;
  headers["scan"] = headers_cpp.scan;
  headers["traces"] = headers_cpp.traces;
  headers["polarity"] = polarity;
  headers["mode"] = mode;
  headers["level"] = headers_cpp.level;
  headers["rt"] = headers_cpp.rt;
  headers["drift"] = headers_cpp.drift;
  headers["bpc_mz"] = headers_cpp.bpcmz;
  headers["bpc_intensity"] = headers_cpp.bpcint;
  headers["tic_intensity"] = headers_cpp.ticint;
  headers["mzlow"] = headers_cpp.mzlow;
  headers["mzhigh"] = headers_cpp.mzhigh;
  headers["pre_scan"] = headers_cpp.pre_scan;
  headers["pre_mz"] = headers_cpp.pre_mz;
  headers["pre_ce"] = headers_cpp.pre_ce;
  headers["pre_loweroffset"] = headers_cpp.pre_loweroffset;
  headers["pre_upperoffset"] = headers_cpp.pre_upperoffset;

  headers.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");

  return headers;
}




xml_utils::spectraHeaders xml_utils::list_to_spectraHeaders(const Rcpp::List& run) {

  xml_utils::spectraHeaders output;

  std::vector<int> index = run["index"];
  int number_spectra = index.size();

  Rcpp::IntegerVector pre_scan = run["pre_scan"];
  Rcpp::NumericVector pre_mz = run["pre_mz"];
  Rcpp::NumericVector pre_ce = run["pre_ce"];

  for (int i = 0; i < number_spectra; i++) {

      if (pre_scan[i] == NA_INTEGER) {
        pre_scan[i] = -1;
      }

      if (pre_mz[i] == NA_REAL) {
        pre_mz[i] = nan("");
      }

      if (pre_ce[i] == NA_REAL) {
      pre_ce[i] = nan("");
      }
  }

  output.index.resize(number_spectra);
  output.scan.resize(number_spectra);
  output.traces.resize(number_spectra);
  // output.polarity.resize(number_spectra);
  output.bpcmz.resize(number_spectra);
  output.bpcint.resize(number_spectra);
  output.ticint.resize(number_spectra);
  output.level.resize(number_spectra);
  // output.mode.resize(number_spectra);
  // output.mzlow.resize(number_spectra);
  // output.mzhigh.resize(number_spectra);
  output.rt.resize(number_spectra);
  // output.drift.resize(number_spectra);
  output.pre_scan.resize(number_spectra);
  output.pre_mz.resize(number_spectra);
  // output.pre_loweroffset.resize(number_spectra);
  // output.pre_upperoffset.resize(number_spectra);
  output.pre_ce.resize(number_spectra);

  output.index = index;
  output.scan = Rcpp::as<std::vector<int>>(run["scan"]);
  output.traces = Rcpp::as<std::vector<int>>(run["traces"]);
  output.level = Rcpp::as<std::vector<int>>(run["level"]);
  output.rt = Rcpp::as<std::vector<double>>(run["rt"]);
  // output.drift = Rcpp::as<std::vector<double>>(run["drift"]);
  output.bpcmz = Rcpp::as<std::vector<double>>(run["bpc_mz"]);
  output.bpcint = Rcpp::as<std::vector<double>>(run["bpc_intensity"]);
  output.ticint = Rcpp::as<std::vector<double>>(run["tic_intensity"]);
  output.pre_scan = Rcpp::as<std::vector<int>>(pre_scan);
  output.pre_mz = Rcpp::as<std::vector<double>>(pre_mz);
  output.pre_ce = Rcpp::as<std::vector<double>>(pre_ce);

  return output;
}




void xml_utils::mzml_chromatograms_headers_parser(
    const pugi::xml_node& node_mzml, xml_utils::chromatogramsHeaders& output) {

  std::string search_run = "//run";

  pugi::xpath_node xps_run = node_mzml.select_node(search_run.c_str());

  pugi::xml_node chrom_list = xps_run.node().child("chromatogramList");

  if (chrom_list != NULL) {

    std::vector<pugi::xml_node> chromatograms;

    for (pugi::xml_node child = chrom_list.first_child(); child; child = child.next_sibling()) {
      chromatograms.push_back(child);
    }

    int number_chromatograms = chromatograms.size();

    if (number_chromatograms > 0) {

      output.index.resize(number_chromatograms);
      output.id.resize(number_chromatograms);
      output.traces.resize(number_chromatograms);
      output.polarity.resize(number_chromatograms);
      output.pre_mz.resize(number_chromatograms);
      output.pre_ce.resize(number_chromatograms);
      output.pro_mz.resize(number_chromatograms);

      for (int i = 0; i < number_chromatograms; i++) {

        const pugi::xml_node chrom = chromatograms[i];

        output.index[i] = chrom.attribute("index").as_int();

        std::string id = chrom.attribute("id").as_string();
        output.id[i] = id;

        output.traces[i] = chrom.attribute("defaultArrayLength").as_int();
        if (std::isnan(output.traces[i])) output.traces[i] = NA_INTEGER;

        pugi::xml_node pol_pos = chrom.find_child_by_attribute("cvParam", "accession", "MS:1000130");
        pugi::xml_node pol_neg = chrom.find_child_by_attribute("cvParam", "accession", "MS:1000129");

        if (pol_pos != NULL) {
          output.polarity[i] = "positive";
        } else if (pol_neg != NULL) {
          output.polarity[i] = "negative";
        } else {
          output.polarity[i] = "NA";
        }

        pugi::xml_node precursor = chrom.child("precursor");

        if (precursor != NULL) {

          pugi::xml_node isolation = precursor.child("isolationWindow");

          pugi::xml_node pre_mz_node = isolation.find_child_by_attribute("cvParam", "name", "isolation window target m/z");
          output.pre_mz[i] = pre_mz_node.attribute("value").as_double();
          if (std::isnan(output.pre_mz[i])) output.pre_mz[i] = NA_REAL;

          pugi::xml_node activation = precursor.child("activation");

          pugi::xml_node pre_ce_node = activation.find_child_by_attribute("cvParam", "name", "collision energy");
          output.pre_ce[i] = pre_ce_node.attribute("value").as_double();
          if (std::isnan(output.pre_ce[i])) output.pre_ce[i] = NA_REAL;

        } else {
          output.pre_mz[i] = NA_REAL;
          output.pre_ce[i] = NA_REAL;
        }

        pugi::xml_node product = chrom.child("product");

        if (product != NULL) {

          pugi::xml_node isolation = product.child("isolationWindow");

          pugi::xml_node pro_mz_node = isolation.find_child_by_attribute("cvParam", "name", "isolation window target m/z");
          output.pro_mz[i] = pro_mz_node.attribute("value").as_double();
          if (std::isnan(output.pro_mz[i])) output.pro_mz[i] = NA_REAL;

        } else {
          output.pro_mz[i] = NA_REAL;
        }
      }
    }
  }
}




xml_utils::chromatogramsHeaders xml_utils::parse_chromatograms_headers(const pugi::xml_node& root) {

  xml_utils::chromatogramsHeaders output;

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {
    xml_utils::mzml_chromatograms_headers_parser(root.child("mzML"), output);

  } else if (strcmp("mzXML", root_name) == 0) {

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
  }

  return output;
}




Rcpp::List xml_utils::chromatogramsHeaders_to_list(const xml_utils::chromatogramsHeaders& headers_cpp) {

  Rcpp::StringVector polarity = Rcpp::wrap(headers_cpp.polarity);

  for (int i = 0; i < polarity.size(); i++) {
    if (polarity(i) == "NA") {
      polarity(i) = NA_STRING;
    }
  }

  Rcpp::List headers;

  headers["index"] = headers_cpp.index;
  headers["id"] = headers_cpp.id;
  headers["traces"] = headers_cpp.traces;
  headers["polarity"] = polarity;
  headers["pre_ce"] = headers_cpp.pre_ce;
  headers["pre_mz"] = headers_cpp.pre_mz;
  headers["pro_mz"] = headers_cpp.pro_mz;

  headers.attr("class") = Rcpp::CharacterVector::create("data.table", "data.frame");

  return headers;
}




xml_utils::runSummary xml_utils::run_summary(
    xml_utils::spectraHeaders& spec_headers,
    xml_utils::chromatogramsHeaders& chrom_headers) {

  xml_utils::runSummary output;

  Rcpp::CharacterVector na_charvec(1, NA_STRING);
  std::vector<int> empty_int_vec(1, NA_INTEGER);

  output.spectra_number = spec_headers.index.size();
  output.chromatograms_number = chrom_headers.index.size();

  if (output.spectra_number <= -1) output.spectra_number = 0;
  if (output.chromatograms_number <= -1) output.chromatograms_number = 0;

  if (output.spectra_number > 0) {

    std::set<std::string> mode(spec_headers.mode.begin(), spec_headers.mode.end());
    std::vector<std::string> mode_v(mode.begin(), mode.end());

    if (mode_v.size() == 1 && mode_v[0] == "NA") {
      output.mode = na_charvec;
    } else {
      output.mode = Rcpp::wrap(mode_v);
    }

    std::set<int> level(spec_headers.level.begin(), spec_headers.level.end());
    std::vector<int> level_v(level.begin(), level.end());
    output.levels = level_v;

    auto mzlow = std::min_element(spec_headers.mzlow.begin(), spec_headers.mzlow.end());
    output.mz_low = *mzlow;

    auto mzhigh = std::max_element(spec_headers.mzhigh.begin(), spec_headers.mzhigh.end());
    output.mz_high = *mzhigh;

    auto rtmin = std::min_element(spec_headers.rt.begin(), spec_headers.rt.end());
    output.rt_start = *rtmin;

    auto rtmax = std::max_element(spec_headers.rt.begin(), spec_headers.rt.end());
    output.rt_end = *rtmax;

    std::set<std::string> polarity(spec_headers.polarity.begin(), spec_headers.polarity.end());
    std::vector<std::string> polarity_v(polarity.begin(), polarity.end());

    if (polarity_v.size() == 1 && polarity_v[0] == "NA") {
      output.polarity = na_charvec;
    } else {
      output.polarity = Rcpp::wrap(polarity_v);
    }

    bool has_im = std::all_of(spec_headers.drift.begin(), spec_headers.drift.end(), [](double d) { return std::isnan(d); });
    output.has_ion_mobility = !has_im;

  } else {

    output.mode = na_charvec;
    output.levels = empty_int_vec;
    output.mz_low = NA_REAL;
    output.mz_high = NA_REAL;
    output.rt_start = NA_REAL;
    output.rt_end = NA_REAL;
    output.polarity = na_charvec;
    output.has_ion_mobility = false;
  }

  return output;
}




std::string xml_utils::encode_base64(const std::string& str) {

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
}




std::string xml_utils::decode_base64(const std::string& encoded_string) {

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

  // int padding = 0;
  // if (encoded_string.size() > 0 && encoded_string[encoded_string.size()-1] == '=') padding++;
  // if (encoded_string.size() > 1 && encoded_string[encoded_string.size()-2] == '=') padding++;
  //
  // const size_t in_len = encoded_string.size();
  // size_t out_len = (in_len * 6) / 8 - padding;
  // std::string out(out_len, '\0');
  //
  // size_t i = 0, j = 0;
  // int bits = 0, value = 0;
  // for (i = 0; i < in_len; i++) {
  //   char c = encoded_string[i];
  //   if (c == '=') {
  //     bits = 0;
  //     value = 0;
  //     continue;
  //   }
  //   value = (value << 6) | xml_utils::base64_chars.find(c);
  //   bits += 6;
  //   if (bits >= 8) {
  //     bits -= 8;
  //     out[j++] = (value >> bits) & 0xFF;
  //   }
  // }
  // return out;

  // size_t in_len = encoded_string.size();
  // size_t i = 0;
  // size_t j = 0;
  // size_t padding = 0;
  // unsigned char char_array_4[4], char_array_3[3];
  // std::string decoded_string;
  //
  // while (in_len-- && (encoded_string[i] != '=') && xml_utils::is_base64(encoded_string[i])) {
  //   char_array_4[j++] = encoded_string[i];
  //   i++;
  //   if (j == 4) {
  //     for (j = 0; j < 4; j++) {
  //       char_array_4[j] = xml_utils::base64_chars.find(char_array_4[j]);
  //     }
  //     char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
  //     char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
  //     char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];
  //
  //     for (j = 0; j < 3; j++) {
  //       decoded_string += char_array_3[j];
  //     }
  //     j = 0;
  //   }
  // }
  //
  // if (j) {
  //   for (size_t i = j; i < 4; i++) {
  //     char_array_4[i] = 0;
  //   }
  //   for (size_t i = 0; i < 4; i++) {
  //     char_array_4[i] = xml_utils::base64_chars.find(char_array_4[i]);
  //   }
  //   char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
  //   char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
  //   char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];
  //
  //   for (size_t i = 0; i < padding; i++) {
  //     decoded_string += char_array_3[i];
  //   }
  // }
  //
  // return decoded_string;
}




std::string xml_utils::compress_zlib(const std::string& str) {

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
}




std::string xml_utils::decompress_zlib(const std::string& compressed_string) {

  // std::string output_string;
  //
  // z_stream zs;
  // memset(&zs, 0, sizeof(zs));
  //
  // inflateInit(&zs);
  //
  // // if (inflateInit(&zs) != Z_OK) {
  // //   throw std::runtime_error("inflateInit failed!");
  // // }
  //
  // zs.next_in = (Bytef*)compressed_string.data();
  // zs.avail_in = compressed_string.size();
  //
  // int ret;
  // char outbuffer[32768];
  //
  // do {
  //   zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
  //   zs.avail_out = sizeof(outbuffer);
  //
  //   ret = inflate(&zs, 0);
  //
  //   if (output_string.size() < zs.total_out) {
  //     output_string.append(outbuffer, zs.total_out - output_string.size());
  //   }
  //
  // } while (ret == Z_OK);
  //
  // inflateEnd(&zs);
  //
  // // if (inflateEnd(&zs) != Z_OK) {
  // //   throw std::runtime_error("inflateEnd failed!");
  // // }
  //
  // return output_string;

  z_stream zs;
  memset(&zs, 0, sizeof(zs));
  inflateInit(&zs);
  zs.next_in = (Bytef*)compressed_string.data();
  zs.avail_in = compressed_string.size();

  int ret;
  char outbuffer[32768];
  std::string outstring;
  do
  {
    zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
    zs.avail_out = sizeof(outbuffer);
    ret = inflate(&zs, 0);
    if (outstring.size() < zs.total_out)
    {
      outstring.append(outbuffer, zs.total_out - outstring.size());
    }
  }
  while (ret == Z_OK);
  inflateEnd(&zs);

  return outstring;
}




std::string xml_utils::encode_little_endian(const Rcpp::NumericVector& input) {

  size_t n = input.size();

  std::string result(n * sizeof(double), '\0');

  const uint8_t* bytes = reinterpret_cast<const uint8_t*>(input.begin());

  uint8_t* output_bytes = reinterpret_cast<uint8_t*>(&result[0]);

  for (size_t i = 0; i < n; ++i) {
    uint64_t little_endian_value;
    std::memcpy(&little_endian_value, bytes + i * sizeof(double), sizeof(double));
    std::memcpy(output_bytes + i * sizeof(uint64_t), &little_endian_value, sizeof(uint64_t));
  }

  return result;
}




Rcpp::NumericVector xml_utils::decode_little_endian(const std::string& str, const int& precision) {

  std::vector<unsigned char> byteVec(str.begin(), str.end());

  int array_len = (byteVec.size() / precision);

  Rcpp::NumericVector result(array_len);

  for (int i = 0; i < array_len; ++i) {
    result[i] = reinterpret_cast<double&>(byteVec[i * precision]);
  }

  return result;
}




std::string xml_utils::encode_big_endian(const Rcpp::NumericVector& input) {

  size_t n = input.size();

  std::vector<uint8_t> result(n * sizeof(double));

  double* values = reinterpret_cast<double*>(result.data());

  for (size_t i = 0; i < n; ++i) {
    uint64_t big_endian_value;
    std::memcpy(&big_endian_value, &input[i], sizeof(double));
    std::reverse(reinterpret_cast<uint8_t*>(&big_endian_value), reinterpret_cast<uint8_t*>(&big_endian_value) + sizeof(uint64_t));
    std::memcpy(values + i, &big_endian_value, sizeof(uint64_t));
  }

  std::string str_result(reinterpret_cast<char*>(result.data()), result.size());

  return str_result;
}




Rcpp::NumericVector xml_utils::decode_big_endian(const std::string& str, const int& precision) {

  Rcpp::NumericVector result;

  const uint8_t* bytes = reinterpret_cast<const uint8_t*>(str.c_str());
  size_t byte_count = str.size();

  size_t precision_t = static_cast<size_t>(precision);

  for (size_t i = 0; i < byte_count; i += precision_t) {

    uint64_t value = 0;

    for (size_t j = 0; j < precision_t; ++j) {
      value = (value << precision) | bytes[i + j];
    }

    uint64_t little_endian_value;
    std::memcpy(&little_endian_value, &value, sizeof(little_endian_value));
    double double_value;
    std::memcpy(&double_value, &little_endian_value, sizeof(double_value));
    result.push_back(double_value);
  }

  return result;
}




std::vector<int> xml_utils::mzml_get_precision(pugi::xml_node& node) {

  pugi::xml_node node_binary_list = node.child("binaryDataArrayList");

  std::vector<int> precision;

  for (pugi::xml_node bin: node_binary_list.children("binaryDataArray")) {

    std::string precision_str;
    int precision_int;

    pugi::xml_node node_float = bin.find_child_by_attribute("cvParam", "accession", "MS:1000523");
    pugi::xml_node node_integer = bin.find_child_by_attribute("cvParam", "accession", "MS:1000522");

    if (node_float != NULL) {
      precision_str = node_float.attribute("name").as_string();
    } else {
      precision_str = node_integer.attribute("name").as_string();
    }

    std::sscanf(precision_str.c_str(), "%d%*c", &precision_int);

    precision_int = precision_int/8;

    precision.push_back(precision_int);
  }

  return precision;
}



int xml_utils::mzxml_get_precision(pugi::xml_node& node) {

  int precision = node.child("peaks").attribute("precision").as_int();

  precision = precision/8;

  return precision;
}




std::vector<std::string> xml_utils::mzml_get_compression(pugi::xml_node& node) {

  pugi::xml_node node_binary_list = node.child("binaryDataArrayList");

  std::vector<std::string> compression_out;

  for (pugi::xml_node bin: node_binary_list.children("binaryDataArray")) {

    std::string compression;

    pugi::xml_node node_comp = bin.find_child_by_attribute("cvParam", "accession", "MS:1000574");

    if (node_comp != NULL) {
      compression = node_comp.attribute("name").as_string();

      if (compression == "zlib" || compression == "zlib compression") {
        compression = "gzip";
      } else {
        compression = "none";
      }

    } else {
      compression = "none";
    }

    compression_out.push_back(compression);
  }

  return compression_out;
}




std::string xml_utils::mzxml_get_compression(pugi::xml_node& node) {

  std::string compression = node.child("peaks").attribute("compressionType").as_string();

  if (compression == "zlib" || compression == "zlib compression") {
    compression = "gzip";
  } else {
    compression = "none";
  }

  return compression;
}




Rcpp::CharacterVector xml_utils::mzml_get_binary_type(pugi::xml_node& node) {
  pugi::xml_node node_binary_list = node.child("binaryDataArrayList");

  Rcpp::CharacterVector type;

  for (pugi::xml_node bin: node_binary_list.children("binaryDataArray")) {
    pugi::xml_node node_float = bin.find_child_by_attribute("cvParam", "accession", "MS:1000523");
    pugi::xml_node node_integer = bin.find_child_by_attribute("cvParam", "accession", "MS:1000522");

    if (node_float != NULL) {
      type.push_back("double");

    } else if (node_integer != NULL) {
      type.push_back("integer");

    } else {
      type.push_back("double");
    }
  }

  return type;
}




Rcpp::CharacterVector xml_utils::mzml_get_binary_names(pugi::xml_node& node) {

  pugi::xml_node node_binary_list = node.child("binaryDataArrayList");

  Rcpp::CharacterVector bin_name_out;

  int counter = 0;

  for (pugi::xml_node bin: node_binary_list.children("binaryDataArray")) {

    pugi::xml_node node_mz = bin.find_child_by_attribute("cvParam", "accession", "MS:1000514");
    pugi::xml_node node_int = bin.find_child_by_attribute("cvParam", "accession", "MS:1000515");
    pugi::xml_node node_rt = bin.find_child_by_attribute("cvParam", "accession", "MS:1000595");
    pugi::xml_node node_other = bin.find_child_by_attribute("cvParam", "accession", "MS:1000786");

    if (node_mz != NULL) {
      bin_name_out.push_back("mz");

    } else if (node_int != NULL) {
      bin_name_out.push_back("intensity");

    } else if (node_rt != NULL) {
      bin_name_out.push_back("rt");

    } else if (node_other != NULL) {
      std::string name = node_other.attribute("value").as_string();
      name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
      bin_name_out.push_back(name);

    } else {
      std::string name = "val_" + std::to_string(counter);
      bin_name_out.push_back(name);
    }

    counter++;
  }

  return bin_name_out;
}




Rcpp::NumericMatrix xml_utils::mzml_parse_binary_data_from_spectrum_node(
  const pugi::xml_node& node,
  const std::vector<int>& precision,
  const std::vector<std::string>& compression,
  const Rcpp::CharacterVector& cols)
{

  int number_traces = node.attribute("defaultArrayLength").as_int();

  pugi::xml_node node_binary_list = node.child("binaryDataArrayList");

  int number_bins = node_binary_list.attribute("count").as_int();

  Rcpp::NumericMatrix output(number_traces, number_bins); //number_bins

  if (number_traces == 0) return  output;

  int counter = 0;

  for (auto it2 = node_binary_list.children("binaryDataArray").begin(); it2 != node_binary_list.children("binaryDataArray").end(); ++it2) {

    const pugi::xml_node& bin = *it2;

    if (counter > cols.size()) {
      break;
    }

    pugi::xml_node node_binary = bin.child("binary");
    const char* encoded_string = node_binary.child_value();

    std::string decoded_string = decode_base64(encoded_string);

    if (compression[counter] == "gzip") {
      decoded_string = xml_utils::decompress_zlib(decoded_string);
    }

    output(Rcpp::_, counter) = xml_utils::decode_little_endian(decoded_string, precision[counter]);

    if (cols[counter] == "rt") {
      pugi::xml_node node_unit = bin.find_child_by_attribute("cvParam", "unitCvRef", "UO");
      std::string unit = node_unit.attribute("unitName").as_string();

      if (unit == "minute") {
        output(Rcpp::_, counter) = output(Rcpp::_, counter) * 60;
      }
    }

    // std::cout << "the n rows: " << output.nrow() << " the n cols: " << output.cols() << " " << output(0, 0) <<  std::endl;

    counter++;
  }

  Rcpp::CharacterVector f_cols = Rcpp::head(cols, number_bins);

  Rcpp::colnames(output) = f_cols;

  return output;
}



Rcpp::NumericMatrix xml_utils::mzxml_parse_binary_data_from_spectrum_node(
    const pugi::xml_node& node,
    const int& precision,
    const std::string& compression,
    const Rcpp::CharacterVector& cols)
{

  int number_traces = node.attribute("peaksCount").as_int();

  if (number_traces == 0) return  Rcpp::NumericMatrix (0, 0);

  const char* encoded_string = node.child("peaks").child_value();

  std::string decoded_string = decode_base64(encoded_string);

  if (compression == "gzip") {
    decoded_string = xml_utils::decompress_zlib(decoded_string);
  }

  std::string endian_ness = node.child("peaks").attribute("byteOrder").as_string();

  Rcpp::NumericVector result;

  if (endian_ness == "network") {
    result = xml_utils::decode_big_endian(decoded_string, precision);

  } else {
    result = xml_utils::decode_little_endian(decoded_string, precision);
  }

  int n_row = result.size() / 2;

  Rcpp::NumericMatrix output(n_row, 2);

  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < 2; j++) {
      output(i, j) = result[i * 2 + j];
    }
  }

  Rcpp::colnames(output) = cols;

  return output;
}




Rcpp::List xml_utils::parse_spectra(const pugi::xml_node& root) {

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {

    pugi::xml_node node = root.child("mzML").child("run").child("spectrumList");

    pugi::xml_node first_spec = node.child("spectrum");
    const std::vector<int> precision = xml_utils::mzml_get_precision(first_spec);
    const std::vector<std::string> compression = xml_utils::mzml_get_compression(first_spec);
    const Rcpp::CharacterVector cols = xml_utils::mzml_get_binary_names(first_spec);

    std::vector<pugi::xml_node> spectra;

    for (pugi::xml_node child = node.first_child(); child; child = child.next_sibling()) {
      spectra.push_back(child);
    }

    int number_spectra = spectra.size();

    Rcpp::List list_output(number_spectra);

    if (number_spectra > 0) {

      for (int i = 0; i < number_spectra; i++) {
        list_output[i] = xml_utils::mzml_parse_binary_data_from_spectrum_node(spectra[i], precision, compression, cols);
      }

    }

    return list_output;

  } else if (strcmp("mzXML", root_name) == 0) {

    const Rcpp::CharacterVector cols = {"mz","intenisty"};

    pugi::xml_node node = root.child("msRun");

    pugi::xml_node first_spec = node.child("scan");
    const int precision = xml_utils::mzxml_get_precision(first_spec);
    const std::string compression = xml_utils::mzxml_get_compression(first_spec);

    std::vector<pugi::xml_node> spectra;

    for (pugi::xml_node child = node.child("scan"); child; child = child.next_sibling()) {
      spectra.push_back(child);
    }

    int number_spectra = spectra.size();

    Rcpp::List list_output(number_spectra);

    if (number_spectra > 0) {

      for (int i = 0; i < number_spectra; i++) {
        list_output[i] = xml_utils::mzxml_parse_binary_data_from_spectrum_node(spectra[i], precision, compression, cols);
      }

      return list_output;
    }

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
  }

  Rcpp::List list_output;
  return list_output;
}




Rcpp::List xml_utils::parse_partial_spectra(
    const pugi::xml_node& root, Rcpp::IntegerVector& index) {

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {

    pugi::xml_node node = root.child("mzML").child("run").child("spectrumList");

    pugi::xml_node first_spec = node.child("spectrum");
    const std::vector<int> precision = xml_utils::mzml_get_precision(first_spec);
    const std::vector<std::string> compression = xml_utils::mzml_get_compression(first_spec);
    const Rcpp::CharacterVector cols = xml_utils::mzml_get_binary_names(first_spec);

    std::vector<pugi::xml_node> spectra;

    for (pugi::xml_node child = node.first_child(); child; child = child.next_sibling()) {
      spectra.push_back(child);
    }

    int number_spectra = spectra.size();

    // remove any larger than the number_spectra
    Rcpp::LogicalVector overhead = index <= number_spectra;
    index = index[overhead];

    // remove duplicated indexes
    std::set<int> index_set(index.begin(), index.end());
    std::vector<int> index_unique(index_set.begin(), index_set.end());

    int number_index = index_unique.size();

    // updates index by reference with index_unique
    index = Rcpp::wrap(index_unique);

    std::sort(index.begin(), index.end());

    Rcpp::List list_output(number_index);

    if (number_index > 0) {

      int counter = 0;

      for (int i : index) {
        list_output[counter] = xml_utils::mzml_parse_binary_data_from_spectrum_node(spectra[i], precision, compression, cols);
        counter++;
      }

      return list_output;
    }

  } else if (strcmp("mzXML", root_name) == 0) {

    const Rcpp::CharacterVector cols = {"mz","intenisty"};

    pugi::xml_node node = root.child("msRun");

    pugi::xml_node first_spec = node.child("scan");
    const int precision = xml_utils::mzxml_get_precision(first_spec);
    const std::string compression = xml_utils::mzxml_get_compression(first_spec);

    std::vector<pugi::xml_node> spectra;

    for (pugi::xml_node child = node.child("scan"); child; child = child.next_sibling()) {
      spectra.push_back(child);
    }

    int number_spectra = spectra.size();

    // remove any larger than the number_spectra
    Rcpp::LogicalVector overhead = index <= number_spectra;
    index = index[overhead];

    // remove duplicated indexes
    std::set<int> index_set(index.begin(), index.end());
    std::vector<int> index_unique(index_set.begin(), index_set.end());

    int number_index = index_unique.size();

    // updates index by reference with index_unique
    index = Rcpp::wrap(index_unique);

    std::sort(index.begin(), index.end());

    Rcpp::List list_output(number_index);

    if (number_index > 0) {
      int counter = 0;
      for (int i : index) {
        list_output[counter] = xml_utils::mzxml_parse_binary_data_from_spectrum_node(spectra[i], precision, compression, cols);
        counter++;
      }

      return list_output;
    }

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
  }

  Rcpp::List list_output;
  return list_output;
}




Rcpp::List xml_utils::parse_chromatograms(const pugi::xml_node& root) {

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {

    pugi::xml_node node = root.child("mzML").child("run").child("chromatogramList");

    pugi::xml_node first_chrom = node.child("chromatogram");
    const std::vector<int> precision = xml_utils::mzml_get_precision(first_chrom);
    const std::vector<std::string> compression = xml_utils::mzml_get_compression(first_chrom);
    const Rcpp::CharacterVector cols = xml_utils::mzml_get_binary_names(first_chrom);

    std::vector<pugi::xml_node> chromatograms;

    for (pugi::xml_node child = node.first_child(); child; child = child.next_sibling()) {
      chromatograms.push_back(child);
    }

    int number_chromatograms = chromatograms.size();

    Rcpp::List list_output(number_chromatograms);

    if (number_chromatograms > 0) {
      for (int i = 0; i < number_chromatograms; i++) {
        list_output[i] = xml_utils::mzml_parse_binary_data_from_spectrum_node(chromatograms[i], precision, compression, cols);
      }

      return list_output;
    }

  } else if (strcmp("mzXML", root_name) == 0) {
    std::cout << "\u2717 Format mzXML does not hold chromatograms!" << std::endl;

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
  }

  Rcpp::List list_output;
  return list_output;
}




Rcpp::List xml_utils::parse_partial_chromatograms(
    const pugi::xml_node& root, Rcpp::IntegerVector& index) {

  const char* root_name = root.name();

  if (strcmp("indexedmzML", root_name) == 0) {

    pugi::xml_node node = root.child("mzML").child("run").child("chromatogramList");

    pugi::xml_node first_chrom = node.child("chromatogram");
    const std::vector<int> precision = xml_utils::mzml_get_precision(first_chrom);
    const std::vector<std::string> compression = xml_utils::mzml_get_compression(first_chrom);
    const Rcpp::CharacterVector cols = xml_utils::mzml_get_binary_names(first_chrom);

    std::vector<pugi::xml_node> chromatograms;

    for (pugi::xml_node child = node.first_child(); child; child = child.next_sibling()) {
      chromatograms.push_back(child);
    }

    int number_chromatograms = chromatograms.size();

    // remove any larger than the number_spectra
    Rcpp::LogicalVector overhead = index <= number_chromatograms;
    index = index[overhead];

    // remove duplicated indexes
    std::set<int> index_set(index.begin(), index.end());
    std::vector<int> index_unique(index_set.begin(), index_set.end());
    int number_index = index_unique.size();

    // updates index by reference with index_unique
    index = Rcpp::wrap(index_unique);

    std::sort(index.begin(), index.end());

    Rcpp::List list_output(number_index);

    if (number_index > 0) {
      int counter = 0;
      for (int i : index) {
        list_output[counter] = xml_utils::mzml_parse_binary_data_from_spectrum_node(chromatograms[i], precision, compression, cols);
        counter++;
      }

      return list_output;
    }

  } else if (strcmp("mzXML", root_name) == 0) {
    std::cout << "\u2717 Format mzXML does not hold chromatograms!" << std::endl;

  } else {
    std::cout << "\u2717 The file must be in valid mzML or mzXML format!" << std::endl;
  }

  Rcpp::List list_output;
  return list_output;
}
