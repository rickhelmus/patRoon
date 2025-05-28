/*
 * SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <fstream>
#include <string>

#include <Rcpp.h>

#include "utils-xml.h"
#include "utils.h"

// startTag: uses partial match (prefix) as tag may contain attributes etc
void parseXMLFile(const char *file, const std::string &startTag,
                  const std::string &endTag,
                  std::function<void(pugi::xml_document &)> func)
{
    std::ifstream ifs(file);
    std::string line, block;
    unsigned tagCounter = 0; // use a counter to supported nested blocks
    
    while (std::getline(ifs, line))
    {
        trim(line);
        if (strStartsWith(line, startTag))
            ++tagCounter;
        
        if (tagCounter != 0)
            block += '\n' + line;
        
        if (line == endTag)
        {
            --tagCounter;
            if (tagCounter == 0)
            {
                pugi::xml_document doc;
                doc.load_string(block.c_str());
                func(doc);
                block.clear();
            }
        }
    }
}

// [[Rcpp::export]]
void addFilesToOpenMSIni(const std::string &file, const std::vector<std::string> &inFiles,
                         const std::vector<std::string> &outFiles)
{
    pugi::xml_document doc;
    const auto result = doc.load_file(file.c_str());
    
    if (!result)
        Rcpp::stop("Failed to parse XML file ('%s'): %s", file, result.description());
    
    auto nodeEl = doc.child("PARAMETERS").child("NODE").child("NODE");

    const auto addFiles = [&](const std::vector<std::string> &files, const char *itemName)
    {
        auto el = nodeEl.find_child_by_attribute("ITEMLIST", "name", itemName);
        for (auto &f : files)
        {
            auto fileElAttr = el.append_child("LISTITEM").append_attribute("value");
            fileElAttr.set_value(f.c_str());
        }
    };
        
    if (!inFiles.empty())
        addFiles(inFiles, "in");
    if (!outFiles.empty())
        addFiles(outFiles, "out");
    
    doc.save_file(file.c_str());
}
