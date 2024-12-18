/*
 * SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#ifndef PATROON_PARSE_XML_H
#define PATROON_PARSE_XML_H

#define PUGIXML_HEADER_ONLY
#define PUGIXML_NO_XPATH
#include "pugixml/pugixml.hpp"

#include <functional>

typedef double numType;
template<typename T> numType getNumericFromXML(const T &x) { return x.as_double(); }

void parseXMLFile(const char *file, const std::string &startTag,
                  const std::string &endTag,
                  std::function<void(pugi::xml_document &)> func);
    
#endif