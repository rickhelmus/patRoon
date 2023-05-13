#include <fstream>
#include <string>

#include "utils-xml.h"
#include "utils.h"

// startTag: uses partial match (prefix) as tag may contain attributes etc
void parseXMLFile(const char *file, const std::string &startTag,
                  const std::string &endTag,
                  std::function<void(pugi::xml_document &)> func)
{
    std::ifstream ifs(file);
    std::string line, block;
    bool inBlock = false;
    
    while (std::getline(ifs, line))
    {
        trim(line);
        if (strStartsWith(line, startTag))
            inBlock = true;
        
        if (inBlock)
            block += '\n' + line;
        
        if (line == endTag)
        {
            pugi::xml_document doc;
            doc.load_string(block.c_str());
            func(doc);
            inBlock = false;
            block.clear();
        }
    }
}
