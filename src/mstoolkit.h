#ifndef PATROON_MSTOOLKIT_H
#define PATROON_MSTOOLKIT_H

#include <memory>
#include <string>

#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

class MSToolkitBackend
{
public:
    using ThreadDataType = std::unique_ptr<MSToolkit::MSReader>;
    
private:
    std::string currentFile;
    
public:
    MSToolkitBackend(void) { };
    
    void open(const std::string &file) { if (!currentFile.empty()) close(); currentFile = file; };
    void close(void) { };
    
    ThreadDataType getThreadData(void); 
};

#endif
