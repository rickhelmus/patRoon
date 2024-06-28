#ifndef PATROON_MSDATA_H
#define PATROON_MSDATA_H

#include <memory>
#include <vector>

#include "spectrum-raw.h"

enum class MSLevel { MS1, MS2 };

struct SpectrumRawMetadata
{
    float time;
    double TIC, BPC, BPMZ;
    MSLevel level;
    std::pair<double, double> isolationRange;
};

class MSReadBackend
{
    std::vector<SpectrumRawMetadata> specMetadata;
    
protected:
    using ThreadDataType = std::shared_ptr<void>;
    
    void addSpecMetadata(const SpectrumRawMetadata &smd) { specMetadata.push_back(smd); }
    void addSpecMetadata(SpectrumRawMetadata &&smd) { specMetadata.emplace_back(std::move(smd)); }
    
public:
    virtual ~MSReadBackend(void) { }
    
    virtual void open(const std::string &file) = 0;
    virtual void close(void) = 0;
    virtual ThreadDataType getThreadData(void) const = 0;
    virtual SpectrumRaw readSpectrum(const ThreadDataType &tdata, int index) const = 0;
    
    const auto &getSpecMetadata(void) const { return specMetadata; }
};

RCPP_EXPOSED_CLASS(MSReadBackend)

#endif
