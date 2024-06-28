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
protected:
    using ThreadDataType = std::shared_ptr<void>;
    
private:
    std::string currentFile;
    std::vector<SpectrumRawMetadata> specMetadata;
    
    virtual void doOpen(const std::string &file) = 0;
    virtual void doClose(void) = 0;
    virtual ThreadDataType doGetThreadData(void) const = 0;
    virtual SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, int index) const = 0;
    
    
protected:
    void addSpecMetadata(const SpectrumRawMetadata &smd) { specMetadata.push_back(smd); }
    void addSpecMetadata(SpectrumRawMetadata &&smd) { specMetadata.emplace_back(std::move(smd)); }
    
public:
    virtual ~MSReadBackend(void) { close(); }
    
    void open(const std::string &file);
    void close(void);
    const std::string &getCurrentFile(void) const { return currentFile; }
    
    ThreadDataType getThreadData(void) const { return doGetThreadData(); }
    SpectrumRaw readSpectrum(const ThreadDataType &tdata, int index) const { return doReadSpectrum(tdata, index); };
    const auto &getSpecMetadata(void) const { return specMetadata; }
};

RCPP_EXPOSED_CLASS(MSReadBackend)

#endif
