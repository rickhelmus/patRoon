#ifndef PATROON_MSDATA_H
#define PATROON_MSDATA_H

#include <memory>
#include <vector>

#include "spectrum-raw.h"

class MSReadBackend
{
protected:
    using ThreadDataType = std::shared_ptr<void>;

private:
    std::string currentFile;
    mutable SpectrumRawMetadata specMetadata; // mutable to allow lazy loading
    
    virtual void doOpen(const std::string &file) = 0;
    virtual void doClose(void) = 0;
    virtual ThreadDataType doGetThreadData(void) const = 0;
    virtual SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::Scan scan) const = 0;
    
protected:
    void emplaceSpecMeta(SpectrumRawMetadata &&msd) const { specMetadata = std::move(msd); }
    
    // NOTE: not private to allow parent method calls
    virtual const SpectrumRawMetadata &doGetSpectrumRawMetadata(void) const { return specMetadata; }
    
public:
    virtual ~MSReadBackend(void) { close(); }
    
    void open(const std::string &file);
    void close(void);
    const std::string &getCurrentFile(void) const { return currentFile; }

    ThreadDataType getThreadData(void) const { return doGetThreadData(); }
    SpectrumRaw readSpectrum(const ThreadDataType &tdata, int index) const { return doReadSpectrum(tdata, index); };
    const SpectrumRawMetadata &getSpecMetadata(void) const { return doGetSpectrumRawMetadata(); }
};

RCPP_EXPOSED_CLASS(MSReadBackend)

#endif
