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
    SpectrumRawMetadata specMetadata;
    
    virtual void doOpen(const std::string &file) = 0;
    virtual void doClose(void) = 0;
    virtual ThreadDataType doGetThreadData(void) const = 0;
    virtual SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                       const SpectrumRawSelection &scanSel,
                                       const SpectrumRawTypes::MobilityRange &mobRange,
                                       SpectrumRawTypes::Intensity minIntensityIMS) const = 0;
    
public:
    MSReadBackend(void) = default;
    MSReadBackend(const MSReadBackend &) = delete;
    virtual ~MSReadBackend(void) { }
    
    void open(const std::string &file);
    void close(void);
    const std::string &getCurrentFile(void) const { return currentFile; }

    ThreadDataType getThreadData(void) const { return doGetThreadData(); }
    SpectrumRaw readSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                             const SpectrumRawSelection &scanSel,
                             const SpectrumRawTypes::MobilityRange &mobRange,
                             SpectrumRawTypes::Intensity minIntensityIMS) const
        { return doReadSpectrum(tdata, MSLevel, scanSel, mobRange, minIntensityIMS); };
    const SpectrumRawMetadata &getSpecMetadata(void) const { return specMetadata; }
    void setSpecMetadata(SpectrumRawMetadata &&msd) { specMetadata = std::move(msd); }
};

RCPP_EXPOSED_CLASS(MSReadBackend)

class MSReadBackendUnavailable: MSReadBackend
{
    void doOpen(const std::string &) override { }
    void doClose(void) override { }
    ThreadDataType doGetThreadData(void) const override { return nullptr; };
    SpectrumRaw doReadSpectrum(const ThreadDataType &, SpectrumRawTypes::MSLevel, const SpectrumRawSelection &,
                               const SpectrumRawTypes::MobilityRange &,
                               SpectrumRawTypes::Intensity) const override { return SpectrumRaw(); };
    
public:
    MSReadBackendUnavailable(const char *n);
    ~MSReadBackendUnavailable(void) { }
};

#endif