#ifdef WITH_MSTK

#include <Rcpp.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "msdata-mstk.h"
#include "spectrum-raw.h"

namespace {

MSToolkit::MSReader getMSTKReader(void)
{
    MSToolkit::MSReader ret;
    ret.addFilter(MSToolkit::MS1);
    ret.addFilter(MSToolkit::MS2);
    return ret;
}

SpectrumRaw getMSTKSpec(MSToolkit::MSReader *msr, const std::string &file, SpectrumRawTypes::Scan scan,
                        const SpectrumRawTypes::MobilityRange &mobRange, SpectrumRawTypes::Intensity minIntensityIMS)
{
    MSToolkit::Spectrum s;
    
    if (!msr->readFile(file.c_str(), s, scan) || s.getScanNumber() == 0)
        Rcpp::stop("Abort: invalid spectrum scan index: %d", scan);
    
    const bool hasMob = s.hasIonMobilityArray();
    SpectrumRaw ret(s.size(), hasMob);
    
    if (hasMob)
    {
        for(size_t i=0; i<ret.size(); ++i)
        {
            if ((!mobRange.isSet() || mobRange.within(s.atIM(i))) && s.at(i).intensity >= minIntensityIMS)
                ret.setPeak(i, s.at(i).mz, s.at(i).intensity, s.atIM(i));
        }
    }
    else
    {
        for(size_t i=0; i<ret.size(); ++i)
            ret.setPeak(i, s.at(i).mz, s.at(i).intensity);
    }
    
    return ret;
}

}


int MSReadBackendMSTK::backends = 0;

// [[Rcpp::interfaces(r, cpp)]]

MSReadBackend::ThreadDataType MSReadBackendMSTK::doGetThreadData(void) const
{
    return std::make_shared<MSToolkit::MSReader>(getMSTKReader());
}

SpectrumRaw MSReadBackendMSTK::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                              const SpectrumRawSelection &scanSel,
                                              const SpectrumRawTypes::MobilityRange &mobRange,
                                              SpectrumRawTypes::Intensity minIntensityIMS) const
{
    auto *msr = reinterpret_cast<MSToolkit::MSReader *>(tdata.get());
    const auto &meta = getSpecMetadata();
    
    if (MSLevel == SpectrumRawTypes::MSLevel::MS1)
        return getMSTKSpec(msr, getCurrentFile(), meta.first.scans[scanSel.index], mobRange, minIntensityIMS);
    if (scanSel.MSMSFrameIndices.empty())
        return getMSTKSpec(msr, getCurrentFile(), meta.second.scans[scanSel.index], mobRange, minIntensityIMS);
    
    // if we are here we need to get MS2 data from an IMS frame...
 
    SpectrumRaw ret;
    for (const auto i : scanSel.MSMSFrameIndices)
        ret.append(getMSTKSpec(msr, getCurrentFile(), meta.second.MSMSFrames[scanSel.index].subScans[i], mobRange,
                               minIntensityIMS));
    
    return ret;
}

void MSReadBackendMSTK::generateSpecMetadata(void)
{
    if (getCurrentFile().empty())
        return;

    auto MSTKReader = getMSTKReader();
    // load first spectrum to get file size and see if we have IMS data
    MSToolkit::Spectrum firstSpec;
    MSTKReader.readFile(getCurrentFile().c_str(), firstSpec, 0);
    const size_t lastScan = MSTKReader.getLastScan();
    const bool haveIMS = firstSpec.hasIonMobilityArray();
    
    if (!haveIMS && firstSpec.getCentroidStatus() != 1)
        Rcpp::stop("Please make sure that file '%s' is centroided!", getCurrentFile().c_str());
    
    if (getNeedIMS() && !haveIMS)
        Rcpp::stop("File '%s' does not contain ion mobility data!", getCurrentFile().c_str());
    
    SpectrumRawMetadata meta;
    
    ThreadExceptionHandler exHandler;
    
    #pragma omp parallel
    {
        auto rd = getMSTKReader();
        MSToolkit::Spectrum spec;
        SpectrumRawMetadata threadMeta;
        
        #pragma omp for schedule(static) nowait
        for (size_t i=1; i<=lastScan; ++i)
        {
            exHandler.run([&]
            {
                rd.readFile(getCurrentFile().c_str(), spec, i);
                const bool isMS1 = spec.getMsLevel() == 1;
                
                SpectrumRawMetadataMS *curMS1MD = (isMS1) ? &threadMeta.first : &threadMeta.second;
                curMS1MD->scans.push_back(spec.getScanNumber());
                curMS1MD->times.push_back(spec.getRTime() * 60); // NOTE: MSTK takes minutes
                curMS1MD->TICs.push_back(spec.getTIC());
                curMS1MD->BPCs.push_back(spec.getBPI());
                curMS1MD->polarities.push_back(spec.getPositiveScan() ? SpectrumRawTypes::MSPolarity::POSITIVE : SpectrumRawTypes::MSPolarity::NEGATIVE);
                if (!isMS1)
                    threadMeta.second.isolationRanges.emplace_back(spec.getSelWindowLower(), spec.getSelWindowUpper());
            });
        }
        
        #pragma omp for schedule(static) ordered
        for(int i=0; i<getOMPNumThreads(); ++i)
        {
            exHandler.run([&]
            {
                #pragma omp ordered
                {
                    meta.first.append(threadMeta.first);
                    meta.second.append(threadMeta.second);
                }
            });
        }
    }

    exHandler.reThrow();
    

    if (haveIMS)
    {
        // For IMS-MS/MS data, the different spectra inside a frame are stored in separate spectra
        // --> move these spectra to MSMSFrames structs, so our main table only contains separate frames
        
        SpectrumRawMetadataMSMS metaMSMSFrames;
        for (size_t i=0; i<meta.second.scans.size(); ++i)
        {
            metaMSMSFrames.scans.push_back(meta.second.scans[i]);
            metaMSMSFrames.times.push_back(meta.second.times[i]);
            metaMSMSFrames.TICs.push_back(meta.second.TICs[i]);
            metaMSMSFrames.BPCs.push_back(meta.second.BPCs[i]);
            metaMSMSFrames.polarities.push_back(meta.second.polarities[i]);
            
            frameMSMSInfo fi;
            const auto curTime = meta.second.times[i];
            bool finished = false;
            while (!finished)
            {
                fi.isolationRanges.push_back(meta.second.isolationRanges[i]);
                fi.subScans.push_back(meta.second.scans[i]);
                ++i;
                if (i > lastScan)
                    break;
                if (!compareTol(curTime, meta.second.times[i]))
                    finished = true; // reached next frame
            }
            metaMSMSFrames.MSMSFrames.push_back(std::move(fi));
            if (finished)
                --i; // we reached one spec to far to find out we're finished
        }
        meta.second = std::move(metaMSMSFrames);
    }

    setSpecMetadata(std::move(meta));
}

#if 0

// the writer was removed from MSTK...
#include "mzMLWriter.h"

void writeMS1SpectraMSTK(std::string path, const std::vector<SpectrumRaw> &spectra, const SpectrumRawMetadataMS &meta)
{
    MSToolkit::MzMLWriter writer;
    
    if (!writer.createMzML(const_cast<char *>(path.data())))
        Rcpp::stop("Failed to create mzML file '%s'!", path.c_str());
    
    if (!writer.createList(true))
        Rcpp::stop("Failed to create spectrum list for '%s'", path.c_str());
    
    for (size_t i=0; i<spectra.size(); ++i)
    {
        MSToolkit::Spectrum MSTKSpec;
        
        MSTKSpec.setRTime(meta.times[i] / 60.0); // NOTE: MSTK writes minutes
        MSTKSpec.setMsLevel(1);
        MSTKSpec.setScanWindow(0.0, 1000.0); // UNDONE: is this important? otherwise we need to add it to the metadata
        
        for (size_t j=0; j<spectra[i].size(); ++j)
            MSTKSpec.add(spectra[i].getMZs()[j], spectra[i].getIntensities()[j]);
        
        if (!writer.writeSpectra(MSTKSpec))
            Rcpp::stop("Failed to write spectrum for file '%s'!", path.c_str());
    }
    
    writer.writeIndex();
    writer.closeList();
    writer.closeMzML();
}
#endif

#endif // WITH_MSTK
