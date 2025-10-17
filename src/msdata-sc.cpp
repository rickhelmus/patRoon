/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <Rcpp.h>

#include "msdata-sc.h"
#include "spectrum-raw.h"

#include "utils-xml.h" // for pugixml
#include "StreamCraft/StreamCraft_lib.cpp"

namespace {

SpectrumRaw getSCSpectrum(sc::MS_FILE *analysis, SpectrumRawTypes::Scan scan,
                          const SpectrumRawTypes::MobilityRange &mobRange, SpectrumRawTypes::Intensity minIntensityIMS)
{
    const auto s = analysis->get_spectrum(scan);
    
    const auto mobArrayIt = std::find(s.binary_names.cbegin(), s.binary_names.cend(), "ion_mobility");
    const bool hasMobArray = mobArrayIt != s.binary_names.cend();
    
    
    // UNDONE: always safe to assume that m/z / intensity are in first two arrays?
    
    if (hasMobArray)
    {
        const size_t mobArrayInd = std::distance(s.binary_names.begin(), mobArrayIt);
        
        // For IMS we don't know yet the size if intensity/mob filters are set: only preallocate if intensity == 0
        
        if (!mobRange.isSet() && minIntensityIMS == 0)
        {
            SpectrumRaw ret(s.array_length, true);
            for (int i=0; i<s.array_length; ++i)
                ret.setPeak(i, s.binary_data[0][i], s.binary_data[1][i], s.binary_data[mobArrayInd][i]);
            return ret;
        }
        else
        {
            SpectrumRaw ret;
            for (size_t i=0; i<s.array_length; ++i)
            {
                const auto mob = s.binary_data[mobArrayInd][i];
                if ((!mobRange.isSet() || mobRange.within(mob)) && s.binary_data[1][i] >= minIntensityIMS)
                    ret.append(s.binary_data[0][i], s.binary_data[1][i], mob);
            }
            return ret;
        }
    }
    
    // No IMS or no IMS array
    SpectrumRaw ret(s.array_length, false);
    for (int i=0; i<s.array_length; ++i)
        ret.setPeak(i, s.binary_data[0][i], s.binary_data[1][i]);
    
    if (s.mobility != 0)
    {
        // IMS spectrum w/out array, eg MS2 produced by TIMSCONVERT
        if (mobRange.isSet() && !mobRange.within(s.mobility))
            return SpectrumRaw(); // UNDONE?
        // HACK: to simplify things, we just set the mobility of all peaks to the same value
        ret.setAllMobilities(s.mobility);
    }
    
    return ret;
}

}

// [[Rcpp::interfaces(r, cpp)]]

void MSReadBackendSC::doOpen(const std::string &file)
{
    handle = std::make_unique<sc::MS_FILE>(file);
}

void MSReadBackendSC::doClose(void)
{
    handle.reset();
}

SpectrumRaw MSReadBackendSC::doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                                            const SpectrumRawSelection &scanSel,
                                            const SpectrumRawTypes::MobilityRange &mobRange,
                                            SpectrumRawTypes::Intensity minIntensityIMS) const
{
    const auto &meta = getSpecMetadata();
    
    if (MSLevel == SpectrumRawTypes::MSLevel::MS1)
        return getSCSpectrum(handle.get(), meta.first.scans[scanSel.index], mobRange, minIntensityIMS);
    if (scanSel.MSMSFrameIndices.empty())
        return getSCSpectrum(handle.get(), meta.second.scans[scanSel.index], mobRange, minIntensityIMS);
    
    // if we are here we need to get MS2 data from an IMS frame...
    
    SpectrumRaw ret;
    for (const auto i : scanSel.MSMSFrameIndices)
        ret.append(getSCSpectrum(handle.get(), meta.second.MSMSFrames[scanSel.index].subScans[i], mobRange,
                                 minIntensityIMS));
    
    return ret;
}

void MSReadBackendSC::generateSpecMetadata(void)
{
    if (getCurrentFile().empty())
        return;

    sc::MS_FILE analysis(getCurrentFile());
    
    const auto hd = analysis.get_spectra_headers();
    
    const bool hasMob = analysis.has_ion_mobility();

    if (!hasMob && analysis.get_spectra_mode({0})[0] != 2)
        Rcpp::stop("Please make sure that file '%s' is centroided!", getCurrentFile().c_str());
    
    if (getNeedIMS() && !hasMob)
        Rcpp::stop("File '%s' does not contain ion mobility data!", getCurrentFile().c_str());
    
    SpectrumRawMetadata meta;
    for (size_t i=0; i<hd.index.size(); ++i)
    {
        const bool isMS1 = hd.level[i] == 1;
        
        SpectrumRawMetadataMS *curMS1MD = (isMS1) ? &meta.first : &meta.second;
        curMS1MD->scans.push_back(hd.index[i]);
        curMS1MD->times.push_back(hd.rt[i]);
        curMS1MD->TICs.push_back(hd.tic[i]);
        curMS1MD->BPCs.push_back(hd.bpint[i]);
        // NOTE: a simple cast is sufficient as the enum values match the int values from SC
        curMS1MD->polarities.push_back(static_cast<SpectrumRawTypes::MSPolarity>(hd.polarity[i]));
        
        if (!isMS1)
        {
            if (!hasMob)
            {
                meta.second.isolationRanges.emplace_back(hd.window_mzlow[i], hd.window_mzhigh[i]);
                meta.second.precursorMZs.push_back(hd.precursor_mz[i]);
            }
            else
            {
                // For IMS-MS/MS data, the different spectra inside a frame are stored in separate spectra
                // --> move these spectra to MSMSFrames structs, so our main table only contains separate frames
                frameMSMSInfo fi;
                const auto curTime = hd.rt[i];
                bool finished = false;
                while (!finished)
                {
                    fi.isolationRanges.emplace_back(hd.window_mzlow[i], hd.window_mzhigh[i]);
                    fi.precursorMZs.push_back(hd.precursor_mz[i]);
                    fi.subScans.push_back(hd.index[i]);
                    ++i;
                    if (i >= hd.index.size())
                        break;
                    if (!compareTol(curTime, hd.rt[i]))
                        finished = true; // reached next frame
                }
                meta.second.MSMSFrames.push_back(std::move(fi));
                if (finished)
                    --i; // we reached one spec to far to find out we're finished
            }
        }
    }
    
    setSpecMetadata(std::move(meta));
}
