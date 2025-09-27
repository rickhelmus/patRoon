#ifndef PATROON_MSDATA_HPP
#define PATROON_MSDATA_HPP

#include <functional>
#include <vector>

#include "msdata.h"
#include "spectrum-raw.h"
#include "utils.hpp"

template<typename OutType, typename FuncType, typename... Args>
std::vector<std::vector<OutType>> applyMSData(const MSReadBackend &backend, SpectrumRawTypes::MSLevel MSLevel,
                                              const std::vector<std::vector<SpectrumRawSelection>> &scanSels,
                                              FuncType func, SpectrumRawTypes::Intensity minIntensityIMS,
                                              SpectrumRawTypes::MSSortType IMSSortType,
                                              std::function<SpectrumRaw(const SpectrumRaw &)> prepFunc = {},
                                              Args... args)
{
    /* This function will apply a callback on selected spectra. Multiple sets of spectra selections are supported, and
     * the function is optimized to avoid reading the same spectra more than once in case of overlap between sets.
     * 
     * The code makes the following assumptions:
     * - each set of SpectrumRawSelection objects is unique and sorted by index from low to high
     * - SpectrumRawSelection objects with the same index but different MSMSFrameIndices may occur across sets. Thus,
     *   across sets the full object instead of merely the index should be compared.
     * 
     * As this function is the main driver for eg EICs and peak lists, this function was heavily optimized and therefore
     * more complex.
     */
    
    const auto entries = scanSels.size();
    
    // to compare spectra selections by just the index, i.e. within a selection set
    const auto selCompInd = [](const SpectrumRawSelection &s, size_t i) { return s.index < i; };
    
    // collect all unique indices so we (1) can avoid >1 reads of a spectrum and (2) can easily use an OpenMP for loop below.
    std::vector<size_t> allScanSelInds;
    const size_t nspecs = (MSLevel == SpectrumRawTypes::MSLevel::MS1) ? backend.getSpecMetadata().first.scans.size() :
        backend.getSpecMetadata().second.scans.size();
    for (size_t i=0; i<nspecs; ++i)
    {
        // check if _any_ of the selection sets has this index
        const bool inSels = std::any_of(scanSels.begin(), scanSels.end(), [&](const auto &sels)
        {
            // NOTE: std::binary_search doesn't seem to support comparator with different types (ie the selCompInd lambda)
            // --> use std::lower_bound instead
            const auto it = std::lower_bound(sels.begin(), sels.end(), i, selCompInd);
            return it != sels.end() && it->index == i;
        });
        if (inSels)
            allScanSelInds.push_back(i);
    }
    
    // pre-allocate ret vector
    std::vector<std::vector<OutType>> ret;
    for (size_t i=0; i<entries; ++i)
        ret.emplace_back(scanSels[i].size());
    
    ThreadExceptionHandler exHandler;
    
    #pragma omp parallel
    {
        auto tdata = backend.getThreadData();
        #pragma omp for
        for (const auto selInd : allScanSelInds)
        {
            // UNDONE: handle mobility range? For now we stick with post-filtering, which makes the mob range argument
            // to readSpectrum currently unused...
            
            /* In case multiple sets contain spectrum selection with index selInd, there is a chance that the
             * SpectrumRawSelection objects still differ (see above). We therefore keep track of the selection of the
             * previous set and if it differs perform a new spectrum read.
             * 
             * UNDONE: this could be optimized further by also storing and re-using spectra before the previous selection set
             */ 
            
            exHandler.run([&] {
                SpectrumRaw curSpec;
                SpectrumRawSelection curSel;
                bool initSel = true;
                
                for (size_t j=0; j<entries; ++j)
                {
                    auto it = std::lower_bound(scanSels[j].begin(), scanSels[j].end(), selInd, selCompInd);
                    while (it != scanSels[j].end() && it->index == selInd)
                    {
                        if (initSel || curSel != *it)
                        {
                            initSel = false;
                            curSel = *it;
                            curSpec = backend.readSpectrum(tdata, MSLevel, curSel, SpectrumRawTypes::MobilityRange(),
                                                           minIntensityIMS);
                            if (curSpec.hasMobilities()) // NOTE: for non-IMS data we assume data is always mz sorted
                                curSpec.sort(IMSSortType);
                            if (prepFunc)
                                curSpec = prepFunc(curSpec);
                        }
                        
                        // UNDONE: optimization could be further pushed by not using a 2d vector here?
                        const auto ind = std::distance(scanSels[j].begin(), it);
                        ret[j][ind] = func(curSpec, curSel, j, args...);
                        ++it;
                    }
                }
            });
        }
    }
    
    exHandler.reThrow();
    
    return ret;
}

template <typename T, typename I> std::vector<I> fillEIXIntensities(const std::vector<T> &allXValues,
                                                                    const std::vector<T> &xvalues,
                                                                    const std::vector<I> &intensities)
{
    if (allXValues.empty() || xvalues.empty())
        return std::vector<I>(allXValues.size(), 0);
    
    std::vector<I> intens(allXValues.size());
    auto it = xvalues.begin();
    auto allIt = std::lower_bound(allXValues.begin(), allXValues.end(), *it);
    while (it != xvalues.end() && allIt != allXValues.end())
    {
        const auto ind = std::distance(allXValues.begin(), allIt);
        intens[ind] = intensities[std::distance(xvalues.begin(), it)];
        ++it;
        if (it != xvalues.end())
            allIt = std::lower_bound(allIt, allXValues.end(), *it);
    }
    
    return intens;
}

#endif
