#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "opentims++/opentims_all.h"

#include <Rcpp.h>

#include "utils.h"

namespace {

struct SpectrumIMS // UNDONE: merge with other struct(s)
{
    std::vector<uint32_t> IDs, intensities;
    std::vector<double> mzs, mobilities;
    
    SpectrumIMS(size_t size) : IDs(size), intensities(size), mzs(size), mobilities(size) { }
    SpectrumIMS() = default;
    
    void addData(unsigned id, double mz, unsigned inten, double mob)
    {
        IDs.push_back(id);
        mzs.push_back(mz);
        intensities.push_back(inten);
        mobilities.push_back(mob);
    }
    void addData(const SpectrumIMS &sp)
    {
        IDs.insert(IDs.end(), sp.IDs.begin(), sp.IDs.end());
        mzs.insert(mzs.end(), sp.mzs.begin(), sp.mzs.end());
        intensities.insert(intensities.end(), sp.intensities.begin(), sp.intensities.end());
        mobilities.insert(mobilities.end(), sp.mobilities.begin(), sp.mobilities.end());
    }
    void resize(size_t s)
    {
        IDs.resize(s);
        mzs.resize(s);
        intensities.resize(s);
        mobilities.resize(s);
    }
    size_t size(void) const { return IDs.size(); }
    bool empty(void) const {return IDs.empty(); }
};

auto getTIMSDecompBuffers(size_t size)
{
    // get buffers for decompression (see TimsDataHandle::extract_frames())
    return std::make_pair(std::unique_ptr<ZSTD_DCtx, decltype(&ZSTD_freeDCtx)>(ZSTD_createDCtx(), &ZSTD_freeDCtx),
                          std::make_unique<char[]>(size));
}

// From http://mochan.info/c++/2019/06/12/returntype-deduction.html
using TIMSDecompBufferPair = decltype(getTIMSDecompBuffers(std::declval<size_t>()));


SpectrumIMS getIMSFrame(TimsFrame &frame)
{
    SpectrumIMS spec(frame.num_peaks);
    frame.save_to_buffs(nullptr, spec.IDs.data(), nullptr, spec.intensities.data(), spec.mzs.data(),
                        spec.mobilities.data(), nullptr);
    return spec;
}

SpectrumIMS getIMSFrame(TimsFrame &frame, const TIMSDecompBufferPair &bufs)
{
    frame.decompress(bufs.second.get(), bufs.first.get());
    SpectrumIMS ret = getIMSFrame(frame);
    frame.close();
    return ret;
}

SpectrumIMS filterSpectrum(const SpectrumIMS &spec, unsigned scanBegin = 0, unsigned scanEnd = 0,
                           double mzStart = 0.0, double mzEnd = 0.0, double mobilityStart = 0.0,
                           double mobilityEnd = 0.0)
{
    if (scanBegin == 0 && scanEnd == 0 && mzStart == 0.0 && mzEnd == 0.0 && mobilityStart == 0.0 &&
        mobilityEnd == 0.0)
        return spec;
    
    SpectrumIMS specFiltered;
    for (size_t i=0; i<spec.size(); ++i)
    {
        if (scanBegin != 0 && spec.IDs[i] < scanBegin)
            continue;
        if (scanEnd != 0 && spec.IDs[i] > scanEnd)
            continue;
        if (mzStart != 0.0 && spec.mzs[i] < mzStart)
            continue;
        if (mzEnd != 0.0 && spec.mzs[i] > mzEnd)
            continue;
        if (mobilityStart != 0.0 && spec.mobilities[i] < mobilityStart)
            continue;
        if (mobilityEnd != 0.0 && spec.mobilities[i] > mobilityEnd)
            continue;
        
        specFiltered.addData(spec.IDs[i], spec.mzs[i], spec.intensities[i], spec.mobilities[i]);
    }
    return specFiltered;
}

SpectrumIMS collapseIMSFrame(const SpectrumIMS &frame, clusterMethod method, double mzWindow)
{
    const std::vector<int> clusts = clusterNums(frame.mzs, method, mzWindow);
    const int maxClust = *(std::max_element(clusts.begin(), clusts.end()));
    SpectrumIMS binnedSpectrum(maxClust + 1);
    std::vector<unsigned> binSizes(maxClust + 1);
    
    // assign unique IDs
    std::iota(binnedSpectrum.IDs.begin(), binnedSpectrum.IDs.end(), 1);
    
    // sum data for each cluster
    for (size_t i=0; i<clusts.size(); ++i)
    {
        const size_t cl = clusts[i];
        binnedSpectrum.mzs[cl] += frame.mzs[i];
        binnedSpectrum.intensities[cl] += frame.intensities[i];
        binnedSpectrum.mobilities[cl] += frame.mobilities[i];
        ++binSizes[cl];
    }

    // average data
    for (size_t i=0; i<binnedSpectrum.size(); ++i)
    {
        const double len = static_cast<double>(binSizes[i]);
        binnedSpectrum.mzs[i] /= len;
        binnedSpectrum.mobilities[i] /= len;
    }

    return binnedSpectrum;
}

SpectrumIMS collapseIMSFrames(TimsDataHandle &TDH, const std::vector<unsigned> &frameIDs,
                              const std::vector<unsigned> &scanStarts, const std::vector<unsigned> &scanEnds,
                              double mobilityStart, double mobilityEnd, clusterMethod method, double mzWindow)
{
    SpectrumIMS sumSpec;
    
    // UNDONE: make num_threads configurable
    #pragma omp parallel num_threads(8)
    {
        auto TBuffers = getTIMSDecompBuffers(TDH.get_decomp_buffer_size());
        SpectrumIMS threadSumSpec;
    
        #pragma omp for nowait
        for (size_t i=0; i<frameIDs.size(); ++i)
        {
            const auto fri = frameIDs[i];
            if (!TDH.has_frame(fri))
                continue;
            auto &fr = TDH.get_frame(fri);
            const SpectrumIMS spec = getIMSFrame(fr, TBuffers);
            const auto specF = filterSpectrum(spec, scanStarts[i], scanEnds[i], 0.0, 0.0, mobilityStart, mobilityEnd);
            threadSumSpec.addData(collapseIMSFrame(specF, method, mzWindow));
        }
        
        #pragma omp critical
        {
            sumSpec.addData(threadSumSpec);
        }
    }
    
    
    // collapse result
    if (!sumSpec.empty())
    {
        sumSpec = collapseIMSFrame(sumSpec, method, mzWindow);
        // average intensities
        for (auto &inten : sumSpec.intensities)
            inten /= frameIDs.size();
    }
    
    return sumSpec;
}

}

// [[Rcpp::export]]
void initBrukerLibrary(const std::string &path)
{
    setup_bruker(path);
}

// [[Rcpp::export]]
Rcpp::DataFrame collapseTIMSFrame(const std::string &file, size_t frameID, const std::string &method, double mzWindow)
{
    TimsDataHandle TDH(file);
    if (!TDH.has_frame(frameID))
        Rcpp::stop("Frame doesn't exist.");
    
    const SpectrumIMS frame = getIMSFrame(TDH.get_frame(frameID));
    const auto spec = collapseIMSFrame(frame, clustMethodFromStr(method), mzWindow);
    
    return Rcpp::DataFrame::create(Rcpp::Named("ID") = spec.IDs,
                                   Rcpp::Named("mz") = spec.mzs,
                                   Rcpp::Named("intensity") = spec.intensities,
                                   Rcpp::Named("mobility") = spec.mobilities);
}

// [[Rcpp::export]]
Rcpp::List getTIMSPeakLists(const std::string &file, Rcpp::List frameIDsList, Rcpp::List scanStartsList,
                            Rcpp::List scanEndsList, const std::vector<double> &mobilityStarts,
                            const std::vector<double> &mobilityEnds, const std::string &method, double mzWindow)
{
    const auto count = frameIDsList.size();
    TimsDataHandle TDH(file);
    std::vector<SpectrumIMS> peakLists(count);
    const auto clMethod = clustMethodFromStr(method);
    
    for (int i=0; i<count; ++i)
    {
        const auto frameIDs = Rcpp::as<std::vector<unsigned>>(frameIDsList[i]);
        const auto scanStarts = Rcpp::as<std::vector<unsigned>>(scanStartsList[i]);
        const auto scanEnds = Rcpp::as<std::vector<unsigned>>(scanEndsList[i]);
        peakLists[i] = std::move(collapseIMSFrames(TDH, frameIDs, scanStarts, scanEnds, mobilityStarts[i],
                                                   mobilityEnds[i], clMethod, mzWindow));
    }
    
    Rcpp::List ret(peakLists.size());
    for (size_t i=0; i<peakLists.size(); ++i)
        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("ID") = peakLists[i].IDs,
                                         Rcpp::Named("mz") = peakLists[i].mzs,
                                         Rcpp::Named("intensity") = peakLists[i].intensities,
                                         Rcpp::Named("mobility") = peakLists[i].mobilities);
    return ret;
}

// [[Rcpp::export]]
Rcpp::List getTIMSEIC(const std::string &file, const std::vector<unsigned> &frameIDs,
                      std::vector<double> mzStarts, std::vector<double> mzEnds, std::vector<double> mobilityStarts,
                      std::vector<double> mobilityEnds)
{
    TimsDataHandle TDH(file);
    struct EIC // UNDONE?
    {
        std::vector<double> times;
        std::vector<uint32_t> intensities;
    };
    std::vector<EIC> EICs(mzStarts.size());
    
    // UNDONE: make num_threads configurable
    #pragma omp parallel num_threads(8)
    {
        auto TBuffers = getTIMSDecompBuffers(TDH.get_decomp_buffer_size());
        std::vector<EIC> threadEICs(mzStarts.size());
        
        #pragma omp for nowait
        for (size_t i=0; i<frameIDs.size(); ++i)
        {
            auto &fr = TDH.get_frame(frameIDs[i]);
            if (fr.msms_type != 0)
                continue; // UNDONE: needed?
            const SpectrumIMS spec = getIMSFrame(fr, TBuffers);
            
            for (size_t j=0; j<EICs.size(); ++j)
            {
                const SpectrumIMS frameF = filterSpectrum(spec, 0, 0, mzStarts[j], mzEnds[j], mobilityStarts[j],
                                                          mobilityEnds[j]);
                threadEICs[j].times.push_back(fr.time);
                threadEICs[j].intensities.push_back(std::accumulate(frameF.intensities.begin(),
                                                                    frameF.intensities.end(), 0));
            }
        }
        
        #pragma omp critical
        {
            for (size_t i=0; i<EICs.size(); ++i)
            {
                EICs[i].times.insert(EICs[i].times.end(), threadEICs[i].times.begin(), threadEICs[i].times.end());
                EICs[i].intensities.insert(EICs[i].intensities.end(), threadEICs[i].intensities.begin(),
                                           threadEICs[i].intensities.end());
            }
        }
    }
    
    // sort and compress
    for (size_t i=0; i<EICs.size(); ++i)
    {
        const auto ord = getSortedInds(EICs[i].times);
        const auto len = ord.size();
        EIC ef;
        
        for (size_t j=0; j<len; ++j)
        {
            const auto ordInd = ord[j];
            if (j > 0 && j < (len-1) && EICs[i].intensities[ordInd] == 0)
            {
                const auto prevInt = EICs[i].intensities[ord[j-1]], nextInt = EICs[i].intensities[ord[j+1]];
                if (prevInt == 0 && nextInt == 0)
                    continue; // skip points with zero intensities that are neighbored by others.
                
            }
            ef.times.push_back(EICs[i].times[ordInd]);
            ef.intensities.push_back(EICs[i].intensities[ordInd]);
        }
        EICs[i] = std::move(ef);
    }
    
    Rcpp::List ret(EICs.size());
    for (size_t i=0; i<EICs.size(); ++i)
        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("time") = EICs[i].times,
                                         Rcpp::Named("intensity") = EICs[i].intensities);
    return ret;
}
