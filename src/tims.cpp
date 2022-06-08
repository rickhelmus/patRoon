#include <cstdint>
#include <vector>

#include "opentims++/opentims_all.h"

#include <Rcpp.h>

#include "utils.h"

namespace {

struct SpectrumIMS // UNDONE: merge with other struct(s)
{
    std::vector<unsigned> IDs, intensities;
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
    void resize(size_t s)
    {
        IDs.resize(s);
        mzs.resize(s);
        intensities.resize(s);
        mobilities.resize(s);
    }
    size_t size(void) const { return IDs.size(); }
};

    
SpectrumIMS getIMSFrame(TimsFrame &frame, unsigned scanBegin = 0, unsigned scanEnd = 0,
                        double mzStart = 0.0, double mzEnd = 0.0, double mobilityStart = 0.0,
                        double mobilityEnd = 0.0)
{
    SpectrumIMS spec(frame.num_peaks);
    frame.save_to_buffs(nullptr, spec.IDs.data(), nullptr, spec.intensities.data(), spec.mzs.data(),
                        spec.mobilities.data(), nullptr);
    
    if (scanBegin == 0 && scanEnd == 0 && mzStart == 0.0 && mzEnd == 0.0 && mobilityStart == 0.0 && mobilityEnd == 0.0)
        return(spec); // no need to filter
    
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

}

// [[Rcpp::export]]
void initBrukerLibrary(const std::string &path)
{
    setup_bruker(path);
}

// [[Rcpp::export]]
Rcpp::List getTIMSFrame(const std::string &file, size_t frameID)
{
    TimsDataHandle TDH(file);
    if (TDH.has_frame(frameID))
    {
        TimsFrame &frame = TDH.get_frame(frameID);
        auto frameIDs = std::vector<uint32_t>(frame.num_peaks);
        auto scanIDs = std::vector<uint32_t>(frame.num_peaks);
        auto TOFs = std::vector<uint32_t>(frame.num_peaks);
        auto intensities = std::vector<uint32_t>(frame.num_peaks);
        auto mzs = std::vector<double>(frame.num_peaks);
        auto invIonMobilities = std::vector<double>(frame.num_peaks);
        auto retentionTimes = std::vector<double>(frame.num_peaks);
        frame.save_to_buffs(frameIDs.data(), scanIDs.data(), TOFs.data(), intensities.data(), mzs.data(),
                            invIonMobilities.data(), retentionTimes.data());
        return Rcpp::List::create(Rcpp::Named("frame") = frameIDs,
                                  Rcpp::Named("scan") = scanIDs,
                                  Rcpp::Named("tof") = TOFs,
                                  Rcpp::Named("intensity") = intensities,
                                  Rcpp::Named("mz") = mzs,
                                  Rcpp::Named("inv_ion_mobility") = invIonMobilities,
                                  Rcpp::Named("retention_time") = retentionTimes);
    }
    return Rcpp::List();
}

// [[Rcpp::export]]
Rcpp::List collapseTIMSFrames(const std::string &file, double retMin, double retMax)
{
    /// UNDONE: double comparisons?
    
    TimsDataHandle TDH(file);
    auto &frameDescs = TDH.get_frame_descs();
    auto it = std::lower_bound(frameDescs.begin(), frameDescs.end(), retMin,
                               [&](const auto &pair, double value) { return pair.second.time < value; });
    
    struct MSFrame
    {
        uint32_t ID;
        double retentionTime;
        std::vector<uint32_t> scanIDs, intensities;
        std::vector<double> mzs, invIonMobilities;
        MSFrame(uint32_t i, double rt, const size_t peaks) : ID(i), retentionTime(rt), scanIDs(peaks),
                                                             intensities(peaks), mzs(peaks), invIonMobilities(peaks) {};
    };
    
    std::vector<MSFrame> MSFrames;
    for (; it != frameDescs.end() && it->second.time < retMax; ++it)
    {
        MSFrame frame(it->second.id, it->second.time, it->second.num_peaks);
        it->second.save_to_buffs(nullptr, frame.scanIDs.data(), nullptr, frame.intensities.data(), frame.mzs.data(),
                                 frame.invIonMobilities.data(), nullptr);
        MSFrames.push_back(std::move(frame));
    }
    
    return Rcpp::List();
}

// [[Rcpp::export]]
Rcpp::List clusterTIMSFrame(const std::string &file, size_t frameID, const std::string &method, double mzWindow)
{
    TimsDataHandle TDH(file);
    if (!TDH.has_frame(frameID))
        Rcpp::stop("Frame doesn't exist.");
    
    SpectrumIMS frame = getIMSFrame(TDH.get_frame(frameID));
    
    clusterMethod clMethod;
    if (method == "bin")
        clMethod = clusterMethod::BIN;
    else if (method == "diff")
        clMethod = clusterMethod::DIFF;
    else // if (method == "hclust")
        clMethod = clusterMethod::HCLUST;
    
    const std::vector<int> clusts = clusterNums(frame.mzs, clMethod, mzWindow);
    
    return Rcpp::List::create(Rcpp::Named("ID") = frame.IDs,
                              Rcpp::Named("mz") = frame.mzs,
                              Rcpp::Named("intensity") = frame.intensities,
                              Rcpp::Named("mobility") = frame.mobilities,
                              Rcpp::Named("clust") = clusts);
}

// [[Rcpp::export]]
Rcpp::List clusterTIMSFrame2(const std::string &file, size_t frameID, const std::string &method, double mzWindow)
{
    TimsDataHandle TDH(file);
    if (!TDH.has_frame(frameID))
        Rcpp::stop("Frame doesn't exist.");
    
    SpectrumIMS frame = getIMSFrame(TDH.get_frame(frameID));
    
    clusterMethod clMethod;
    if (method == "bin")
        clMethod = clusterMethod::BIN;
    else if (method == "diff")
        clMethod = clusterMethod::DIFF;
    else // if (method == "hclust")
        clMethod = clusterMethod::HCLUST;
    
    auto spec = collapseIMSFrame(frame, clMethod, mzWindow);
    
    return Rcpp::List::create(Rcpp::Named("ID") = spec.IDs,
                              Rcpp::Named("mz") = spec.mzs,
                              Rcpp::Named("intensity") = spec.intensities,
                              Rcpp::Named("mobility") = spec.mobilities);
}
