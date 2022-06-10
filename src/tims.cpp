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
};

    
SpectrumIMS getIMSFrame(TimsFrame &frame)
{
    SpectrumIMS spec(frame.num_peaks);
    frame.save_to_buffs(nullptr, spec.IDs.data(), nullptr, spec.intensities.data(), spec.mzs.data(),
                        spec.mobilities.data(), nullptr);
    return spec;
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

SpectrumIMS collapseIMSFrames(TimsDataHandle &TDH, const std::vector<unsigned> &frameIDs, clusterMethod method,
                              double mzWindow)
{
    SpectrumIMS ret;
    int specN = 0;
    
    for (auto i : frameIDs)
    {
        if (!TDH.has_frame(i))
            continue;
        auto &fr = TDH.get_frame(i);
        const SpectrumIMS spec = getIMSFrame(fr); // UNDONE: filter args
        ret.addData(collapseIMSFrame(spec, method, mzWindow));
        ++specN;
    }
    
    // collapse result
    ret = collapseIMSFrame(ret, method, mzWindow);
    // average intensities
    for (auto &inten : ret.intensities)
        inten /= specN;
    
    return ret;
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
Rcpp::DataFrame collapseTIMSFrames(const std::string &file, const std::vector<unsigned> &frameIDs,
                                   const std::string &method, double mzWindow)
{
    TimsDataHandle TDH(file);
    const auto spec = collapseIMSFrames(TDH, frameIDs, clustMethodFromStr(method), mzWindow);
    return Rcpp::DataFrame::create(Rcpp::Named("ID") = spec.IDs,
                                   Rcpp::Named("mz") = spec.mzs,
                                   Rcpp::Named("intensity") = spec.intensities,
                                   Rcpp::Named("mobility") = spec.mobilities);
}

// [[Rcpp::export]]
Rcpp::List getTIMSEIC(const std::string &file, const std::vector<unsigned> &frameIDs,
                           std::vector<double> mzStarts, std::vector<double> mzEnds, std::vector<double> mobilityStarts,
                           std::vector<double> mobilityEnds)
{
    TimsDataHandle TDH(file);
    struct EIC // UNDONE?
    {
        std::vector<double> times, intensities;
    };
    std::vector<EIC> EICs(mzStarts.size());

    for (auto i : frameIDs)
    {
        auto &fr = TDH.get_frame(i);
        if (fr.msms_type != 0)
            continue; // UNDONE?
        const SpectrumIMS spec = getIMSFrame(fr);
        for (size_t j=0; j<EICs.size(); ++j)
        {
            const SpectrumIMS frameF = filterSpectrum(spec, 0, 0, mzStarts[j], mzEnds[j], mobilityStarts[j],
                                                      mobilityEnds[j]);
            EICs[j].times.push_back(fr.time);
            EICs[j].intensities.push_back(std::accumulate(frameF.intensities.begin(), frameF.intensities.end(), 0));
        }
    }
    
    Rcpp::List ret(EICs.size());
    for (size_t i=0; i<EICs.size(); ++i)
        ret[i] = Rcpp::DataFrame::create(Rcpp::Named("time") = EICs[i].times,
                                         Rcpp::Named("intensity") = EICs[i].intensities);
    return ret;
}
