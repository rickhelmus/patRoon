#include "opentims++/opentims_all.h"

#include <Rcpp.h>

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
