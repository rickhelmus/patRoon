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