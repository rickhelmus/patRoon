#ifndef PATROON_UTILS_SPECTRA_H
#define PATROON_UTILS_SPECTRA_H

#include <vector>

struct Spectrum
{
    std::vector<double> mzs, intensities;
};

double doCalcSpecSimilarity(Spectrum sp1, Spectrum sp2, const std::string &method,
                            const std::string &shift, double precDiff,
                            double mzWeight, double intWeight, double mzWindow);

#endif
