/*
 * SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#ifndef PATROON_UTILS_SPECTRA_H
#define PATROON_UTILS_SPECTRA_H

#include <vector>

struct Spectrum
{
    std::vector<int> IDs;
    std::vector<double> mzs, intensities;
};

double doCalcSpecSimilarity(Spectrum sp1, Spectrum sp2, const std::string &method,
                            const std::string &shift, double precDiff,
                            double mzWeight, double intWeight, double mzWindow);

#endif
