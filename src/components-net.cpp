/*
 * SPDX-FileCopyrightText: 2016-2026 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <vector>

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix getComponNetCoMatrix(Rcpp::List componGroupsList, Rcpp::List featGroupsList,
                                         Rcpp::CharacterVector gNames)
{
    Rcpp::NumericMatrix coCount(gNames.size(), gNames.size());
    Rcpp::NumericMatrix coPossible(gNames.size(), gNames.size());
    coCount.attr("dimnames") = Rcpp::List::create(gNames, gNames);
    
    const auto getGrpIndex = [&gNames](const char *grp)
    {
        const auto it = std::find(gNames.begin(), gNames.end(), grp);
        return std::distance(gNames.begin(), it);
    };
    
    for (size_t i=0; i<componGroupsList.size(); ++i)
    {
        const Rcpp::List cmps = componGroupsList[i];
        for (size_t j=0; j<cmps.size(); ++j)
        {
            const Rcpp::CharacterVector cmpGrps = cmps[j];
            for (size_t k=0; k<cmpGrps.size(); ++k)
            {
                const auto indk = getGrpIndex(cmpGrps[k]);
                for (size_t l=0; l<cmpGrps.size(); ++l)
                {
                    const auto indl = getGrpIndex(cmpGrps[l]);
                    ++coCount(indk, indl);
                    if (k != l)
                        ++coCount(indl, indk);
                }
                
            }
        }
        
        const Rcpp::CharacterVector featGrps = featGroupsList[i];
        for (size_t j=0; j<featGrps.size(); ++j)
        {
            const auto indj = getGrpIndex(featGrps[j]);
            for (size_t k=0; k<featGrps.size(); ++k)
            {
                const auto indk = getGrpIndex(featGrps[k]);
                ++coPossible(indj, indk);
                if (j != k)
                    ++coPossible(indk, indj);
            }
        }
    }
    
    // normalize
    for (int i=0; i<coCount.rows(); ++i)
    {
        for (int j=0; j<coCount.cols(); ++j)
            coCount(i, j) /= coPossible(i, j);
    }
    
    return coCount;
}
