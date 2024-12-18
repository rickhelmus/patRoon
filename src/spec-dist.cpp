/*
 * SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <Rcpp.h>

#include "utils-spectra.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix specDistMatrix(Rcpp::List specList, Rcpp::CharacterVector method,
                                   Rcpp::CharacterVector shift, Rcpp::NumericVector precMZs,
                                   Rcpp::NumericVector mzWeight, Rcpp::NumericVector intWeight,
                                   Rcpp::NumericVector mzWindow)
{
    Rcpp::NumericMatrix ret(specList.length(), specList.length());
    const size_t len = specList.length();
    const std::string meth = Rcpp::as<std::string>(method), sh = Rcpp::as<std::string>(shift);
    const std::vector<double> pmzs = Rcpp::as<std::vector<double>>(precMZs);
    const double mzw = Rcpp::as<double>(mzWeight), intw = Rcpp::as<double>(intWeight), mzwin = Rcpp::as<double>(mzWindow);
    
    std::vector<Spectrum> spectra;
    for (size_t i=0; i<len; ++i)
    {
        const Rcpp::DataFrame sp = Rcpp::as<Rcpp::DataFrame>(specList[i]);
        const Rcpp::NumericVector ids = sp["ID"], mzs = sp["mz"], ints = sp["intensity"];
        spectra.push_back(Spectrum{ Rcpp::as<std::vector<int>>(ids), Rcpp::as<std::vector<double>>(mzs),
                                    Rcpp::as<std::vector<double>>(ints) });
    }
    
    for (size_t i=0; i<len; ++i)
    {
        for (size_t j=i+1; j<=len; ++j)
        {
            // write to output matrix
            const double d = doCalcSpecSimilarity(spectra[i], spectra[j-1], meth, sh, pmzs[j-1] - pmzs[i], mzw,
                                                  intw, mzwin);
            ret(i, j-1) = d;
            ret(j-1, i) = d;
        }
    }
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix specDistRect(Rcpp::List specList1, Rcpp::List specList2, Rcpp::CharacterVector method,
                                 Rcpp::CharacterVector shift, Rcpp::NumericVector precMZs1,
                                 Rcpp::NumericVector precMZs2, Rcpp::NumericVector mzWeight,
                                 Rcpp::NumericVector intWeight, Rcpp::NumericVector mzWindow)
{
    Rcpp::NumericMatrix ret(specList1.length(), specList2.length());
    const std::string meth = Rcpp::as<std::string>(method), sh = Rcpp::as<std::string>(shift);
    const std::vector<double> pmzs1 = Rcpp::as<std::vector<double>>(precMZs1);
    const std::vector<double> pmzs2 = Rcpp::as<std::vector<double>>(precMZs2);
    const double mzw = Rcpp::as<double>(mzWeight), intw = Rcpp::as<double>(intWeight), mzwin = Rcpp::as<double>(mzWindow);
    
    const size_t len1 = specList1.length(), len2 = specList2.length();
    std::vector<Spectrum> spectra1, spectra2;
    
    for (size_t i=0; i<len1; ++i)
    {
        const Rcpp::DataFrame sp = Rcpp::as<Rcpp::DataFrame>(specList1[i]);
        const Rcpp::NumericVector ids = sp["ID"], mzs = sp["mz"], ints = sp["intensity"];
        spectra1.push_back(Spectrum{ Rcpp::as<std::vector<int>>(ids), Rcpp::as<std::vector<double>>(mzs),
                                     Rcpp::as<std::vector<double>>(ints) });
    }
    for (size_t i=0; i<len2; ++i)
    {
        const Rcpp::DataFrame sp = Rcpp::as<Rcpp::DataFrame>(specList2[i]);
        const Rcpp::NumericVector ids = sp["ID"], mzs = sp["mz"], ints = sp["intensity"];
        spectra2.push_back(Spectrum{ Rcpp::as<std::vector<int>>(ids), Rcpp::as<std::vector<double>>(mzs),
                                     Rcpp::as<std::vector<double>>(ints) });
    }

    for (size_t i=0; i<len1; ++i)
    {
        for (size_t j=0; j<len2; ++j)
        {
            // write to output matrix
            const double d = doCalcSpecSimilarity(spectra1[i], spectra2[j], meth, sh, pmzs2[j] - pmzs1[i], mzw,
                                                  intw, mzwin);
            ret(i, j) = d;
        }
    }
    
    return ret;
}
