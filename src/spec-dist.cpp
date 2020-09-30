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
    for (int i=0; i<specList.length(); ++i)
    {
        const Rcpp::DataFrame sp = Rcpp::as<Rcpp::DataFrame>(specList[i]);
        const Rcpp::NumericVector mzs = sp["mz"], ints = sp["intensity"];
        spectra.push_back(Spectrum{ Rcpp::as<std::vector<double>>(mzs), Rcpp::as<std::vector<double>>(ints) });
    }
    
    for (size_t i=0; i<len; ++i)
    {
        for (size_t j=i+1; j<=len; ++j)
        {
            // write to output matrix
            const double d = doCalcSpecSimilarity(spectra[i], spectra[j-1], meth, sh,
                                                  pmzs[j-1] - pmzs[i], mzw, intw, mzwin);
            ret(i, j-1) = d;
            ret(j-1, i) = d;
        }
    }
    
    return ret;
}
