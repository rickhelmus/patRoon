#include <Rcpp.h>

#include "utils-spectra.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix specDistMatrix(Rcpp::List specList, Rcpp::CharacterVector method,
                                   Rcpp::CharacterVector shift, Rcpp::NumericVector mzWeight,
                                   Rcpp::NumericVector intWeight, Rcpp::NumericVector mzWindow)
{
    Rcpp::NumericMatrix ret(specList.length(), specList.length());
    const size_t len = specList.length();
    const std::string meth = Rcpp::as<std::string>(method), sh = Rcpp::as<std::string>(shift);
    const double mzw = Rcpp::as<double>(mzWeight), intw = Rcpp::as<double>(intWeight), mzwin = Rcpp::as<double>(mzWindow);
    const double precDiff = 0.0; // UNDONE
    
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
            // // rows we will operate on
            // Rcpp::NumericMatrix::Row row1 = mat.row(i);
            // Rcpp::NumericMatrix::Row row2 = mat.row(j);
            // 
            // // compute the average using std::tranform from the STL
            // std::vector<double> avg(row1.size());
            // std::transform(row1.begin(), row1.end(), // input range 1
            //                row2.begin(),             // input range 2
            //                avg.begin(),              // output range 
            //                average);                 // function to apply
            // 
            // // calculate divergences
            // double d1 = kl_divergence(row1.begin(), row1.end(), avg.begin());
            // double d2 = kl_divergence(row2.begin(), row2.end(), avg.begin());
            
            // write to output matrix
            const double d = doCalcSpecSimularity(spectra[i], spectra[j-1], meth, sh,
                                                  precDiff, mzw, intw, mzwin);
            ret(i, j-1) = d;
            ret(j-1, i) = d;
        }
    }
    
    return ret;
}
