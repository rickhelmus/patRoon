#include "./qalgorithms_qpeaks.h"
#include "./qalgorithms_datatypes.h"
#include "./qalgorithms_utils.h"

#include <Rcpp.h>

#include <vector>

// [[Rcpp::export]]
Rcpp::List doFindPeaksQalgo(Rcpp::NumericVector x_values = Rcpp::NumericVector::create(),
                            Rcpp::NumericVector y_values = Rcpp::NumericVector::create(),
                            size_t maxScale = 20)
{
    // Return empty result if inputs are empty or incorrect
    if (y_values.length() == 0 || x_values.length() == 0)
        return Rcpp::List();
    if (y_values.length() != x_values.length())
        return Rcpp::List();

    size_t minScale = 2;
    maxScale = qAlgorithms::max(maxScale, minScale);
    size_t size = y_values.length();

    // since the default for R is a double, we have to do an inefficient copy here (sad!)
    std::vector<float> y_f(size);
    std::vector<float> x_f(size);
    for (size_t i = 0; i < size; i++)
    {
        y_f[i] = y_values[i];
        x_f[i] = x_values[i];
    }

    std::vector<qAlgorithms::PeakFit> peaks;
    int status = qAlgorithms::qpeaks_find(
        y_f.data(),
        x_f.data(),
        nullptr, // this is temporary
        y_values.size(),
        maxScale,
        &peaks);

    // move from array of structures to structure of arrays
    std::vector<double> ret, retmin, retmax, intensity, area, width;
    for (size_t i = 0; i < peaks.size(); i++)
    {
        ret.push_back(peaks[i].position);
        retmin.push_back(0);
        retmax.push_back(0);
        intensity.push_back(peaks[i].height);
        area.push_back(peaks[i].area);
        width.push_back(peaks[i].fwhm);
    }

    return Rcpp::List::create(
        Rcpp::Named("ret") = ret,
        Rcpp::Named("retmin") = retmin,
        Rcpp::Named("retmax") = retmax,
        Rcpp::Named("area") = area,
        Rcpp::Named("intensity") = intensity,
        Rcpp::Named("fwhm") = width);
}
