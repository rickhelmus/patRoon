/*
 * SPDX-FileCopyrightText: 2016 - 2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <algorithm>
#include <cctype>
#include <locale>
#include <memory>
#include <cmath>

#include "utils.h"
#include "utils.hpp"

// HACK
#include "hclust-cpp/fastcluster.cpp"


// ---
// Following three functions were taken from https://stackoverflow.com/a/217605

// trim from start (in place)
void ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
    {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
void rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch)
    {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
void trim(std::string &s)
{
    ltrim(s);
    rtrim(s);
}
// ---

bool hasWS(const std::string &s)
{
    return std::find_if(s.begin(), s.end(), [](unsigned char ch) { return std::isspace(ch); }) != s.end();
}

bool strStartsWith(const std::string &str, const std::string &pref)
{
    return(str.compare(0, pref.size(), pref) == 0);
}

bool compareTol(double x, double y, double tol)
{
    return std::fabs(x - y) <= tol;
}

bool numberLTE(double x, double y, double tol)
{
    return x < y || compareTol(x, y, tol);
}

bool numberGTE(double x, double y, double tol)
{
    return x > y || compareTol(x, y, tol);
}

void normalizeNums(std::vector<double> &v)
{
    double m = 0;
    for (double d : v)
        m = std::max(m, d);
    for (double &d : v)
        d /= m;
}

clusterMethod clustMethodFromStr(const std::string &str)
{
    if (str == "bin")
        return clusterMethod::BIN;
    else if (str == "distance")
        return clusterMethod::DISTANCE;
    else if (str == "hclust")
        return clusterMethod::HCLUST;
    
    Rcpp::stop("Unknown cluster method.");
    return clusterMethod::BIN; // avoid warning
}

int getOMPNumThreads()
{
#ifdef _OPENMP
    return omp_get_num_threads();
#else
    return 1;
#endif
}

// [[Rcpp::export]]
int getOMPMaxNumThreads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// [[Rcpp::export]]
void setOMPNumThreads(int n)
{
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
}

// [[Rcpp::export]]
double calcCenterOfMass(const Rcpp::NumericVector &x, const Rcpp::NumericVector &y)
{
    const auto size = x.size();
    
    if (size == 0)
        return 0.0;
    if (x.size() == 1)
        return x[0];
    
    double area = 0.0, wArea = 0.0;
    
    for (int i=1; i<x.size(); ++i)
    {
        const double dx = x[i] - x[i-1];
        const double avgY = (std::max(0.0, y[i]) + std::max(0.0, y[i-1])) / 2.0;
        const double avgX = (x[i] + x[i-1]) / 2.0;
        area += (dx * avgY);
        wArea += (dx * avgY * avgX);
    }
    
    if (area == 0.0)
        return 0.0; // avoid division by zero
    
    return wArea / area;
}

// [[Rcpp::export]]
std::vector<double> testMovingAverage(const std::vector<double> &data, int window)
{
    return movingAverage(data, window);
}