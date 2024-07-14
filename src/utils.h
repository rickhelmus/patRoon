#ifndef PATROON_UTILS_H
#define PATROON_UTILS_H

#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
inline int getOMPNumThreads(void) { return omp_get_num_threads(); }
#else
inline int getOMPNumThreads(void) { return 1; }
#endif

void ltrim(std::string &s);
void rtrim(std::string &s);
void trim(std::string &s);
bool hasWS(const std::string &s);
bool strStartsWith(const std::string &str, const std::string &pref);
bool compareTol(double x, double y, double tol = 1E-8);
bool numberLTE(double x, double y, double tol = 1E-8);
bool numberGTE(double x, double y, double tol = 1E-8);
void normalizeNums(std::vector<double> &v);

enum class clusterMethod { BIN, DIFF, HCLUST };
clusterMethod clustMethodFromStr(const std::string &str);

#include "utils.hpp"

#endif
