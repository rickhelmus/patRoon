#include <string>
#include <algorithm>
#include <cctype>
#include <locale>
#include <cmath>

// UNDONE: for some reason we have to include Rcpp.h here to make clusterNums() work properly, what's going on here?!?!
#include <Rcpp.h>

#include "hclust-cpp/fastcluster.h"

#include "utils.h"

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

std::vector<int> clusterNums(const std::vector<double> &nums, clusterMethod method, double window)
{
    std::vector<int> ret(nums.size());
    
    if (method == clusterMethod::BIN || method == clusterMethod::DIFF)
    {
        const auto sortedNumInds = getSortedInds(nums);
        
        int curCluster = 0;
        if (method == clusterMethod::BIN)
        {
            int curBin = -1;
            for (auto i : sortedNumInds)
            {
                const int bin = static_cast<int>(nums[i] / window);
                if (curBin == -1 || bin != curBin)
                {
                    if (curBin != -1)
                        ++curCluster;
                    curBin = bin;
                }
                ret[i] = curCluster;
            }
        }
        else // clusterMethod::DIFF
        {
            int curBinSize = 0;
            double binSum = 0;
            for (auto i : sortedNumInds)
            {
                if (curBinSize == 0 || (nums[i] - (binSum / static_cast<double>(curBinSize))) > window)
                {
                    if (curBinSize > 0)
                        ++curCluster;
                    curBinSize = 1;
                    binSum = nums[i];
                }
                ret[i] = curCluster;
                binSum += nums[i];
                ++curBinSize;
            }
        }
    }
    else // clusterMethod::HCLUST
    {
        // get distance matrix, derived from hclust-cpp example
        const auto n = nums.size();
        auto distm = std::make_unique<double[]>((n * (n-1)) / 2);
        for (size_t i=0, k=0; i<n; ++i)
        {
            for (size_t j=i+1; j<n; ++j)
            {
                distm[k] = abs(nums[i] - nums[j]);
                ++k;
            }
        }
        auto merge = std::make_unique<int[]>(2 * (n-1));
        auto height = std::make_unique<double[]>(n - 1);
        hclust_fast(n, distm.get(), HCLUST_METHOD_COMPLETE, merge.get(), height.get());
        cutree_cdist(n, merge.get(), height.get(), window, ret.data());
    }

    return ret;
}

clusterMethod clustMethodFromStr(const std::string &str)
{
    if (str == "bin")
        return clusterMethod::BIN;
    else if (str == "diff")
        return clusterMethod::DIFF;
    else if (str == "hclust")
        return clusterMethod::HCLUST;
    
    Rcpp::stop("Unknown cluster method.");
    return clusterMethod::BIN; // avoid warning
}
