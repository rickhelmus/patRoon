#ifndef PATROON_UTILS_HPP
#define PATROON_UTILS_HPP

// UNDONE: for some reason we have to include Rcpp.h here to make clusterNums() work properly, what's going on here?!?!
#include <Rcpp.h>

#include "hclust-cpp/fastcluster.h"

#include "utils.h"

#include <algorithm>
#include <exception>
#include <memory>
#include <numeric>
#include <vector>

template <typename C> std::vector<size_t> getSortedInds(const C &cont)
{
    // get sorted indices
    // inspired from https://stackoverflow.com/a/40183830
    std::vector<size_t> ret(cont.size());
    std::iota(ret.begin(), ret.end(), 0);
    std::sort(ret.begin(), ret.end(), [&](size_t i, size_t j) { return (cont[i] < cont[j]); });
    return ret;
}

// inspired from eg https://stackoverflow.com/questions/11828539/elegant-exceptionhandling-in-openmp
class ThreadExceptionHandler
{
    std::exception_ptr exPtr = nullptr;
    bool gotEx = false; // use separate bool flag which can be used with OpenMP atomic
    
    bool gotException(void) const
    {
        bool ret;
        #pragma omp atomic read
        ret = gotEx;
        return ret;
    }
    
public:
    template <typename F> void run(F func)
    {
        if (gotException())
            return;
        
        try
        {
            func();
        }
        catch (...)
        {
            #pragma omp critical
            {
                if (!exPtr)
                {
                    exPtr = std::current_exception();
                    gotEx = true;
                }
            }
        }
    }
    void reThrow(void) const
    {
        if (exPtr)
            std::rethrow_exception(exPtr);
    }
};

template <typename T> struct NumRange
{
    T start, end;
    NumRange(void) = default;
    NumRange(T s, T e) : start(s), end(e) { }
    void set(T s, T e) { start = s; end = e; }
    bool isSet(void) const { return start != 0 || end != 0; }
    bool overlap(const NumRange &o, T tol = 1E-8) const { return numberLTE(start, o.end, tol) && numberGTE(end, o.start, tol); }
    bool within(T v, T tol = 1E-8) const { return numberLTE(start, v, tol) && numberGTE(end, v, tol); }
};
template <typename T> NumRange<T> makeNumRange(T s, T e) { return NumRange<T>(s, e); }

template <typename NumType> std::vector<int> clusterNums(const std::vector<NumType> &nums, clusterMethod method,
                                                         NumType window)
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
            NumType binSum = 0;
            for (auto i : sortedNumInds)
            {
                if (curBinSize == 0 || (nums[i] - (binSum / static_cast<NumType>(curBinSize))) > window)
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

#endif
