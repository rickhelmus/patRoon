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

// based on https://stackoverflow.com/a/55779158
template <typename C> C median(std::vector<C> v)
{
    const auto n = v.size();
    const auto middleIt = v.begin() + n / 2;
    std::nth_element(v.begin(), middleIt, v.end());
    
    if (n % 2 == 0)
    {
        // get 2nd highest element from first half
        const auto middlePrevIt = std::max_element(v.begin(), middleIt);
        return (*middleIt + *middlePrevIt) / 2.0;
    }
    
    return *middleIt;
}

template <typename C, typename F> std::vector<size_t> getSortedInds(const C &cont, F func)
{
    // get sorted indices
    // inspired from https://stackoverflow.com/a/40183830
    std::vector<size_t> ret(cont.size());
    std::iota(ret.begin(), ret.end(), 0);
    std::sort(ret.begin(), ret.end(), [&](size_t i, size_t j) { return func(i, j); });
    return ret;
}

template <typename C> std::vector<size_t> getSortedInds(const C &cont)
{
    // get sorted indices
    // inspired from https://stackoverflow.com/a/40183830
    std::vector<size_t> ret(cont.size());
    std::iota(ret.begin(), ret.end(), 0);
    std::sort(ret.begin(), ret.end(), [&](size_t i, size_t j) { return cont[i] < cont[j]; });
    return ret;
}

template <typename IT> std::vector<size_t> getSortedInds(const IT &start, const IT &end)
{
    // get sorted indices
    // inspired from https://stackoverflow.com/a/40183830
    std::vector<size_t> ret(std::distance(start, end+1));
    std::iota(ret.begin(), ret.end(), 0);
    std::sort(ret.begin(), ret.end(), [&](size_t i, size_t j) { return *(start + i) < *(start + j); });
    return ret;
}

template <typename C> C movingAverage(const C &data, unsigned window)
{
    using T = typename std::remove_reference<typename C::value_type>::type;
    
    if (data.size() < window)
        window = data.size();
    if (window % 2 == 0)
        --window; // make window size odd
    if (window < 3)
        return data;
    
    const unsigned flank = (window - 1) / 2;
    C ret(data.size());
    
    // Process each element
    for (size_t i=0; i<data.size(); ++i)
    {
        // Determine the range of elements to average
        const size_t start = (i >= flank) ? i - flank : 0;
        const size_t end = ((i + flank) < data.size()) ? i + flank : data.size() - 1;
        
        T sum = T{};
        for (size_t j=start; j<=end; ++j)
            sum += data[j];
        
        ret[i] = sum / static_cast<T>(end - start + 1);
    }
    
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
    NumRange(void) : start(), end() { }
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
    
    if (nums.size() < 2)
        return ret; // no need to assign any clusters (for size==1 a one sized vector with value 0 will be returned)
    
    if (method == clusterMethod::BIN || method == clusterMethod::DISTANCE_MEAN || method == clusterMethod::DISTANCE_POINT)
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
        else // if (method == clusterMethod::DISTANCE_MEAN || method == clusterMethod::DISTANCE_POINT)
        {
            int curBinSize = 0;
            NumType binSum = 0;
            for (auto i : sortedNumInds)
            {
                if (curBinSize > 0)// && (nums[i] - (binSum / static_cast<NumType>(curBinSize))) > window)
                {
                    const auto diff = (method == clusterMethod::DISTANCE_MEAN)
                        ? nums[i] - (binSum / static_cast<NumType>(curBinSize))
                        : nums[i] - nums[i - 1];
                    if (diff > window)
                    {
                        ++curCluster;
                        curBinSize = 0;
                        binSum = 0;
                    }
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
                distm[k] = std::fabs(nums[i] - nums[j]);
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
