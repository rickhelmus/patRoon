#ifndef PATROON_UTILS_H
#define PATROON_UTILS_H

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

void ltrim(std::string &s);
void rtrim(std::string &s);
void trim(std::string &s);
bool hasWS(const std::string &s);
bool strStartsWith(const std::string &str, const std::string &pref);
bool compareTol(double x, double y, double tol);
bool numberLTE(double x, double y, double tol);
bool numberGTE(double x, double y, double tol);
void normalizeNums(std::vector<double> &v);

template <typename C> std::vector<size_t> getSortedInds(const C &cont)
{
    // get sorted indices
    // inspired from https://stackoverflow.com/a/40183830
    std::vector<size_t> ret(cont.size());
    std::iota(ret.begin(), ret.end(), 0);
    std::sort(ret.begin(), ret.end(), [&](size_t i, size_t j) { return (cont[i] < cont[j]); });
    return ret;
}

enum class clusterMethod { BIN, DIFF, HCLUST };
std::vector<int> clusterNums(const std::vector<double> &nums, clusterMethod method, double window);
clusterMethod clustMethodFromStr(const std::string &str);

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

#endif
