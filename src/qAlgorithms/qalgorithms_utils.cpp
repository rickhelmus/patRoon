#include "qalgorithms_utils.h"
#include <cstdint> // uint64_t
#define _USE_MATH_DEFINES
#include <math.h> // std::abs()
#include <cassert>
// #include "../external/CDFlib/cdflib.hpp"
#include "cephes.h"

#include <unordered_map>

namespace qAlgorithms
{
    bool F_test_regs(const double RSS_complex, const double RSS_simple,
                     const double params_complex, const double params_simple,
                     const double n, const double alpha)
    {
        const double fval = F_value(RSS_complex, RSS_simple, params_complex, params_simple, n);
        const double F_refval = F_stat(alpha, params_complex, params_simple, n);
        const bool f_ok = fval > F_refval; // reject H0, complex regression is significantly better than simple model
        return f_ok;
    }

    size_t hashm(int a, int b)
    {
        size_t a2 = a;
        a2 = a2 << __INT_WIDTH__;
        return a2 | b;
    }
    // this global hashmap is used to avoid recalculating the f value every time - better, thread-safe solution possible? @todo
    std::unordered_map<size_t, double> global_fhash_5perc;

    double F_stat(double alpha, size_t params_complex, size_t params_simple, size_t numPoints)
    {
        assert(params_complex > params_simple);
        assert(alpha > 0);
        assert(alpha < 1);

        int dfn_int = params_complex - params_simple;
        int dfd_int = numPoints - params_complex;
        size_t key = hashm(dfn_int, dfd_int);
        if (/*global_fhash_5perc.contains(key) && */ (alpha == 0.05)) // @todo generalise for all alpha
        {
            double Fhash = global_fhash_5perc[key];
            if (Fhash != 0)
                return Fhash;
        }

        double F = cephes::F_density(dfn_int, dfd_int, alpha);

#if false
        double F_old = 0; // return value
        {
            double q = alpha;
            double dfn = double(params_complex - params_simple); // numerator degrees of freedom
            double dfd = double(numPoints - params_complex);     // denominator degrees of freedom
            double p = 1 - alpha;                                // area of the covered distribution
            int which = 2;
            double bound = 0;
            int status = 1;
            cdff(&which, &p, &q, &F_old, &dfn, &dfd, &status, &bound); // library function, see https://people.math.sc.edu/Burkardt/cpp_src/cdflib/cdflib.html
            assert(status == 0);
        }
        assert(float(F) == float(F_old));
#endif

        global_fhash_5perc[key] = F;

        return F;
    }

    double F_value(const double RSS_complex, const double RSS_simple,
                   const double params_complex, const double params_simple,
                   const double n)
    {
        // Calculate F value of two models by their residual sum of squares (RSS) and number of regression
        // parameters (params). H0 is the model being compared against. n is the number of real points
        // both regressions are applied to. Refer to https://en.wikipedia.org/wiki/F-test#Regression_problems
        assert(params_complex > params_simple);
        assert(RSS_complex > 0);
        assert(RSS_simple > 0);
        assert(n > 1);
        double RSS_ratio = (RSS_simple - RSS_complex) / RSS_complex;
        double params_ratio = (n - params_complex) / (params_complex - params_simple);
        return RSS_ratio * params_ratio;
    }

    void linReg_intx(const float *yvals,
                     const size_t length,
                     double *slope, double *intercept)
    {
        assert(length > 2); // regression through two points is nonsensical;

        double s = 0, s_x = 0, s_y = 0, s_xx = 0, s_xy = 0; // accumulators for the different sums

        s = double(length);
        for (size_t x = 0; x < length; x++)
        {
            double y = yvals[x];

            s_x += double(x);
            s_y += y;
            s_xx += double(x * x);
            s_xy += double(x) * y;
        }

        double delta = 1 / (s * s_xx - s_x * s_x);
        *intercept = (s_xx * s_y - s_x * s_xy) * delta;
        *slope = (s * s_xy - s_x * s_y) * delta;
    }

    void weightedLinReg(const double *xvals,
                        const double *yvals,
                        const double *variance,
                        const size_t length,
                        double *slope, double *intercept)
    {
        /*
            The weighted linear regression as discussed in chapter 15.2 of numerical recepies (ISBN 0-521-43108-5).
            Variable names are taken from the text, not the included code example. The control for significance
            using the gamma function is skipped since there should always be a significant correlation. This could be
            replaced by a check against the mse when just using the mean.
        */
        assert(length > 2); // regression through two points is nonsensical;

        double s = 0, s_x = 0, s_y = 0, s_xx = 0, s_xy = 0; // accumulators for the different sums

        if (variance == nullptr)
        {
            s = double(length);
            for (size_t i = 0; i < length; i++)
            {
                double y = yvals[i];
                double x = xvals[i];

                s_x += x;
                s_y += y;
                s_xx += x * x;
                s_xy += x * y;
            }
        }
        else
        {
            for (size_t i = 0; i < length; i++)
            {
                double inverse = 1 / (variance[i] * variance[i]);

                double y = yvals[i];
                double x = xvals[i];

                s += inverse;
                s_x += x * inverse;
                s_y += y * inverse;
                s_xx += x * x * inverse;
                s_xy += x * y * inverse;
            }
        }

        // now, the two linear systems a * s + b * s_x = s_y and a * s_x + b * s_xx = s_xy can be solved for a and b
        // a is the inclination and b the intercept

        double delta = 1 / (s * s_xx - s_x * s_x);
        *intercept = (s_xx * s_y - s_x * s_xy) * delta;
        *slope = (s * s_xy - s_x * s_y) * delta;

        // calculate the uncertainties of a and b @todo
    }

    void coeffsQuadratic(const double x1, const double x2, const double x3,
                         const double y1, const double y2, const double y3,
                         double *b0, double *b1, double *b2)
    {
        // y1 = b2 x1^2 + b1 * x1 + b0 etc.
        const double div_x_23 = 1 / (x2 - x3);
        const double x22 = x2 * x2;
        const double x33 = x3 * x3;
        *b2 = ((y1 - y2) - (y2 - y3) * (x1 - x2) * div_x_23) / ((x1 * x1 - x22) - (x22 - x33) * (x1 - x2) * div_x_23);
        *b1 = (y2 - y3 + *b2 * (x33 - x22)) * div_x_23;
        *b0 = y3 - x33 * *b2 - x3 * *b1;
    }

    double quadraticAt(const double b0, const double b1, const double b2,
                       const double x)
    {
        return b0 + (b1 + x * b2) * x;
    }

    int solveQuadratic(const double a, const double b, const double c,
                       double *x1, double *x2)
    {
        // reference: https://www.av8n.com/physics/quadratic-formula.htm
        assert(a != 0), assert(b != 0), assert(c != 0);

        double sign = b <= 0 ? -1 : 1;
        double root = b * b - 4 * a * c;
        bool invalid = root <= 0;
        if (invalid)
        {
            *x1 = INFINITY;
            *x2 = INFINITY;
            return 1;
        }
        double x_big = (-b - sign * sqrt(root)) / (2 * a);
        double x_sml = c / (a * x_big);

        *x1 = min(x_sml, x_big);
        *x2 = max(x_sml, x_big);

        return 0;
    }

    double calcRSS(const float *predict,
                   const float *observed,
                   const Range_i *range)
    {
        double RSS = 0;
        for (size_t i = range->startIdx; i <= range->endIdx; i++)
        {
            double diff = predict[i] - observed[i];
            RSS += diff * diff;
        }
        return RSS;
    }

    double exp_approx_d(const double x)
    {
        return exp(x);
        // assert(x > 0); // @todo this is specified in the header but not respected throughout the code
        // assert(x < 26);
        // constexpr double LOG2E = 1.44269504088896340736;
        // constexpr double OFFSET = 1022.9329329329329;
        // constexpr uint64_t EXP_OFFSET = 1LL << 52;
        // union
        // {
        //     uint64_t i;
        //     double d;
        // } v = {(uint64_t)((x * LOG2E + OFFSET) * EXP_OFFSET)};
        // return v.d;
    }

    float erf_approx_f(const float x)
    {
        return erf(x);
        // float sign = x < 0 ? -1.0f : 1.0f; // get sign as the approximation is only valid for positive x
        // float t = std::abs(x);             // get the absolute value of x
        // constexpr float a1 = 0.278393f;    // empirically determined
        // constexpr float a2 = 0.230389f;    // empirically determined
        // constexpr float a3 = 0.000972f;    // empirically determined
        // constexpr float a4 = 0.078108f;    // empirically determined
        // t = 1.0f + t * (a1 + t * (a2 + t * (a3 + t * a4)));
        // t = t * t * t * t;               // t^4
        // return sign * (1.0f - 1.0f / t); // return the final approximation
    }

    double dawson5(double x)
    {
        double y, p, q;
        y = x * x;
        p = 1.0 + y * (0.1049934947 + y * (0.0424060604 + y * (0.0072644182 + y * (0.0005064034 + y * (0.0001789971)))));
        q = 1.0 + y * (0.7715471019 + y * (0.2909738639 + y * (0.0694555761 + y * (0.0140005442 + y * (0.0008327945 + 2 * 0.0001789971 * y)))));
        return x * (p / q);
    }

    double experfc(double x, double sign) // @todo this was updated in dev_G - understand, change and sync it
    {
        constexpr double a = 0.978795604954049; // empirically determined
        constexpr double b = 1.25731022692317;  // empirically determined
        double t = -x * x;
        return SQRTPI_2 * exp_approx_d(t) + sign * a * x * exp_approx_d(t * b); // @todo t is always negative, which goes against the definition given in exp_approx_d
    }

    double erfi_qalgo(const double x)
    {
        /* This function uses the Dawson Integral, i.e.
        erfi(x) = 2 * Dawson * exp(x^2) / sqrt(pi)

        The Dawson Integral is calculated by Stan Sykora's rational function approximations. URL : http://www.ebyte.it/library/codesnippets/DawsonIntegralApproximations.html (Dawson5())
        */
        // calculate the Dawson Integral:
        double y, p, q;
        y = x * x;
        p = 1.0 + y * (0.1049934947 + y * (0.0424060604 + y * (0.0072644182 + y * (0.0005064034 + y * (0.0001789971)))));
        q = 1.0 + y * (0.7715471019 + y * (0.2909738639 + y * (0.0694555761 + y * (0.0140005442 + y * (0.0008327945 + 2 * 0.0001789971 * y)))));
        double D = x * (p / q);

        return D * exp_approx_d(x * x);
    }

    // critical order space of two normally distributed populations
    double binningCritVal(const int n, const double stdDev)
    {
        // these values are determined empirically (see https://github.com/GeRe87/OS_critVal)
        const double OS_CRIT_A = 0.1443340625173891;
        const double OS_CRIT_B = 3.2412322699344687;
        return (OS_CRIT_A + (OS_CRIT_B / std::sqrt(std::log(n + 1)))) * stdDev;
    }

    size_t min(size_t a, size_t b)
    {
        return a < b ? a : b;
    }
    size_t max(size_t a, size_t b)
    {
        return a > b ? a : b;
    }
    int min(int a, int b)
    {
        return a < b ? a : b;
    }
    int max(int a, int b)
    {
        return a > b ? a : b;
    }
    double min(double a, double b)
    {
        return a < b ? a : b;
    }
    double max(double a, double b)
    {
        return a > b ? a : b;
    }

    float *minVal(float *const arrayStart, const size_t length)
    {
        assert(length > 0);
        float *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret < *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }
    const float *minVal(const float *const arrayStart, const size_t length)
    {
        assert(length > 0);
        const float *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret < *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }
    double *minVal(double *const arrayStart, const size_t length)
    {
        assert(length > 0);
        double *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret < *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }
    const double *minVal(const double *const arrayStart, const size_t length)
    {
        assert(length > 0);
        const double *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret < *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }

    float *maxVal(float *const arrayStart, const size_t length)
    {
        assert(length > 0);
        float *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret > *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }
    const float *maxVal(const float *const arrayStart, const size_t length)
    {
        assert(length > 0);
        const float *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret > *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }
    double *maxVal(double *const arrayStart, const size_t length)
    {
        assert(length > 0);
        double *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret > *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }
    const double *maxVal(const double *const arrayStart, const size_t length)
    {
        assert(length > 0);
        const double *ret = arrayStart;
        for (size_t i = 1; i < length; i++) // no need to check the first element
        {
            ret = *ret > *(arrayStart + i) ? ret : arrayStart + i;
        }
        return ret;
    }

    double meanOfCumulative(double *const cumArray, const size_t startIdx, const size_t endIdx)
    {
        assert(startIdx <= endIdx);
        double subtractor = startIdx == 0 ? 0 : cumArray[startIdx - 1]; // account for index 0
        double totalSum = cumArray[endIdx] - subtractor;
        return totalSum / (endIdx - startIdx + 1);
    }
    double meanOfCumulative(const double *const cumArray, const size_t startIdx, const size_t endIdx)
    {
        assert(startIdx <= endIdx);
        double subtractor = startIdx == 0 ? 0 : cumArray[startIdx - 1]; // account for index 0
        double totalSum = cumArray[endIdx] - subtractor;
        return totalSum / (endIdx - startIdx + 1);
    }

    unsigned int sumOfCumulative(const unsigned int *const cumArray, const Range_i *r)
    {
        // if the cumulative array does not exist (== null), assume that
        // the toal df is the length
        if (cumArray == nullptr)
            return rangeLen(r);

        // it is assumed that the range does not violate array bounds
        assert(r->startIdx <= r->endIdx);
        unsigned int subtractor = r->startIdx == 0 ? 0 : cumArray[r->startIdx - 1];
        unsigned int totalSum = cumArray[r->endIdx] - subtractor;
        return totalSum;
    }

    double sdev(double *const array, const size_t n)
    {
        assert(n > 2); // while standard deviation of two numbers is possible, it makes no sense
        double mean = 0;
        double sdev = 0;
        for (size_t i = 0; i < n; i++)
        {
            mean += array[i];
        }
        mean /= n;
        for (size_t i = 0; i < n; i++)
        {
            sdev += (array[i] - mean) * (array[i] - mean);
        }
        return sqrt(sdev / (n - 1));
    }

    double sdev(const double *const array, const size_t n)
    {
        assert(n > 2); // while standard deviation of two numbers is possible, it makes no sense
        double mean = 0;
        double sdev = 0;
        for (size_t i = 0; i < n; i++)
        {
            mean += array[i];
        }
        mean /= n;
        for (size_t i = 0; i < n; i++)
        {
            sdev += (array[i] - mean) * (array[i] - mean);
        }
        return sqrt(sdev / (n - 1));
    }

    double calcJaccardIdx(const float *const array1, const float *const array2, const size_t length)
    {
        // the jaccard index is defined as the intersection of two shapes divided by the union
        // The shapes are the two peak profiles here
        // We assume that both share the same x-axis, in which case every point can
        // be assumed to be at a distance of 1 to every other point.
        double intersect = 0;
        double union_val = 0;

        // calculate area as mean(array1[1], array[2]) * x etc. -> sum(array1) - array1[0] - array1[length - 1]
        for (size_t i = 0; i < length; i++)
        {
            double a1 = array1[i];
            double a2 = array2[i];
            intersect += min(a1, a2);
            union_val += max(a1, a2);
        }
        size_t end = length - 1;
        intersect -= (min(array1[0], array2[0]) + min(array1[end], array2[end])) / 2;
        union_val -= (max(array1[0], array2[0]) + max(array1[end], array2[end])) / 2;

        return intersect / union_val;
    }
}