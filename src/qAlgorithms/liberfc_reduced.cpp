#include "liberfc_reduced.h"
#include "auto_cheb_erfcx.h"

#define _USE_MATH_DEFINES
#include <math.h>

namespace liberfc
{
    const double spi2 = 0.88622692545275801364908374167057; // sqrt(pi)/2

    // internally used functions
    double im_w_of_x(double x);

    double erfi(double x)
    {
        // Compute erfi(x) = -i erf(ix),
        // the imaginary error function.

        return x * x > 720 ? (x > 0 ? INFINITY : -INFINITY) : exp(x * x) * im_w_of_x(x);
    }

    double dawson(double x)
    {

        // Compute dawson(x) = sqrt(pi)/2 * exp(-x^2) * erfi(x),
        // Dawson's integral for a real argument.

        return spi2 * im_w_of_x(x);
    }

    //! Simpler replacement for frexp from math.h, assuming that 0 < value < inf.
    //!
    //! Adapted from https://github.com/dioptre/newos/blob/master/lib/libm/arch/sh4/frexp.c.
    //! However, the mantissa must _not_ be broken into two variables to prevent errors
    //! on architectures like MIPS that do not revert the byte order of simple types.
    inline double frexp2(double value, int *eptr)
    {
        union
        {
            double v;
            struct
            {
#if ENDIAN_IS_LITTLE
                unsigned long long mantissa : 52;
                unsigned long long exponent : 11;
                unsigned long long sign : 1;
#else
                unsigned long long sign : 1;
                unsigned long long exponent : 11;
                unsigned long long mantissa : 52;
#endif
            } s;
        } u;

        u.v = value;
        *eptr = u.s.exponent - 1022;
        u.s.exponent = 1022;
        return u.v;
    }

    /* File im_w_of_x.c:
     *   Compute scaled Dawson integral im_w_of_x(x) = 2*dawson(x)/sqrt(pi),
     *   equivalent to the imaginary part of the Faddeeva function w(x) for real x.
     *
     * Website:
     *   http://apps.jcns.fz-juelich.de/libcerf
     *
     * Revision history:
     *   libcerf-2.5, April 2025:
     *      - For large |x|, continuous fractions replaced with asymptotic expansions.
     *      - Expanded the small |x| range where Taylor series is used.
     *      - Code for intermediate |x| entirely rewritten, using methods described in:
     *        Joachim Wuttke and Alexander Kleinsorge:
     *        "Code generation for piecewise Chebyshev approximation."
     *
     *   See also ../CHANGELOG
     *
     * Manual page:
     *   man 3 im_w_of_x
     */

    //! Returns polynomial approximation to f(x), using auto-tabulated expansion coefficients.
    //! Code taken from https://jugit.fz-juelich.de/mlz/ppapp.

    static double chebApproximant(double x)
    {
        // Application-specific constants:
        static const int loff = (ppapp_j0 + 1) * (1 << ppapp_M) + ppapp_l0; // precomputed offset

        // For given x, obtain mantissa xm and exponent je:
        int je;                           // will be set in next line
        const double xm = frexp2(x, &je); // sets xm and je

        // Integer arithmetics to obtain reduced coordinate t:
        const int ip = (int)((1 << (ppapp_M + 1)) * xm); // index in octave + 2^M
        const int lij = je * (1 << ppapp_M) + ip - loff; // index in lookup table
        const double t = (1 << (ppapp_M + 2)) * xm - (1 + 2 * ip);

        const double *const P = ppapp_Coeffs0 + lij * 8;
        const double *const Q = ppapp_Coeffs1 + lij * 4;
        return ((((((((((P[0]) * t + P[1]) * t + P[2]) * t + P[3]) * t + P[4]) * t + P[5]) * t + P[6]) * t + P[7]) * t + Q[0]) * t + Q[1]) * t + Q[2];
    }

    /******************************************************************************/
    /*  Library function im_w_of_z                                                */
    /******************************************************************************/

    double im_w_of_x(double x)
    {
        // Uses three different methods:
        // - asymptotic expansion for large |x|,
        // - Chebyshev polynomials for medium |x|,
        // - Taylor (Maclaurin) series for small |x|.

        const double ax = fabs(x); // very fast

        if (ax < .51)
        {
            // Use Taylor expansion (2/sqrt(pi)) * (x - 2/3 x^3  + 4/15 x^5  - 8/105 x^7 ...)

            const double x2 = x * x;

            if (ax < .083)
            {
                if (ax < 0.003)
                {
                    return ((((-0.085971746064420005629) * x2 // x^7
                              + 0.30090111122547001970) *
                                 x2 // x^5
                             - 0.75225277806367504925) *
                                x2 // x^3
                            + 1.1283791670955125739) *
                           x;
                }

                return (((((((+0.00053440090793734269229) * x2 // x^13
                             - 0.0034736059015927275001) *
                                x2 // x^11
                            + 0.019104832458760001251) *
                               x2 // x^9
                           - 0.085971746064420005629) *
                              x2 // x^7
                          + 0.30090111122547001970) *
                             x2 // x^5
                         - 0.75225277806367504925) *
                            x2 // x^3
                        + 1.1283791670955125739) *
                       x;
            }

            if (ax < .272)
            {
                return ((((((((((-8.82395720020380130481012927e-7) * x2 // x^19
                                + 8.38275934019361123956e-6) *
                                   x2 // x^17
                               - 7.1253454391645686483238e-5) *
                                  x2 // x^15
                              + 0.00053440090793734269229) *
                                 x2 // x^13
                             - 0.0034736059015927275001) *
                                x2 // x^11
                            + 0.019104832458760001251) *
                               x2 // x^9
                           - 0.085971746064420005629) *
                              x2 // x^7
                          + 0.30090111122547001970) *
                             x2 // x^5
                         - 0.75225277806367504925) *
                            x2 // x^3
                        + 1.1283791670955125739) *
                       x;
            }

            return (((((((((((((+5.8461000084165966602290712e-10) * x2 // x^25
                               - 7.30762501052074563638866034e-9) *
                                  x2 // x^23
                              + 8.40376876209885782941868884e-8) *
                                 x2 // x^21
                             - 8.82395720020380130481012927e-7) *
                                x2 // x^19
                            + 8.38275934019361123956e-6) *
                               x2 // x^17
                           - 7.1253454391645686483238e-5) *
                              x2 // x^15
                          + 0.00053440090793734269229) *
                             x2 // x^13
                         - 0.0034736059015927275001) *
                            x2 // x^11
                        + 0.019104832458760001251) *
                           x2 // x^9
                       - 0.085971746064420005629) *
                          x2 // x^7
                      + 0.30090111122547001970) *
                         x2 // x^5
                     - 0.75225277806367504925) *
                        x2 // x^3
                    + 1.1283791670955125739) *
                   x;
        }

        if (ax < 12.)
        {
            // Intermediate range: Use Chebyshev approximants.

            return copysign(chebApproximant(ax), x);
        }

        /* else */ {
            // Use asymptotic expansion up to N = 0, 3, 6, or 10
            //
            // With N=15 or 20 we could extend the range down to 7.73 or 6.72,
            // but we expect Chebyshev to be faster.
            //
            // Coefficient are a_0 = 1/sqrt(pi), a_N = (2N-1)!!/2^N/sqrt(pi).

            const double r = 1 / x;
            const double r2 = r * r;

            if (ax < 150)
            {
                if (ax < 23.2)
                {
                    return (((((((((((+3.607337150008375756e+05) * r2 +
                                     3.797197000008816394e+04) *
                                        r2 +
                                    4.467290588245667095e+03) *
                                       r2 +
                                   5.956387450994221808e+02) *
                                      r2 +
                                  9.163673001529572559e+01) *
                                     r2 +
                                 1.666122363914467641e+01) *
                                    r2 +
                                3.702494142032150659e+00) *
                                   r2 +
                               1.057855469152042982e+00) *
                                  r2 +
                              4.231421876608172372e-01) *
                                 r2 +
                             2.820947917738781396e-01) *
                                r2 +
                            5.641895835477562793e-01) *
                           r;
                }

                return (((((((+9.163673001529572559e+01) * r2 +
                             1.666122363914467641e+01) *
                                r2 +
                            3.702494142032150659e+00) *
                               r2 +
                           1.057855469152042982e+00) *
                              r2 +
                          4.231421876608172372e-01) *
                             r2 +
                         2.820947917738781396e-01) *
                            r2 +
                        5.641895835477562793e-01) *
                       r;
            }

            if (ax < 6.9e7)
            {
                return ((((+1.057855469152042982e+00) * r2 + 4.231421876608172372e-01) * r2 +
                         2.820947917738781396e-01) *
                            r2 +
                        5.641895835477562793e-01) *
                       r;
            }

            // 1-term expansion, important to avoid overflow
            return 5.641895835477562793e-01 / x;
        }

    } // im_w_of_z
}

/* Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013, 2025 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, 2025
 */