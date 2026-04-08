// this is a wrapper around the cephes mathematical library found here:
// https://www.moshier.net/doubldoc.html

// function definitions have been changed to being c++ compatible
// the parts of the library oriented towards outdated systems have been removed.
// error handling for function input has been replaced by assertions.

#ifndef CEPHES
#define CEPHES

namespace cephes
{
    /**
     * @brief Evaluate polynomial.
     *
     * @param x The variable for the polynomial.
     * @param coef Array of coefficients in reverse order.
     * @param N Degree of the polynomial.
     * @return The evaluated polynomial value y.
     *
     * @details
     * This function evaluates a polynomial of degree N:
     *
     * y = C0 + C1*x + C2*x^2 + ... + CN*x^N
     *
     * Coefficients are stored in reverse order:
     *
     * coef[0] = CN, ..., coef[N] = C0.
     *
     * The function p1evl() assumes that coef[N] = 1.0 and is omitted from
     * the array. Its calling arguments are otherwise the same as for polevl().
     *
     * @section speed Speed
     *
     * In the interest of speed, there are no checks for out-of-bounds
     * arithmetic. This routine is used by most of the functions in the library.
     * Depending on available equipment features, the user may wish to rewrite
     * the program in microcode or assembly language.
     *
     * @note
     * Cephes Math Library Release 2.1: December, 1988
     * Copyright 1984, 1987, 1988 by Stephen L. Moshier
     */
    double polevl(double x, const double coef[], int N);

    /* p1evl()	*/
    /*
     * Evaluate polynomial when coefficient of x^N is 1.0.
     * Otherwise same as polevl.
     */
    double p1evl(double x, const double coef[], int N);

    /**
     * @brief Inverse of the complemented F distribution.
     *
     * @param df1 Degrees of freedom for the numerator, must be a positive integer
     * @param df2 Degrees of freedom for the denominator must be a positive integer
     * @param p The probability value.
     * @return The F density argument x such that the integral from x to infinity
     *         of the F density is equal to the given probability p.
     *
     * @details
     * This function finds the F density argument x by using the inverse beta
     * integral function and the following relations:
     *
     * - z = incbi(df2/2, df1/2, p)
     * - x = df2 * (1 - z) / (df1 * z)
     *
     * Note: The following relations hold for the inverse of the uncomplemented
     * F distribution:
     *
     * - z = incbi(df1/2, df2/2, p)
     * - x = df2 * z / (df1 * (1 - z)).
     *
     * @section Accuracy
     * Tested at random points (a,b,p).
     *
     * | a,b       | # trials | Peak        | RMS          |
     * |-----------|----------|-------------|--------------|
     * | IEEE 1,100        | 100000   | 8.3e-15 | 4.7e-16  |
     * | IEEE 1,10000      | 100000   | 2.1e-11 | 1.4e-13  |
     * | For p between 10^-6 and 10^-3:                    |
     * | IEEE 1,100        | 50000    | 1.3e-12 | 8.4e-15  |
     * | IEEE 1,10000      | 50000    | 3.0e-12 | 4.8e-14  |
     *
     * See also fdtrc.c.
     *
     * @note
     * Cephes Math Library Release 2.8: June, 2000
     * Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
     */
    double F_density(const int df1, const int df2, const double alpha);
    double fdtri(int df1, int df2, double p);

    /**
     * @brief Incomplete beta integral.
     *
     * @param a First parameter of the integral.
     * @param b Second parameter of the integral.
     * @param x Argument at which to evaluate the integral.
     * @return The incomplete beta integral of the arguments, evaluated
     *         from zero to x.
     *
     * @details
     * This function returns the incomplete beta integral defined as:
     *
     *                  x
     *     -            -
     *    | (a+b)      | |  a-1     b-1
     *  -----------    |   t   (1-t)   dt.
     *   -     -     | |
     *  | (a) | (b)   -
     *                0
     *
     * The domain of definition is 0 <= x <= 1. In this
     * implementation, a and b are restricted to positive values.
     * The integral from x to 1 may be obtained by the symmetry relation:
     *
     * 1 - incbet(a, b, x)  =  incbet(b, a, 1 - x ).
     *
     * The integral is evaluated using either a continued fraction expansion
     * or, when b * x is small, by a power series.
     *
     * @section Accuracy
     * @todo
     *
     * @note
     * zlib License
     *
     * Regularized Incomplete Beta Function
     *
     * Copyright (c) 2016, 2017 Lewis Van Winkle
     * http://CodePlea.com
     */
    double incbeta(double a, double b, double x);

    double student_t_cdf(double t, double v);

    /**
     * @brief Inverse of incomplete beta integral.
     *
     * @param a First parameter of the integral.
     * @param b Second parameter of the integral.
     * @param y The probability value for which to find x.
     * @return The value of x such that incbet(a, b, x) = y.
     *
     * @details
     * Given y, this function finds x such that:
     *
     * incbeta(a, b, x) = y.
     *
     * The routine performs interval halving or Newton iterations to find the
     * root of incbeta(a, b, x) - y = 0.
     *
     * @section Accuracy
     *
     * | arithmetic | x    | a,b   | domain   | # trials | Peak        | RMS         |
     * |------------|------|-------|----------|----------|-------------|-------------|
     * | IEEE       | 0,1  | .5,10000 | 50000    | 5.8e-12   | 1.3e-13     |
     * | IEEE       | 0,1  | .25,100   | 100000   | 1.8e-13   | 3.9e-15     |
     * | IEEE       | 0,1  | 0,5      | 50000    | 1.1e-12   | 5.5e-15     |
     * | VAX        | 0,1  | .5,100   | 25000    | 3.5e-14   | 1.1e-15     |
     *
     * With a and b constrained to half-integer or integer values:
     * | IEEE       | 0,1  | .5,10000 | 50000    | 5.8e-12   | 1.1e-13     |
     * | IEEE       | 0,1  | .5,100  | 100000   | 1.7e-14   | 7.9e-16     |
     *
     * With a = .5, b constrained to half-integer or integer values:
     * | IEEE       | 0,1  | .5,10000 | 10000    | 8.3e-11   | 1.0e-11     |
     *
     * @note
     * Cephes Math Library Release 2.8: June, 2000
     * Copyright 1984, 1996, 2000 by Stephen L. Moshier
     */
    double incbi(double aa, double bb, double yy0);

    /**
     * @brief Inverse of the Normal distribution function.
     *
     * @param y The area under the Gaussian probability density function
     *          (integrated from minus infinity to x).
     * @return The corresponding value of x for which the area is equal to y.
     *
     * @details
     * This function returns the argument, x, for which the area under the
     * Gaussian probability density function (integrated from minus infinity
     * to x) is equal to y.
     *
     * For small arguments (0 < y < exp(-2)), the program computes:
     * z = sqrt(-2.0 * log(y)); then the approximation is:
     *
     * x = z - log(z)/z - (1/z) P(1/z) / Q(1/z).
     *
     * There are two rational functions P/Q, one for 0 < y < exp(-32)
     * and the other for y up to exp(-2). For larger arguments, w = y - 0.5,
     * and x/sqrt(2pi) = w + w^3 R(w^2)/S(w^2).
     *
     * @section Accuracy
     *
     * | arithmetic | domain               | # trials | Peak         | RMS         |
     * |------------|----------------------|----------|--------------|-------------|
     * | DEC        | 0.125, 1             | 5500     | 9.5e-17      | 2.1e-17     |
     * | DEC        | 6e-39, 0.135         | 3500     | 5.7e-17      | 1.3e-17     |
     * | IEEE       | 0.125, 1             | 20000    | 7.2e-16      | 1.3e-16     |
     * | IEEE       | 3e-308, 0.135        | 50000    | 4.6e-16      | 9.8e-17     |
     *
     * @section Error Messages
     *
     * | Message          | Condition              | Value Returned |
     * |-------------------|-----------------------|-----------------|
     * | ndtri domain      | x <= 0                | -MAXNUM         |
     * | ndtri domain      | x >= 1                | MAXNUM          |
     */
    double ndtri(double y0);
}
#endif