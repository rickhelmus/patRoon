#include "cephes.h"
#include <cassert>
#include <math.h>

namespace cephes
{
    // constants from const.h
    const double MACHEP = 1.11022302462515654042E-16; // 2**-53
    double MAXLOG = 7.09782712893383996732E2;         // log(MAXNUM)
    const double MINLOG = -7.451332191019412076235E2; // log(2**-1075)
#define MAXGAM 171.624376956302725

    double polevl(double x, const double coef[], int N)
    {
        double ans;
        int i;
        const double *p;

        p = coef;
        ans = *p++;
        i = N;

        do
            ans = ans * x + *p++;
        while (--i);

        return (ans);
    }

    /* p1evl()
     * Evaluate polynomial when coefficient of x^N  is 1.0.
     * Otherwise same as polevl.
     */

    double p1evl(double x, const double coef[], int N)
    {
        double ans;
        const double *p;
        int i;

        p = coef;
        ans = x + *p++;
        i = N - 1;

        do
            ans = ans * x + *p++;
        while (--i);

        return (ans);
    }

    double F_density(const int dfn, const int dfd, const double alpha)
    {
        return fdtri(dfn, dfd, alpha);
    }

    double fdtri(int df1, int df2, double p)
    {
        assert(df1 > 0);
        assert(df2 > 0);
        assert((0 < p) && (p <= 1));

        double df1_d = df1;
        double df2_d = df2;
        double ret;

        /* Compute probability for x = 0.5.  */
        double w = incbeta(0.5 * df2_d, 0.5 * df1_d, 0.5);
        /* If that is greater than y, then the solution w < .5.
           Otherwise, solve at 1-y to remove cancellation in (b - b*w).  */
        if (w > p || p < 0.001)
        {
            w = incbi(0.5 * df2_d, 0.5 * df1_d, p);
            ret = (df2_d - df2_d * w) / (df1_d * w);
        }
        else
        {
            w = incbi(0.5 * df1_d, 0.5 * df2_d, 1.0 - p);
            ret = df2_d * w / (df1_d * (1.0 - w));
        }
        return (ret);
    }

    /*
     * zlib License
     *
     * Regularized Incomplete Beta Function
     *
     * Copyright (c) 2016, 2017 Lewis Van Winkle
     * http://CodePlea.com
     *
     * This software is provided 'as-is', without any express or implied
     * warranty. In no event will the authors be held liable for any damages
     * arising from the use of this software.
     *
     * Permission is granted to anyone to use this software for any purpose,
     * including commercial applications, and to alter it and redistribute it
     * freely, subject to the following restrictions:
     *
     * 1. The origin of this software must not be misrepresented; you must not
     *    claim that you wrote the original software. If you use this software
     *    in a product, an acknowledgement in the product documentation would be
     *    appreciated but is not required.
     * 2. Altered source versions must be plainly marked as such, and must not be
     *    misrepresented as being the original software.
     * 3. This notice may not be removed or altered from any source distribution.
     */

#define STOP 1.0e-8
#define TINY 1.0e-30

    double incbeta(double a, double b, double x)
    {
        // modification due to https://github.com/codeplea/incbeta/issues/3
        assert(a > 0);
        assert(a > 0);

        if (x < 0.0 || x > 1.0)
            return INFINITY; // modification due to https://github.com/codeplea/incbeta/issues/2

        /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
        if (x > (a + 1.0) / (a + b + 2.0))
        {
            return (1.0 - incbeta(b, a, 1.0 - x)); /*Use the fact that beta is symmetrical.*/
        }

        /*Find the first part before the continued fraction.*/
        const double lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
        const double front = exp(log(x) * a + log(1.0 - x) * b - lbeta_ab) / a;

        /*Use Lentz's algorithm to evaluate the continued fraction.*/
        double f = 1.0, c = 1.0, d = 0.0;

        int i, m;
        for (i = 0; i <= 200; ++i)
        {
            m = i / 2;

            double numerator;
            if (i == 0)
            {
                numerator = 1.0; /*First numerator is 1.0.*/
            }
            else if (i % 2 == 0)
            {
                numerator = (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m)); /*Even term.*/
            }
            else
            {
                numerator = -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1)); /*Odd term.*/
            }

            /*Do an iteration of Lentz's algorithm.*/
            d = 1.0 + numerator * d;
            if (fabs(d) < TINY)
                d = TINY;
            d = 1.0 / d;

            c = 1.0 + numerator / c;
            if (fabs(c) < TINY)
                c = TINY;

            const double cd = c * d;
            f *= cd;

            /*Check for stop.*/
            if (fabs(1.0 - cd) < STOP)
            {
                return front * (f - 1.0);
            }
        }

        return INFINITY; /*Needed more loops, did not converge.*/
    }

    double student_t_cdf(double t, double v)
    {
        /*The cumulative distribution function (CDF) for Student's t distribution*/
        double x = (t + sqrt(t * t + v)) / (2.0 * sqrt(t * t + v));
        double prob = incbeta(v / 2.0, v / 2.0, x);
        return prob;
    }

    double incbi(double aa, double bb, double yy0)
    {
        double a, b, y0, d, y, x, x0, x1, lgm, yp, di, dithresh, yl, yh, xt;
        int i, dir;
        bool rflag, nflag;

        i = 0;
        if (yy0 <= 0)
            return (0.0);
        if (yy0 >= 1.0)
            return (1.0);
        x0 = 0.0;
        yl = 0.0;
        x1 = 1.0;
        yh = 1.0;
        nflag = false;

        if (aa <= 1.0 || bb <= 1.0)
        {
            dithresh = 1.0e-6;
            rflag = false;
            a = aa;
            b = bb;
            y0 = yy0;
            x = a / (a + b);
            y = incbeta(a, b, x);
            goto ihalve;
        }
        else
        {
            dithresh = 1.0e-4;
        }
        /* approximation to inverse function */

        yp = -ndtri(yy0);

        if (yy0 > 0.5)
        {
            rflag = true;
            a = bb;
            b = aa;
            y0 = 1.0 - yy0;
            yp = -yp;
        }
        else
        {
            rflag = false;
            a = aa;
            b = bb;
            y0 = yy0;
        }

        lgm = (yp * yp - 3.0) / 6.0;
        x = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
        d = yp * sqrt(x + lgm) / x - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0)) * (lgm + 5.0 / 6.0 - 2.0 / (3.0 * x));
        d = 2.0 * d;
        if (d < MINLOG)
        {
            x = 1.0;
            goto under;
        }
        x = a / (a + b * exp(d));
        y = incbeta(a, b, x);
        yp = (y - y0) / y0;
        if (fabs(yp) < 0.2)
            goto newt;

    /* Resort to interval halving if not close enough. */
    ihalve:

        dir = 0;
        di = 0.5;
        for (i = 0; i < 100; i++)
        {
            if (i != 0)
            {
                x = x0 + di * (x1 - x0);
                if (x == 1.0)
                    x = 1.0 - MACHEP;
                if (x == 0.0)
                {
                    di = 0.5;
                    x = x0 + di * (x1 - x0);
                    if (x == 0.0)
                        goto under;
                }
                y = incbeta(a, b, x);
                yp = (x1 - x0) / (x1 + x0);
                if (fabs(yp) < dithresh)
                    goto newt;
                yp = (y - y0) / y0;
                if (fabs(yp) < dithresh)
                    goto newt;
            }
            if (y < y0)
            {
                x0 = x;
                yl = y;
                if (dir < 0)
                {
                    dir = 0;
                    di = 0.5;
                }
                else if (dir > 3)
                    di = 1.0 - (1.0 - di) * (1.0 - di);
                else if (dir > 1)
                    di = 0.5 * di + 0.5;
                else
                    di = (y0 - y) / (yh - yl);
                dir += 1;
                if (x0 > 0.75)
                {
                    if (rflag)
                    {
                        rflag = false;
                        a = aa;
                        b = bb;
                        y0 = yy0;
                    }
                    else
                    {
                        rflag = true;
                        a = bb;
                        b = aa;
                        y0 = 1.0 - yy0;
                    }
                    x = 1.0 - x;
                    y = incbeta(a, b, x);
                    x0 = 0.0;
                    yl = 0.0;
                    x1 = 1.0;
                    yh = 1.0;
                    goto ihalve;
                }
            }
            else
            {
                x1 = x;
                if (rflag && x1 < MACHEP)
                {
                    x = 0.0;
                    goto done;
                }
                yh = y;
                if (dir > 0)
                {
                    dir = 0;
                    di = 0.5;
                }
                else if (dir < -3)
                    di = di * di;
                else if (dir < -1)
                    di = 0.5 * di;
                else
                    di = (y - y0) / (yh - yl);
                dir -= 1;
            }
        }
        // mtherr("incbi", PLOSS);
        if (x0 >= 1.0)
        {
            x = 1.0 - MACHEP;
            goto done;
        }
        if (x <= 0.0)
        {
        under:
            // mtherr("incbi", UNDERFLOW);
            x = 0.0;
            goto done;
        }

    newt:

        if (nflag)
            goto done;
        nflag = true;
        lgm = lgamma(a + b) - lgamma(a) - lgamma(b);

        for (i = 0; i < 8; i++)
        {
            /* Compute the function at this point. */
            if (i != 0)
                y = incbeta(a, b, x);
            if (y < yl)
            {
                x = x0;
                y = yl;
            }
            else if (y > yh)
            {
                x = x1;
                y = yh;
            }
            else if (y < y0)
            {
                x0 = x;
                yl = y;
            }
            else
            {
                x1 = x;
                yh = y;
            }
            if (x == 1.0 || x == 0.0)
                break;
            /* Compute the derivative of the function at this point. */
            d = (a - 1.0) * log(x) + (b - 1.0) * log(1.0 - x) + lgm;
            if (d < MINLOG)
                goto done;
            if (d > MAXLOG)
                break;
            d = exp(d);
            /* Compute the step to the next approximation of x. */
            d = (y - y0) / d;
            xt = x - d;
            if (xt <= x0)
            {
                y = (x - x0) / (x1 - x0);
                xt = x0 + 0.5 * y * (x - x0);
                if (xt <= 0.0)
                    break;
            }
            if (xt >= x1)
            {
                y = (x1 - x) / (x1 - x0);
                xt = x1 - 0.5 * y * (x1 - x);
                if (xt >= 1.0)
                    break;
            }
            x = xt;
            if (fabs(d / x) < 128.0 * MACHEP)
                goto done;
        }
        /* Did not converge.  */
        dithresh = 256.0 * MACHEP;
        goto ihalve;

    done:

        if (rflag)
        {
            if (x <= MACHEP)
                x = 1.0 - MACHEP;
            else
                x = 1.0 - x;
        }
        return (x);
    }

    const double LS2PI = 0.91893853320467274178; // log( sqrt( 2*pi ) )
    const double LOGPI = 1.14472988584940017414;

#define MAXLGM 2.556348e305

    // inverse of normal distribution

    /* approximation for 0 <= |y - 0.5| <= 3/8 */
    const double P0[5] = {
        -5.99633501014107895267E1,
        9.80010754185999661536E1,
        -5.66762857469070293439E1,
        1.39312609387279679503E1,
        -1.23916583867381258016E0,
    };
    const double Q0[8] = {
        /* 1.00000000000000000000E0,*/
        1.95448858338141759834E0,
        4.67627912898881538453E0,
        8.63602421390890590575E1,
        -2.25462687854119370527E2,
        2.00260212380060660359E2,
        -8.20372256168333339912E1,
        1.59056225126211695515E1,
        -1.18331621121330003142E0,
    };

    /* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
     * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
     */
    const double P1[9] = {
        4.05544892305962419923E0,
        3.15251094599893866154E1,
        5.71628192246421288162E1,
        4.40805073893200834700E1,
        1.46849561928858024014E1,
        2.18663306850790267539E0,
        -1.40256079171354495875E-1,
        -3.50424626827848203418E-2,
        -8.57456785154685413611E-4,
    };
    const double Q1[8] = {
        /*  1.00000000000000000000E0,*/
        1.57799883256466749731E1,
        4.53907635128879210584E1,
        4.13172038254672030440E1,
        1.50425385692907503408E1,
        2.50464946208309415979E0,
        -1.42182922854787788574E-1,
        -3.80806407691578277194E-2,
        -9.33259480895457427372E-4,
    };

    /* Approximation for interval z = sqrt(-2 log y ) between 8 and 64
     * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
     */
    const double P2[9] = {
        3.23774891776946035970E0,
        6.91522889068984211695E0,
        3.93881025292474443415E0,
        1.33303460815807542389E0,
        2.01485389549179081538E-1,
        1.23716634817820021358E-2,
        3.01581553508235416007E-4,
        2.65806974686737550832E-6,
        6.23974539184983293730E-9,
    };
    const double Q2[8] = {
        /*  1.00000000000000000000E0,*/
        6.02427039364742014255E0,
        3.67983563856160859403E0,
        1.37702099489081330271E0,
        2.16236993594496635890E-1,
        1.34204006088543189037E-2,
        3.28014464682127739104E-4,
        2.89247864745380683936E-6,
        6.79019408009981274425E-9,
    };
    const double s2pi = 2.50662827463100050242E0;

    double ndtri(double y0)
    {
        assert((0 <= y0) && (y0 <= 1));

        double x, y, z, y2, x0, x1;
        int code = 1;

        y = y0;
        if (y > (1.0 - 0.13533528323661269189)) /* 0.135... = exp(-2) */
        {
            y = 1.0 - y;
            code = 0;
        }

        if (y > 0.13533528323661269189)
        {
            y = y - 0.5;
            y2 = y * y;
            x = y + y * (y2 * polevl(y2, P0, 4) / p1evl(y2, Q0, 8));
            x = x * s2pi;
            return (x);
        }

        x = sqrt(-2.0 * log(y));
        x0 = x - log(x) / x;

        z = 1.0 / x;
        if (x < 8.0) /* y > exp(-32) = 1.2664165549e-14 */
            x1 = z * polevl(z, P1, 8) / p1evl(z, Q1, 8);
        else
            x1 = z * polevl(z, P2, 8) / p1evl(z, Q2, 8);
        x = x0 - x1;
        if (code != 0)
            x = -x;
        return (x);
    }
}