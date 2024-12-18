/*
 * SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#include <string>
#include <algorithm>
#include <cctype>
#include <locale>
#include <cmath>

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
