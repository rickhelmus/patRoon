/*
 * SPDX-FileCopyrightText: 2016-2024 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#ifndef PATROON_UTILS_H
#define PATROON_UTILS_H

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

#endif
