/*
 *  mselementintervalmaputil.cpp, part of GenForm by M. Meringer Copyright (C) 2015
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
*/

#ifndef MS_ELEMENT_INTERVAL_MAP_UTIL_H
#define MS_ELEMENT_INTERVAL_MAP_UTIL_H

#include "mselementintervalmap.h"

const ElementIntervalMap& Init(
		ElementIntervalMap& EIM,
		const std::string& strFormula);

const ElementIntervalMap& Init(
		ElementIntervalMap& EIM,
		int nLowerBound,
		const ElementMultMap& EMM);

#endif
