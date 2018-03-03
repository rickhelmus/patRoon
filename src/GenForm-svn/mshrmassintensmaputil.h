/*
 *  mshrmassintensmaputil.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_HR_MASS_INTENS_MAP_UTIL_H
#define MS_HR_MASS_INTENS_MAP_UTIL_H

#include "mshrmassintensmap.h"
#include <string>

void WeightSpectrum(
	HrMassIntensMap& MIM,
	const std::string strWeightMass,
	const std::string strWeightIntens,
	int iExponent);

void NormalizeLog(
	HrMassIntensMap& MIM,
	bool bMassWeight);

HrMassIntensMap& ShiftMass(
	const HrMassIntensMap& mapSource, 
	HrMassIntensMap& mapTarget, 
	double dMassShift);

#endif
