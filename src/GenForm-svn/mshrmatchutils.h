/*
 *  mshrmatchutils.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_HR_MATCH_UTILS_H
#define MS_HR_MATCH_UTILS_H

#include <set>
#include "mshrmassintensmap.h"
#include "mselementmultmap.h"

//is called from ComputeHrMatchvalue as used for
//calculation of MS/MS MV for molecular formulas (mshrmatchbruttoformula.cpp)
//if sub-formulas explaining peak masses must be returned
double ComputeHrMsMsExplainedIntens(
	const HrMassIntensMap& mapObserved,
	const std::set<ElementMultMap>& setEMM,
	int iPpmAccept, int iPpmReject,
	bool bPositiveIonMode=true, //positive or negative charged ions?
	std::map<double, std::list<ElementMultMap> >* pMapMassEMM=NULL);

//is called from ComputeHrMatchvalue as used for
//calculation of MS/MS MV for molecular formulas (mshrmatchbruttoformula.cpp)
//if sub-formulas explaining peak masses must not be returned
double ComputeHrMsMsExplainedIntens(
	const HrMassIntensMap& mapObserved,
	const HrMassIntensMap& mapTheoretical,
	int iPpmAccept, int iPpmReject,
	HrMassIntensMap* pMapExplained=NULL);

/*
	helper function that returns value 1, if dTheoretical is within
	an interval of +/- iPpmAccept around dObserved, value 0 outside
	the +/- iPppmReject interval and and an interpolated value otherwise
*/
double CalcFuzzyMassAccept(
	double dObserved,
	double dTheoretical,
	int iPpmAccept,
	int iPpmReject);

#endif
