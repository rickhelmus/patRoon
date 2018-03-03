/*
 *  mshrmatchformula.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_HR_MATCH_FORMULA_H
#define MS_HR_MATCH_FORMULA_H

#include "mshrmassintensmap.h"
#include "mselementmultmap.h"

double ComputeHrMatchvalue(
	const HrMassIntensMap& mapObserved,
	const ElementMultMap& EMM,
	int iPpmAccept,//=10, 
	int iPpmReject,//=0, 
	int iExcessDBE,//=4,
	int iMinDiffMaxVal,
	int iMinDiffAtomCount,
	bool bAllowRadicalIons,//=true,
	bool bPositiveIonMode=true, //positive or negative charged ions?
	HrMassIntensMap* pMapExplained=NULL,
	std::map<double,std::list<ElementMultMap> >* pMapMassEMM=NULL);

#endif
