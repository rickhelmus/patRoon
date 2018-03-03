/*
 *  msisotope.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_ISOTOPE_H
#define MS_ISOTOPE_H

#include "mschemelement.h"
#include "mslrmassintensmap.h"
#include "mselementmultmap.h"

const LrMassIntensMap& Init(LrMassIntensMap& MIM,const ElementMultMap& EMM);
const LrMassIntensMap& Init(
	LrMassIntensMap& MIM,
	const std::vector<int>& vecCode,
	const std::vector<int>& vecCount);

//void ReadIsotopeTableSIS(const std::string& strFileName);

template<class MASS> 
MASS GetNominalMass(
	const std::vector<int>& vecCode,
	const std::vector<int>& vecCount)
{
	MASS Sum=0;

	for(int iElem=0;iElem<vecCode.size();iElem++)
		Sum+=vecCount[iElem]*ChemElement(vecCode[iElem]).GetElementInfo().GetNominalMass(MASS());

	return Sum;
}

#endif
