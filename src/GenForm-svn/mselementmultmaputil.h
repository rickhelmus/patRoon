/*
 *  mselementmultmaputil.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef MS_ELEMENT_MULT_MAP_UTIL_H
#define MS_ELEMENT_MULT_MAP_UTIL_H

#include "mschemelement.h"
#include "mselementmultmap.h"
#include "mselementinfo.h"

const ElementMultMap& Init(ElementMultMap& EMM,const std::string& strFormula);

bool LessChem(const ElementMultMap& EMM1,const ElementMultMap& EMM2);

bool HasHeteroAtom(const ElementMultMap& EMM);

std::string ChemName(const ElementMultMap& EMM);

std::string TexName(const ElementMultMap& EMM);

int GetDBE(const ElementMultMap& EMM,bool* pInteger=NULL);

double CalcDBE(const ElementMultMap& EMM);

bool IsGraphicalMultVal(const ElementMultMap& EMM);

template<class MASS>
MASS GetNominalMass(const ElementMultMap& EMM)
{
	MASS Sum=0;
	
	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
	{
		const ElementInfo& EI=ChemElement(it->first).GetElementInfo();
		Sum+=EI.GetNominalMass(MASS())*it->second;
	}	

	return Sum;
}

#endif
