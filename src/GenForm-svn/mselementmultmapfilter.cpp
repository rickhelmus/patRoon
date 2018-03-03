/*
 *  mselementmultmapfilter.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mselementmultmaputil.h"
#include "mselementmultmapfilter.h"

#include <cmath>

using namespace std;

bool HeuerdingClercFilter(const ElementMultMap& EMM)
{
	int iMass=GetNominalMass<int>(EMM);

	int iC=EMM.GetMult(e_C);
	int iF=EMM.GetMult(e_F);
	int iCl=EMM.GetMult(e_Cl);
	int iI=EMM.GetMult(e_I);
	int iBr=EMM.GetMult(e_Br);

	if(iF+iCl+iI+iBr==0&&iC<0.07*(double)iMass-6.9) return false;
	if(iF+iCl+iI+iBr!=0&&iC<0.07*(double)iMass-14.0) return false;

	int iN=EMM.GetMult(e_N);
	int iO=EMM.GetMult(e_O);
	int iP=EMM.GetMult(e_P);
	int iS=EMM.GetMult(e_S);

	if(iN+iO+iP+iS>sqrt(0.2*(double)iMass)) return false;

	int iH=EMM.GetMult(e_H);

	if(iH+iF+iCl+iBr+iI<((double)(2*iC+iN+iP))/3.0-2.0/3.0) return false;

	return true;
}

