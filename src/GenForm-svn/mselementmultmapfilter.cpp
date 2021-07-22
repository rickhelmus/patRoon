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

bool KindFiehnElementRatios(const ElementMultMap& EMM, bool bExtended)
{
	// ranges for element ratios corresponding to Table 2 of
	// Kind, Tobias, and Oliver Fiehn. Seven Golden Rules for heuristic filtering of molecular formulas
	// obtained by accurate mass spectrometry. BMC bioinformatics 8.1 (2007): 105.
	// https://link.springer.com/article/10.1186/1471-2105-8-105/tables/2

	double dC=EMM.GetMult(e_C);
	double dRatio=EMM.GetMult(e_H)/dC;

	if(!bExtended)
	{
		if(dRatio<0.2 || dRatio>3.1) return false;
	}
	else
	{
		if(dRatio<0.1 || dRatio>6.0) return false;
	}

	const int nEl=8;
	const int pEl[nEl]={e_F,e_Cl,e_Br,e_N,e_O,e_P,e_S,e_Si};
	const double pMaxCommon[nEl]={1.5,0.8,0.8,1.3,1.2,0.3,0.8,0.5};
	const double pMaxExtended[nEl]={6.0,2.0,2.0,4.0,3.0,2.0,3.0,1.0};

	for(int iEl=0;iEl<nEl;iEl++)
	{
		dRatio=EMM.GetMult(pEl[iEl])/dC;

		if(!bExtended)
		{
			if(dRatio>pMaxCommon[iEl]) return false;
		}
		else
		{
			if(dRatio>pMaxExtended[iEl]) return false;
		}
	}

	return true;
}


