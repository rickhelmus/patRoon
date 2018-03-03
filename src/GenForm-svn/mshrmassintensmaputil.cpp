/*
 *  mshrmassintensmaputil.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mshrmassintensmaputil.h"
#include <cmath>

using namespace std;

void WeightSpectrum(HrMassIntensMap& MIM,const string strWeightMass,const string strWeightIntens,int iExponent)
{
	for(HrMassIntensMap::iterator it=MIM.begin();it!=MIM.end();it++)
	{
		double dMass=1.0,dIntens=1.0;

		if(strWeightMass=="lin"||strWeightMass.empty()) dMass=it->first;
		else
		if(strWeightMass=="sqrt") dMass=sqrt(it->first);
		else
		if(strWeightMass=="log") dMass=log10(1.0+it->first);

		if(strWeightIntens=="lin"||strWeightIntens.empty()) dIntens=it->second;
		else
		if(strWeightIntens=="sqrt") dIntens=sqrt(it->second);
		else
		if(strWeightIntens=="log") dIntens=log10(1.0+it->second*pow(10.0,(double)iExponent));

		it->second=dMass*dIntens;
	}	
}

void NormalizeLog(HrMassIntensMap& MIM,bool bMassWeight)
{
	MIM.Normalize(10000.0);
	for(HrMassIntensMap::iterator it=MIM.begin();it!=MIM.end();it++)
	{
		it->second*=(bMassWeight?it->first:100.0);

		if(it->second>10) 
			it->second=log(it->second);
		else 
			it->second=1.0;
	}
}

HrMassIntensMap& ShiftMass(const HrMassIntensMap& mapSource, HrMassIntensMap& mapTarget, double dMassShift)
{
	mapTarget.clear();
	HrMassIntensMap::iterator itTarget=mapTarget.begin();

	for(HrMassIntensMap::const_iterator itSource=mapSource.begin();itSource!=mapSource.end();itSource++)
		itTarget=mapTarget.insert(itTarget,make_pair(itSource->first+dMassShift,itSource->second));

	return mapTarget;
}

