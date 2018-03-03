/*
 *  msconvert.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "msconvert.h"

using namespace std;

HrMassIntensMap& Init(HrMassIntensMap& HR,const LrMassIntensMap& LR)
{
	HR.clear();
	HrMassIntensMap::iterator itHR=HR.begin();

	for(LrMassIntensMap::const_iterator itLR=LR.begin();itLR!=LR.end();itLR++)
		itHR=HR.insert(itHR,pair<double,double>(itLR->first,itLR->second));

	return HR;
}

LrMassIntensMap& Init(LrMassIntensMap& LR,const HrMassIntensMap& HR, double dShift)
{
	LR.clear();
	LrMassIntensMap::iterator itLR=LR.begin();

	for(map<double,double>::const_iterator itHR=HR.begin();itHR!=HR.end();itHR++)
	{
		unsigned int iPeakCount=LR.size(), iMass=(unsigned int)(itHR->first+dShift);
		itLR=LR.insert(itLR,make_pair(iMass,itHR->second));
		if(iPeakCount==LR.size()) //no new integer peak mass inserted; add intensity
			itLR->second+=itHR->second;
	}
/*
 	equivalent, but probably slower code:
	for(HrMassIntensMap::const_iterator itHR=HR.begin();itHR!=HR.end();itHR++)
		LR[(int)(itHR->first+dShift)]+=itHR->second;
*/
	return LR;
}

