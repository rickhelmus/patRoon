/*
 *  mscomparelr.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include <cmath>
#include "mscomparelr.h"

using namespace std;

double CompareEuklidDistance(const LrMassIntensMap& mapObserved,const LrMassIntensMap& mapTheoretical)
{
	LrMassIntensMap MIM(mapTheoretical);
	MIM*=-1.0;
	MIM+=mapObserved;

	return sqrt(MIM.GetTotalIntensQuad());
}

double CompareNormalizedSumSquares(const LrMassIntensMap& mapObserved,const LrMassIntensMap& mapTheoretical)
{
	LrMassIntensMap mapSum(mapTheoretical),mapDiff(mapTheoretical);
	mapDiff*=-1.0;
	mapDiff+=mapObserved;
	mapSum+=mapObserved;

	return mapDiff.GetTotalIntensQuad()/mapSum.GetTotalIntensQuad();
}

double CompareNormalizedDotProduct(const LrMassIntensMap& mapObserved,const LrMassIntensMap& mapTheoretical)
{
	double dSumProduct=0.0;

	for(LrMassIntensMap::const_iterator it=mapObserved.begin();it!=mapObserved.end();it++)
		dSumProduct+=it->second*mapTheoretical.GetIntens(it->first);

	return 1.0-dSumProduct/sqrt(mapObserved.GetTotalIntensQuad()*mapTheoretical.GetTotalIntensQuad());
}


double CompareNormalizedDotProductSpecial(LrMassIntensMap mapObserved,LrMassIntensMap mapTheoretical)
//removes basepeak
{
	double dSumProduct=0.0;

//	mapObserved.Normalize();
	LrMassIntensMap::iterator itObs=mapObserved.GetBasePeak();
	LrMassIntensMap::iterator itTheo=mapTheoretical.GetBasePeak();

//	mapTheoretical*=(itObs->second/itTheo->second);

	mapObserved.erase(itObs);
	mapTheoretical.erase(itTheo);

	for(LrMassIntensMap::const_iterator it=mapObserved.begin();it!=mapObserved.end();it++)
		dSumProduct+=it->second*mapTheoretical.GetIntens(it->first);

	return 1.0-dSumProduct/sqrt(mapObserved.GetTotalIntensQuad()*mapTheoretical.GetTotalIntensQuad());
}


double CompareShimadzuSimilarityIndex(const LrMassIntensMap& mapObserved,const LrMassIntensMap& mapTheoretical)
{
	double dSum=0,dDif=0;

	LrMassIntensMap mapDif=mapObserved+(mapTheoretical*(-1.0));
	LrMassIntensMap mapSum=mapObserved+mapTheoretical;
	LrMassIntensMap::const_iterator it;

	for(it=mapDif.begin();it!=mapDif.end();it++) dDif+=abs(it->second);
	for(it=mapSum.begin();it!=mapSum.end();it++) dSum+=it->second;

	return dDif/dSum;
}


double CompareShimadzuSimilarityIndexSpecial(LrMassIntensMap mapObserved,LrMassIntensMap mapTheoretical)
//normalizes to basepeak and removes basepeak
{
	double dSum=0,dDif=0;

	mapObserved.Normalize();
	LrMassIntensMap::iterator itObs=mapObserved.GetBasePeak();
	LrMassIntensMap::iterator itTheo=mapTheoretical.find(itObs->first);
	
	if(itTheo==mapObserved.end()) return 0.0;

	mapTheoretical*=(itObs->second/itTheo->second);

	mapObserved.erase(itObs);
	mapTheoretical.erase(itTheo);

	//mapTheoretical.erase(mapTheoretical.find(mapObserved.rbegin()->first+1),mapTheoretical.end());

	LrMassIntensMap mapDif=mapObserved+(mapTheoretical*(-1.0));
	LrMassIntensMap mapSum=mapObserved+mapTheoretical;
	LrMassIntensMap::const_iterator it;

	for(it=mapDif.begin();it!=mapDif.end();it++) dDif+=fabs(it->second);
	for(it=mapSum.begin();it!=mapSum.end();it++) dSum+=it->second;
	//for(it=mapDif.begin();it!=mapDif.end();it++) dDif+=it->second*it->second;
	//for(it=mapSum.begin();it!=mapSum.end();it++) dSum+=it->second*it->second;

	return dDif/dSum;
}
