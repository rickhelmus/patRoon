/*
 *  mshrmatchutils.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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
#include "mshrmatchutils.h"
#include "mselementmultmaputil.h"

using namespace std;

double ComputeHrMsMsExplainedIntens(
	const HrMassIntensMap& mapObserved,
	const set<ElementMultMap>& setEMM,
	int iPpmAccept, 
	int iPpmReject,
	bool bPositiveIonMode,
	map<double, list<ElementMultMap> >* pMapMassEMM)
{
	int iIonCharge=bPositiveIonMode?+1:-1;

	if(pMapMassEMM) pMapMassEMM->clear();

	HrMassIntensMap mapExplained;
	for(set<ElementMultMap>::const_iterator itEMM=setEMM.begin();itEMM!=setEMM.end();itEMM++)
	{
		double dMassCalc=GetNominalMass<double>(*itEMM)-iIonCharge*c_dMassElectron;
		HrMassIntensMap::const_iterator	itPeak=mapObserved.GetNearestPeak(dMassCalc);//GetMaxPeak(dMassCalc,iPPM);
		if(itPeak==mapObserved.end()) continue;

		/* old version
			if(fabs(CalcDiffPPM(dMassCalc,itPeak->first))>iPPM) continue;
			mapExplained[itPeak->first]=itPeak->second;
			if(pMapMassEMM) (*pMapMassEMM)[itPeak->first].push_back(*itEMM);
		*/

		double dFuzzyMassAccept=CalcFuzzyMassAccept(itPeak->first,dMassCalc,iPpmAccept,iPpmReject);

		if(dFuzzyMassAccept<=0) continue;

		if(pMapMassEMM) (*pMapMassEMM)[itPeak->first].push_back(*itEMM);

		//quick and dirty
		if(mapExplained[itPeak->first]<dFuzzyMassAccept*itPeak->second)
			mapExplained[itPeak->first]=dFuzzyMassAccept*itPeak->second;
	}

	// obsolete return bMassWeighted?
	//	mapExplained.GetTotalIntensMassProd()/mapObserved.GetTotalIntensMassProd():
	//  if Mass weighting is desired, it should already have been applied to mapObserved
	return mapExplained.GetTotalIntens()/mapObserved.GetTotalIntens();
}


double ComputeHrMsMsExplainedIntens(
	const HrMassIntensMap& mapObserved,
	const HrMassIntensMap& mapTheoretical,
	int iPpmAccept, int iPpmReject,
	HrMassIntensMap* pMapExplained)
{
	//obsolete: double dRelDiff=(double)iPPM/1000000.0;
	HrMassIntensMap mapOptional; //to be used, if pMapExplained==NULL
	HrMassIntensMap& mapExplained=pMapExplained?*pMapExplained:mapOptional;

	mapExplained.clear();

	for(HrMassIntensMap::const_iterator itObserved=mapObserved.begin();itObserved!=mapObserved.end();itObserved++)
	{
		HrMassIntensMap::const_iterator	itNearestPeak=mapTheoretical.GetNearestPeak(itObserved->first);

		if(itObserved==mapTheoretical.end()) continue;

		//if(fabs(itNearestPeak->first-itObserved->first)/itNearestPeak->first<=dRelDiff)
		//	mapExplained[itObserved->first]=itObserved->second;

		double dFuzzyMassAccept=CalcFuzzyMassAccept(itObserved->first,itNearestPeak->first,iPpmAccept,iPpmReject);

		mapExplained[itObserved->first]=itObserved->second*dFuzzyMassAccept;
	}

	// obsolete return bMassWeighted?
	//	mapExplained.GetTotalIntensMassProd()/mapObserved.GetTotalIntensMassProd():
	//  if Mass weighting is desired, it should already have been applied to mapObserved
	return mapExplained.GetTotalIntens()/mapObserved.GetTotalIntens();
}

double CalcFuzzyMassAccept(
	double dObserved,
	double dTheoretical,
	int iPpmAccept,
	int iPpmReject)
{
	double dMultAccept=(double)iPpmAccept/1000000.0;
	double dLeftAccept=dObserved*(1.0-dMultAccept);
	double dRightAccept=dObserved*(1.0+dMultAccept);

	if(dTheoretical<dLeftAccept)
	{
		double dLeftReject=dObserved*(1.0-(double)iPpmReject/1000000.0);

		if(dTheoretical<dLeftReject)
			return 0.0;
		else
			return (dTheoretical-dLeftReject)/(dLeftAccept-dLeftReject);
	}
	else	
	if(dTheoretical>dRightAccept)
	{
		double dRightReject=dObserved*(1.0+(double)iPpmReject/1000000.0);

		if(dTheoretical>dRightReject)
			return 0.0;
		else
			return (dRightReject-dTheoretical)/(dRightReject-dRightAccept);
	}

	return 1.0;
}
