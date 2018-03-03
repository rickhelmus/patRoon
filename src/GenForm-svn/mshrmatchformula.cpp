/*
 *  mshrmatchformula.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mshrmatchformula.h"
#include "msgenformbyspecalgo.h"
#include "msconvert.h"
#include "mselementmultmaputil.h"
#include "mselementintervalmaputil.h"
#include "mshrmatchutils.h"

using namespace std;

double ComputeHrMatchvalue(
	const HrMassIntensMap& mapObserved,
	const ElementMultMap& EMM,
	int iPpmAccept, 
	int iPpmReject, 
	int iExcessDBE,
	int iMinDiffMaxVal,
	int iMinDiffAtomCount,
	bool bAllowRadicalIons,
	bool bPositiveIonMode, 
	HrMassIntensMap* pMapExplained,
	std::map<double,std::list<ElementMultMap> >* pMapMassEMM)
{
	HrMassIntensMap mapTheoretical;
	LrMassIntensMap MIM;
	Init(MIM,mapObserved);
	// There is a problem for peaks that belong to  molecular formulae when accurate masses differ >> 0.5 form integer masses, 
	// shift argument like GetNominalMass<int>(EMM)-GetNominalMass<double>(EMM)+0.5) might be more appropriate, but won't work
	// in general, because small peaks are then shifted too much
	set<ElementMultMap> setEMM;
	int iMaxDBE=GetDBE(EMM,NULL)+iExcessDBE;
	int iIonCharge=bPositiveIonMode?+1:-1;
	bool3 bEvenValSum=bAllowRadicalIons?BOOL3_PERH:BOOL3_FALSE;

	if(iPpmReject<iPpmAccept)
		iPpmReject=iPpmAccept;

	ElementIntervalMap EIM;
	Init(EIM,0,EMM);
	GenFormBySpecAlgo BFBSA(EIM,MIM);
	//BrFormByMassAlgo<double> BFBSA(BFWI,mapObserved.begin()->first-0.5,mapObserved.rbegin()->first+0.5);
	//This would be a workaround for the above problem; more time intensive, but correct

	if(pMapMassEMM) pMapMassEMM->clear();

	ElementMultMap EMMFrag;

	for(BFBSA.begin();!BFBSA.end();++BFBSA)
		if(BFBSA.IsValid(bEvenValSum,iMinDiffMaxVal,iMinDiffAtomCount))
		{
			BFBSA.GetElementMultMap(EMMFrag);
			if(GetDBE(EMMFrag,NULL)>iMaxDBE) continue;

			if(pMapMassEMM) 
				setEMM.insert(EMMFrag);
			else
			{
				double dMass=GetNominalMass<double>(EMMFrag)-iIonCharge*c_dMassElectron;
				mapTheoretical[dMass]=1.0;
			}
		}

	return pMapMassEMM?
		ComputeHrMsMsExplainedIntens(mapObserved,setEMM,iPpmAccept,iPpmReject,bPositiveIonMode,pMapMassEMM):
		ComputeHrMsMsExplainedIntens(mapObserved,mapTheoretical,iPpmAccept,iPpmReject,pMapExplained);
}
