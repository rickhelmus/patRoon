/*
 *  msisotope.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "msisotope.h"

using namespace std;

const LrMassIntensMap& Init(LrMassIntensMap& MIM,const ElementMultMap& EMM)
{
	bool bFirst=true;
	MIM.clear();

	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
	{
		const ElementInfo& EI=ChemElement(it->first).GetElementInfo();
		for(int iCount=0;iCount<it->second;iCount++)
			if(bFirst){ MIM=EI.m_mapLrIsotope; bFirst=false; }
			else MIM*=EI.m_mapLrIsotope;
	}

	return MIM;
}

const LrMassIntensMap& Init(LrMassIntensMap& MIM,const std::vector<int>& vecCode,const std::vector<int>& vecCount)
{
	bool bFirst=true;
	MIM.clear();

	for(unsigned int iElem=0;iElem<vecCode.size();iElem++)
	{
		const ElementInfo& EI=ChemElement(vecCode[iElem]).GetElementInfo();
		for(int iCount=0;iCount<vecCount[iElem];iCount++)
			if(bFirst){ MIM=EI.m_mapLrIsotope; bFirst=false; }
			else MIM*=EI.m_mapLrIsotope;
	}

	return MIM;	
}

/*
 * Code to read table from SIS
 */
/*
#include <fstream>
#include <sstream>
#include <iomanip>
#include "msexception.h"
#include "mshrmassintensmap.h"

void ReadIsotopeTableSIS(const string& strFileName)
{
	ifstream istFile(strFileName.c_str());
	map<Element, HrMassIntensMap> mapElMim;

	while(istFile)
	{
		string strLine,strLineStart,strElement,strIsotope,strSymbol;
		double dMass,dRatio;

		istFile >> strIsotope;

		if(strIsotope.empty())
			continue;

		if(*strIsotope.rbegin()!=')')
			continue;

		istringstream istLine(strIsotope);

		if(!getline(istLine,strSymbol,'('))
			continue;

		Element E;
		try{E=ChemElement(strSymbol).GetCode();}
		catch(MsError& ME)
		{
			cerr << "unknown element: " << strSymbol << endl;
			continue;//skip unknown element
		}

		//new isotope found

		istFile >> dMass >> dRatio;

		// cout << strSymbol << ' ' << dMass << ' ' << dRatio << endl;

		mapElMim[E].AddPeak(dMass,dRatio);

	}

	map<Element, HrMassIntensMap>::const_iterator itElMim=mapElMim.begin();
	for(;itElMim!=mapElMim.end();itElMim++)
	{
		const HrMassIntensMap& MIM=itElMim->second;
		HrMassIntensMap::const_iterator itMIM=MIM.begin();
		double dSumRatios=0.;
		for(;itMIM!=MIM.end();itMIM++)
		{
			cout << '{' << setw(10) << setiosflags(ios::fixed) << setprecision(6) << itMIM->first
			     << ',' << setw(8) << setiosflags(ios::fixed) << setprecision(5) << itMIM->second/100.0 << "}, ";
			dSumRatios+=itMIM->second;
		}
		cout << "// " << ChemElement(itElMim->first).GetElementInfo().m_strSymbol
			 << (dSumRatios!=100.0?" (*)":"") << endl;
 	}

	for(Element E=e_H;E<e_ElementCount;E++)
	{
		map<Element, HrMassIntensMap>::const_iterator itElMim=mapElMim.find(E);
		cout << (itElMim==mapElMim.end()?0:itElMim->second.size()) << ','
			 << ((int)E%5==4?'\n':' ');
	}
}
*/
