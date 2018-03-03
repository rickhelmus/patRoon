/*
 *  mselementintervalmaputil.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mselementintervalmaputil.h"
#include "msstringconversion.h"
#include "msexception.h"
#include <sstream>
#include <stdlib.h>
#include <limits>
#include "mschemelement.h"

using namespace std;

const ElementIntervalMap& Init(ElementIntervalMap& EIM,const std::string& strFormula)
{
	unsigned int iPos=0;
	string strAlpha,strDigitMin,strDigitMax;
	bool bOpenIntervalFound,bClosedIntervalFound;
	Element E;
	Interval I;

	EIM.clear();

	while(iPos<strFormula.length())
	{
		bOpenIntervalFound=bClosedIntervalFound=false;

		while(iPos<strFormula.length())
		{
			if(isalnum(strFormula[iPos])) break;
			iPos++;
		}

		strAlpha="";
		while(iPos<strFormula.length())
		{	
			if(!isalpha(strFormula[iPos])||
				(isupper(strFormula[iPos])&&!strAlpha.empty())||
				strAlpha.length()>=2) 
				break;
			strAlpha+=strFormula[iPos];
			iPos++;
		}

		if(iPos+1<strFormula.length())
		{
			if(strFormula[iPos]=='>'&&strFormula[iPos+1]=='=')
			{
				bOpenIntervalFound=true;
				iPos+=2;
			}
		}

		strDigitMin="";
		while(iPos<strFormula.length())
		{
			if(!isdigit(strFormula[iPos])) break;
			strDigitMin+=strFormula[iPos];
			iPos++;
		}

		if(iPos<strFormula.length()&&!bOpenIntervalFound)
		{
			if(strFormula[iPos]=='-')
			{
				bClosedIntervalFound=true;
				iPos++;
			}
		}

		strDigitMax="";
		while(iPos<strFormula.length()&&bClosedIntervalFound)
		{
			if(!isdigit(strFormula[iPos]))
				break;
			strDigitMax+=strFormula[iPos];
			iPos++;
		}

		if(!strAlpha.empty())
		{
			strAlpha[0]=toupper(strAlpha[0]);

			if(strDigitMin.empty())
				I=Interval(1,1);
			else
			if(strDigitMax.empty()&&!bOpenIntervalFound)
				I=Interval(atoi(strDigitMin.c_str()));
			else
			if(strDigitMax.empty()&&bOpenIntervalFound)
				I=Interval(atoi(strDigitMin.c_str()),numeric_limits<ElementMult>::max());
			else
				I=Interval(atoi(strDigitMin.c_str()),atoi(strDigitMax.c_str()));

			if(I.IsValid()&&I.first>=0&&I.second>=1)
			{
				try{E=ChemElement(strAlpha).GetCode();}
				catch(MsError& ME)
				{
					throw(ME);
				}

				EIM[E]=I;
			}
		}
	}

	return EIM;
}

const ElementIntervalMap& Init(ElementIntervalMap& EIM,int nLowerBound,const ElementMultMap& EMM)
{
	EIM.clear();

	ElementIntervalMap::iterator itEIM=EIM.begin();

	for(ElementMultMap::const_iterator itEMM=EMM.begin();itEMM!=EMM.end();itEMM++)
		if(nLowerBound<=itEMM->second)
			itEIM=EIM.insert(itEIM,make_pair(itEMM->first,Interval(nLowerBound,itEMM->second)));

	return EIM;
}

