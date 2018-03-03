/*
 *  mselementmultmaputil.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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
#include "msstringconversion.h"
#include "msexception.h"

#include <sstream>
#include <stdlib.h>

using namespace std;

const ElementMultMap& Init(ElementMultMap& EMM,const std::string& strFormula)
{
	unsigned int iPos=0;
	string strAlpha,strDigit;
	Element E;

	EMM.clear();

	while(iPos<strFormula.length())
	{
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

		strDigit="";
		while(iPos<strFormula.length())
		{
			if(!isdigit(strFormula[iPos])) break;
			strDigit+=strFormula[iPos];
			iPos++;
		}

		if(!strAlpha.empty())
		{
			strAlpha[0]=toupper(strAlpha[0]);

			int iMult=strDigit.empty()?1:atoi(strDigit.c_str());

			if(iMult)
			{
				try{E=ChemElement(strAlpha).GetCode();}
				catch(MsError& ME)
				{
					throw(ME);
				}

				EMM[E]+=iMult;
			}
		}
	}

	return EMM;
}

bool LessChem(const ElementMultMap& EMM1,const ElementMultMap& EMM2)
{
	//alphabetical with the exception that C and H are first.
	//smallest formula is the one with least C atoms,
	//if C counts are equal, then the one with least H atoms,
	//further elements follow in the sequence of their atomic number.

	ElementMultMap::const_iterator itC1=EMM1.find(e_C);
	ElementMultMap::const_iterator itC2=EMM2.find(e_C);

	if(itC1!=EMM1.end()&&itC2!=EMM2.end())
	{
		if(itC1->second<itC2->second) return true;
		else
		if(itC1->second>itC2->second) return false;
	}
	else
	if(itC1==EMM1.end()&&itC2!=EMM2.end()) return true;
	else
	if(itC1!=EMM1.end()&&itC2==EMM2.end()) return false;

	ElementMultMap::const_iterator itH1=EMM1.find(e_H);
	ElementMultMap::const_iterator itH2=EMM2.find(e_H);

	if(itH1!=EMM1.end()&&itH2!=EMM2.end())
	{
		if(itH1->second<itH2->second) return true;
		else
		if(itH1->second>itH2->second) return false;
	}
	else
	if(itH1==EMM1.end()&&itH2!=EMM2.end()) return true;
	else
	if(itH1!=EMM1.end()&&itH2==EMM2.end()) return false;


	return EMM1<EMM2;
}

bool HasHeteroAtom(const ElementMultMap& EMM)
{
	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
		if(it->first!=e_H&&it->first!=e_C) return true;

	return false;
}

string ChemName(const ElementMultMap& EMM)
{
	//MM, 20.04.10: type conversion to short introduced before printing
	//in order to get a better output in case negative element multiplicities occur

	ostringstream out;

	ElementMultMap::const_iterator itC=EMM.find(e_C);
	if(itC!=EMM.end())
		out << "C" << IntToFormulaString((short)itC->second);

	ElementMultMap::const_iterator itH=EMM.find(e_H);
	if(itH!=EMM.end())
		out << "H" << IntToFormulaString((short)itH->second);

	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
		if(it!=itC&&it!=itH)
			out << ChemElement(it->first).GetElementInfo().m_strSymbol << IntToFormulaString((short)it->second);

	return out.str();
}

string TexName(const ElementMultMap& EMM)
{
	ostringstream out;

	out << "$\\chem{";

	ElementMultMap::const_iterator itC=EMM.find(e_C);
	if(itC!=EMM.end())
		out << "C_{" << IntToFormulaString((short)itC->second) << "}";

	ElementMultMap::const_iterator itH=EMM.find(e_H);
	if(itH!=EMM.end())
		out << "H_{" << IntToFormulaString(itH->second) << "}";

	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
		if(it!=itC&&it!=itH)
			out << ChemElement(it->first).GetElementInfo().m_strSymbol << "_{" << IntToFormulaString(it->second) << "}";

	out << "}$";

	return out.str();
}

int GetDBE(const ElementMultMap& EMM,bool* pInteger)
{
	int iDBE=2;

	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
		iDBE+=it->second*(ChemElement(it->first).GetElementInfo().m_iValency-2);

	if(pInteger) *pInteger=iDBE%2?false:true;

	return iDBE/2;
}

double CalcDBE(const ElementMultMap& EMM)
{
	int iDBE=2;

	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
		iDBE+=it->second*(ChemElement(it->first).GetElementInfo().m_iValency-2);

	return iDBE/2.0;
}

bool IsGraphicalMultVal(const ElementMultMap& EMM)
{
	int iSumVal=0,iMaxVal=0;
	int nP=0,nS=0,iNumAtoms=0;

	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
	{
		iNumAtoms+=it->second;

		if(it->first==e_P)
			nP=it->second;
		else
		if(it->first==e_S)
			nS=it->second;
		else
		{
			int iVal=ChemElement(it->first).GetElementInfo().m_iValency;
			iSumVal+=it->second*iVal;
			iMaxVal=max<int>(iMaxVal,iVal);
		}
	}

	if( (iSumVal+nP*3)%2 != 0 ) return false; // Check criterion (i) at first

	if(nP==0)
	{
		if(nS==0)
		{
			if( iSumVal - 2*iMaxVal >= 0 &&
				iSumVal - 2*iNumAtoms +2 >= 0 )
				return true;
		}
		else
		{
			for(int iMaxValS=2;iMaxValS<=6;iMaxValS+=2) // loop over the maximum valencies for S
				for(int iSumValS=iMaxValS+(nS-1)*2;iSumValS<=nS*iMaxValS;iSumValS+=2) // loop over sum of valencies for S
				{
					int iSumValAll=iSumVal+iSumValS;
					int iMaxValAll=max<int>(iMaxVal,iMaxValS);

					if( iSumValAll - 2*iMaxValAll >= 0 &&
						iSumValAll - 2*iNumAtoms +2 >= 0 )
						return true;
				}
		}
	}
	else//(nP>=1)
	{
		if(nS==0)
		{
			for(int iMaxValP=3;iMaxValP<=5;iMaxValP+=2) // loop over the maximum valencies for P
				for(int iSumValP=iMaxValP+(nP-1)*2;iSumValP<=nP*iMaxValP;iSumValP+=2) // loop over sum of valencies for P
				{
					int iSumValAll=iSumVal+iSumValP;
					int iMaxValAll=max<int>(iMaxVal,iMaxValP);

					if( iSumValAll - 2*iMaxValAll >= 0 &&
						iSumValAll - 2*iNumAtoms +2 >= 0 )
						return true;
				}
		}
		else
		{
			for(int iMaxValP=3;iMaxValP<=5;iMaxValP+=2) // loop over the maximum valencies for P
				for(int iSumValP=iMaxValP+(nP-1)*2;iSumValP<=nP*iMaxValP;iSumValP+=2) // loop over sum of valencies for P
					for(int iMaxValS=2;iMaxValS<=6;iMaxValS+=2) // loop over the maximum valencies for S
						for(int iSumValS=iMaxValS+(nS-1)*2;iSumValS<=nS*iMaxValS;iSumValS+=2) // loop over sum of valencies for S
						{
							int iSumValAll=iSumVal+iSumValP+iSumValS;
							int iMaxValAll=max<int>(iMaxVal,max<int>(iMaxValP,iMaxValS));

							if( iSumValAll - 2*iMaxValAll >= 0 &&
								iSumValAll - 2*iNumAtoms +2 >= 0 )
								return true;
						}
		}
	}

	return false;
}


