/*
 *  msaddion.cpp, part of GenForm by M. Meringer, Copyright (C) 2015
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

#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "mselementmultmaputil.h"
#include "msstringconversion.h"
#include "msexception.h"
#include "msaddion.h"

using namespace std;

AddIon::AddIon() :
	m_iCharge(0),
	m_iMolMult(0)
{
}

AddIon::AddIon(
	const ElementMultMap& emmPlus,
	const ElementMultMap& emmMinus,
	int iCharge,int iMolMult,const string strName) :
	m_emmPlus(emmPlus),
	m_emmMinus(emmMinus),
	m_iCharge(iCharge),
	m_iMolMult(iMolMult),
	m_strName(strName)
{
}

AddIon::AddIon(const AddIon& AI)
{
	operator=(AI);
}

AddIon::AddIon(const string strName)
{
	const int nCount=53;
	// adduct info adapted from
	// http://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/
	const char* pAdductInfo[nCount][5] =
	{
		{"M+3H", "H3", "", "+3", "1"},
		{"M+2H+Na", "H2Na", "", "+3", "1"},
		{"M+H+2Na", "HNa2", "", "+3", "1"},
		{"M+3Na", "Na3", "", "+3", "1"},
		{"M+2H", "H2", "", "+2", "1"},
		{"M+H+NH4", "H5N", "", "+2", "1"},
		{"M+H+Na", " HNa", "", "+2", "1"},
		{"M+H+K", "HK", "", "+2", "1"},
		{"M+ACN+2H", "C2H5N", "", "+2", "1"},//ACN=C2H3N
		{"M+2Na", "Na2", "", "+2", "1"},
		{"M+2ACN+2H", "C4H8N2", "", "+2", "1"},
		{"M+3ACN+2H", "C6H11N3", "", "+2", "1"},
		{"M+H", "H", "", "+1", "1"},
		{"M+NH4", "H4N", "", "+1", "1"},
		{"M+Na", "Na", "", "+1", "1"},
		{"M+CH3OH+H", "CH5O", "", "+1", "1"},
		{"M+K", "K", "", "+1", "1"},
		{"M+ACN+H", "C2H4N", "", "+1", "1"},
		{"M+2Na-H", "Na2", "H", "+1", "1"},
		{"M+IsoProp+H", "C3H10O", "", "+1", "1"}, //IsoProp=C3H9O(?)
		{"M+ACN+Na", "C2H3NNa", "", "+1", "1"},
		{"M+2K-H", "K2", "H", "+1", "1"},
		{"M+DMSO+H", "C2H7OS", "", "+1", "1"}, //DMSO=C2H6OS
		{"M+2ACN+H", "C4H7N2", "", "+1", "1"},
		{"M+IsoProp+Na+H", "C3H10ONa", "", "+1", "1"},
		{"2M+H", "H", "", "+1", "2"},
		{"2M+NH4", "H4N", "", "+1", "2"},
		{"2M+Na", "Na", "", "+1", "2"},
	    {"2M+K", "K", "", "+1", "2"},
		{"2M+ACN+H", "C2H4N", "", "+1", "2"},
		{"2M+ACN+Na", "C2H3NNa", "", "+1", "2"},

		{"M-3H", "", "H3", "-3", "1"},
		{"M-2H", "", "H2", "-2", "1"},
		{"M-H2O-H", "", "H3", "-1", "1"},
		{"M-H", "", "H", "-1", "1"},
		{"M+Na-2H", "Na", "H2", "-1", "1"},
		{"M+Cl", "Cl", "", "-1", "1"},
		{"M+K-2H", "K", "H2", "-1", "1"},
		{"M+FA-H", "CH2O2", "H", "-1", "1"}, //FA=CH2O2
		{"M+Hac-H", "C2H4O2", "H", "-1", "1"}, //HAc=C2H4O2
		{"M+Br", "Br", "", "-1", "1"},
		{"M+TFA-H", "C2HF3O2", "H", "-1", "1"}, //TFA=C2HF3O2
		{"2M-H", "", "H", "-1", "2"},
		{"2M+FA-H", "", "H", "-1", "2"},
		{"2M+Hac-H", "C2H4O2", "H", "-1", "2"},
		{"3M-H", "", "H", "-1", "3"},

		{"M-e", "", "", "+1", "1"},
		{"M+e", "", "H", "-1", "1"},

		{"-e", "", "", "+1", "1"},
		{"+e", "", "", "-1", "1"},
		{"+H", "H", "", "+1", "1"},
		{"-H", "", "H", "-1", "1"},
		{"+Na", "Na", "", "+1", "1"},
	};

	for(int nIdx=0;nIdx!=nCount;nIdx++)
		if(strName==pAdductInfo[nIdx][0])
		{
			Init(m_emmPlus,pAdductInfo[nIdx][1]);
			Init(m_emmMinus,pAdductInfo[nIdx][2]);
			m_iCharge=atoi(pAdductInfo[nIdx][3]);
			m_iMolMult=atoi(pAdductInfo[nIdx][4]);
			m_strName=strName;
			return;
		}

	string s="error: unknown adduct type '" + strName +"'";
	MsError Me(s); throw(Me);
}

const AddIon& AddIon::operator=(const AddIon& AI)
{
	m_emmPlus=AI.m_emmPlus;
	m_emmMinus=AI.m_emmMinus;
	m_iCharge=AI.m_iCharge;
	m_iMolMult=AI.m_iMolMult;
	m_strName=AI.m_strName;

	return *this;
}

bool AddIon::operator==(const AddIon& AI) const
{
	return
		m_emmPlus==AI.m_emmPlus&&
		m_emmMinus==AI.m_emmMinus&&
		m_iCharge==AI.m_iCharge&&
		m_iMolMult==AI.m_iMolMult&&
		m_strName==AI.m_strName;
}

bool AddIon::operator!=(const AddIon& AI) const
{
	return !operator==(AI);
}

double AddIon::CalcMass(double dMoleculeMass) const
{
	return dMoleculeMass*(double)m_iMolMult
		+GetNominalMass<double>(m_emmPlus)
		-GetNominalMass<double>(m_emmMinus);
}

// method that also takes electron masses into account
double AddIon::CalcIonMass(double dMoleculeMass) const
{
	return dMoleculeMass*(double)m_iMolMult
		+GetNominalMass<double>(m_emmPlus)
		-GetNominalMass<double>(m_emmMinus)
		-m_iCharge*c_dMassElectron;
}

// method used for generation of molecular formulas
// from accurate mass measurements
double AddIon::CalcMolMass(double dRatioMz) const
{
	return (dRatioMz*fabs((double)m_iCharge)
		+m_iCharge*c_dMassElectron
		-GetNominalMass<double>(m_emmPlus)
		+GetNominalMass<double>(m_emmMinus))/m_iMolMult;
}

double AddIon::CalcRatioMZ(double dMoleculeMass) const
{
	return CalcMass(dMoleculeMass)/(double)m_iCharge;
}

double AddIon::CalcRatioMZ(const ElementMultMap& EMM) const
{
	return CalcRatioMZ(GetNominalMass<double>(EMM));
}

const ElementMultMap& AddIon::CalcEMM(const ElementMultMap& emmSource, ElementMultMap& emmTarget) const
{
	emmTarget=emmSource;
	emmTarget+=m_emmPlus;
	emmTarget-=m_emmMinus;
	return emmTarget;
}

string ChemName(const AddIon& AI)
{
	ostringstream out;

	out << '[' 
		<< IntToFormulaString(AI.GetMolMult())
		<< 'M' 
		<< (AI.GetFormulaPlus().empty()?"":("+"))
		<< ChemName(AI.GetFormulaPlus())
		<< (AI.GetFormulaMinus().empty()?"":("-"))
		<< ChemName(AI.GetFormulaMinus())
		<< ']'
		<< IntToChargeString(AI.GetCharge());

	return out.str();
}
