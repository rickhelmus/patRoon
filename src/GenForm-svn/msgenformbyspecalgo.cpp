/*
 *  msgenformbyspecalgo.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "msgenformbyspecalgo.h"

using namespace std;

GenFormBySpecAlgo::GenFormBySpecAlgo(
	const ElementIntervalMap& EIM, const LrMassIntensMap& MIM)
{
	Init(EIM,MIM);
}

void GenFormBySpecAlgo::Init(
	const ElementIntervalMap& EIM, const LrMassIntensMap& MIM)
{
	GenFormByMassAlgo<int>::Init(EIM,MIM.GetMinMass(),MIM.GetMaxMass());
	
	m_vecNullIntens.assign(m_iMaxMass+1,true);
	for(LrMassIntensMap::const_iterator it=MIM.begin();it!=MIM.end();it++)
		if(it->second) m_vecNullIntens[it->first]=false;
}

void GenFormBySpecAlgo::Init(
	const ElementIntervalMap& EIM,int iMinMass, int iMaxMass)
{
	GenFormByMassAlgo<int>::Init(EIM,iMinMass,iMaxMass);
	
	m_vecNullIntens.assign(m_iMaxMass+1,true);
	for(int iMass=m_iMinMass;iMass<=m_iMaxMass;iMass++) 
		m_vecNullIntens[iMass]=false;
}

void GenFormBySpecAlgo::begin()
{
	InitBegin();
	if(m_vecNullIntens[m_iMass]) operator++();
}


void GenFormBySpecAlgo::operator++()
{
	do ElementaryStep();
	while(m_vecNullIntens[m_iMass]&&!m_bEnd);
}
