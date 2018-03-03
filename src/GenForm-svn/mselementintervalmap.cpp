/*
 *  mselementintervalmap.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mselementintervalmap.h"
#include "mselementintervalmaputil.h"
#include "mschemelement.h"

using namespace std;

ElementIntervalMap::ElementIntervalMap() : map<Element,Interval>()
{
}

ElementIntervalMap::ElementIntervalMap(const ElementIntervalMap& EIM) : map<Element,Interval>(EIM)
{
}

ElementIntervalMap::ElementIntervalMap(const ElementMultMap &EMM)
{
	ElementIntervalMap::iterator itThis=begin();
	for(ElementMultMap::const_iterator itEMM=EMM.begin();itEMM!=EMM.end();itEMM++)
		itThis=this->insert(itThis,make_pair(itEMM->first,Interval(itEMM->second)));
}

void ElementIntervalMap::SetInterval(const Interval& I)
{
	for(ElementIntervalMap::iterator it=begin();it!=end();it++)
		it->second=I;
}

bool ElementIntervalMap::RaiseLowerLimits(const ElementMultMap& EMM)
{
	ElementMultMap::const_iterator itEMM=EMM.begin();

	for(;itEMM!=EMM.end();itEMM++)
	{
		Interval& I=(*this)[itEMM->first];
		I.first=itEMM->second;
		if(!I.IsValid()) return false;
	}

	return true;
}

istream& operator>>(istream &in, ElementIntervalMap &EIM)
{
	string strEIM;

	in >> strEIM;

	if(!in) return in;

	// try(...) here error resulting from  invalid element symbols should be catched
	Init(EIM,strEIM);

	return in;
}

ostream& operator<<(ostream &out, const ElementIntervalMap &EIM)
{
	for(ElementIntervalMap::const_iterator it=EIM.begin();it!=EIM.end();it++)
		out << ChemElement(it->first).GetElementInfo().m_strSymbol
		    << it->second.first << '-' << it->second.second;

	return out;
}

