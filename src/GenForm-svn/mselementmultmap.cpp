/*
 *  mselementmultmap.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mselementmultmap.h"
#include "mselementmultmaputil.h"

using namespace std;

ElementMultMap::ElementMultMap() : map<Element,unsigned short>()
{
}

ElementMultMap::ElementMultMap(const ElementMultMap &EMM) : map<Element,unsigned short>(EMM)
{
}

bool ElementMultMap::operator<(const ElementMultMap& EMM) const
{
    const_iterator itThis,itEMM;
	for(itThis=begin(),itEMM=EMM.begin();
		itThis!=end()&&itEMM!=EMM.end();itThis++,itEMM++)
	{
		if(itThis->first<itEMM->first) return true;
		else
		if(itThis->first>itEMM->first) return false;
		
		if(itThis->second<itEMM->second) return true;
		else
		if(itThis->second>itEMM->second) return false;
	}

	return itThis==end()&&itEMM!=EMM.end();
}


bool ElementMultMap::operator<=(const ElementMultMap& EMM) const
{
	const_iterator itThis,itEMM;
	for(itThis=begin(),itEMM=EMM.begin();
		itThis!=end()&&itEMM!=EMM.end();itThis++,itEMM++)
	{
		if(itThis->first<itEMM->first) return true;
		else
		if(itThis->first>itEMM->first) return false;
		
		if(itThis->second<itEMM->second) return true;
		else
		if(itThis->second>itEMM->second) return false;
	}

	return itThis==end();			
};

bool ElementMultMap::operator>(const ElementMultMap& EMM) const
{
	return !operator<=(EMM);
}

bool ElementMultMap::operator>=(const ElementMultMap& EMM) const
{
	return !operator<(EMM);
}

//returns true, if *this is subset of EMM, read as 'IsSubsetOf'
bool ElementMultMap::IsSubset(const ElementMultMap& EMM) const
{	//not yet tested
    const_iterator itThis=begin(),itEMM=EMM.begin();
	while(itThis!=end()&&itEMM!=EMM.end())
	{
		if(itThis->first<itEMM->first) return false;

		if(itThis->first>itEMM->first)
		{ itEMM++; continue; }

		if(itThis->second>itEMM->second) return false;

		itThis++; itEMM++;
	}

	return itThis==end();
}

ElementMult ElementMultMap::GetAtomCount() const
{
	ElementMult nCount=0;

	for(const_iterator it=begin();it!=end();it++)
		nCount+=it->second;

	return nCount;
}

const ElementMultMap& ElementMultMap::operator+=(const ElementMultMap& EMM)
{	// not yet optimized
	for(const_iterator it=EMM.begin();it!=EMM.end();it++)
		operator[](it->first)+=it->second;

	return *this;
}

const ElementMultMap& ElementMultMap::operator-=(const ElementMultMap& EMM)
{	// not yet optimized
	for(const_iterator it=EMM.begin();it!=EMM.end();it++)
		operator[](it->first)-=it->second;

	for(iterator itTest=begin();itTest!=end();)
	{
		iterator itDel=itTest++;
		if(!itDel->second) erase(itDel);
	}

	return *this;
}
	
const ElementMultMap& ElementMultMap::operator|=(const ElementMultMap& EMM)
{	// not yet optimized
	for(const_iterator it=EMM.begin();it!=EMM.end();it++)
		operator[](it->first)=std::max(operator[](it->first),it->second);

	return *this;
}

ostream& operator<<(ostream &out, const ElementMultMap &EMM)
{
	for(ElementMultMap::const_iterator it=EMM.begin();it!=EMM.end();it++)
		out << ChemElement(it->first).GetElementInfo().m_strSymbol << it->second;

	return out;
}

istream& operator>>(istream &in, ElementMultMap &EMM)
{
	string strEMM;

	in >> strEMM;

	if(!in) return in;

	// try(...) here error resulting from  invalid element symbols should be catched
	Init(EMM,strEMM);

	return in;
}

ElementMultMap operator-(const ElementMultMap& first,const ElementMultMap& second)
{
	return ElementMultMap(first)-=second;
}

