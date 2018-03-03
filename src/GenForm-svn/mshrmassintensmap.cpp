/*
 *  mshrmassintensmap.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mshrmassintensmap.h"
#include "msdef.h"
#include <iomanip>
#include <cmath>

using namespace std;

HrMassIntensMap::HrMassIntensMap() : 
	GenericMassIntensMap<double,double>()
{
}

HrMassIntensMap::HrMassIntensMap(const GenericMassIntensMap<double,double> &MIM) : 
	GenericMassIntensMap<double,double>(MIM)
{
}
	
const HrMassIntensMap& HrMassIntensMap::operator=(const std::map<double,double> &MIM)
{
	map<double,double>::operator=(MIM);

	return *this;
}

const HrMassIntensMap& HrMassIntensMap::operator=(const std::map<int,double> &MIM)
{
	clear();
	iterator itTarget=begin();

	for(map<int,double>::const_iterator itSource=MIM.begin();itSource!=MIM.end();itSource++)					
		itTarget=insert(itTarget,make_pair((double)(itSource->first),itSource->second));

	return *this;
}
		
bool HrMassIntensMap::IsHighRes() const
{
	for(const_iterator it=begin();it!=end();it++)
		if((double)((int)it->first)!=it->first) return true;

	return false;
}

double HrMassIntensMap::GetSumIntens(double dMinMass,double dMaxMass) const
{
	double dSum=0.0;

	for(const_iterator itPeak=lower_bound(dMinMass);itPeak!=end()&&itPeak->first<=dMaxMass;itPeak++)
	{
		dSum+=itPeak->second;
	//	cout << "\n" << itPeak->first << " " << itPeak->second << " " << dSum << "\n";
	}

	return dSum;
}

double HrMassIntensMap::GetSumIntens(double dMass,int iPPM) const
{
	double dDiff=dMass*(double)iPPM/1000000.0;

	return GetSumIntens(dMass-dDiff,dMass+dDiff);
}

double HrMassIntensMap::GetMaxIntens(double dMinMass, double dMaxMass) const
{
	const_iterator itPeak=GetMaxPeak(dMinMass,dMaxMass);

	return itPeak==end()?0.0:itPeak->second;
}

double HrMassIntensMap::GetMaxIntens(double dMass, int iPPM) const
{
	double dDiff=dMass*(double)iPPM/1000000.0;

	return GetMaxIntens(dMass-dDiff,dMass+dDiff);
}


HrMassIntensMap::const_iterator	HrMassIntensMap::GetMaxPeak(double dMinMass, double dMaxMass) const
{
	const_iterator itMaxPeak=lower_bound(dMinMass);

	if(itMaxPeak->first>dMaxMass) return end();

	for(const_iterator itPeak=lower_bound(dMinMass);itPeak!=end()&&itPeak->first<=dMaxMass;itPeak++)
	{
		if(itMaxPeak->second<itPeak->second)
			itMaxPeak=itPeak;
	//	cout << "\n" << itPeak->first << " " << itPeak->second << " " << itMaxPeak->first << itMaxPeak->second << "\n";
	}

	return itMaxPeak;	
}

HrMassIntensMap::const_iterator	HrMassIntensMap::GetMaxPeak(double dMass, int iPPM) const
{
	double dDiff=dMass*(double)iPPM/1000000.0;

	return GetMaxPeak(dMass-dDiff,dMass+dDiff);
}

double HrMassIntensMap::GetNearestIntens(double dMass, double dMaxDiff)	const
{
	if(empty()) return 0.0;

	const_iterator itPeak=GetNearestPeak(dMass);

	return fabs(itPeak->first-dMass)<=dMaxDiff?itPeak->second:0.0;
}

double HrMassIntensMap::GetNearestIntens(double dMass, int iPPM)		const
{
	double dDiff=dMass*(double)iPPM/1000000.0;

	return GetNearestIntens(dMass,dDiff);
}

HrMassIntensMap::const_iterator	HrMassIntensMap::GetNearestPeak(double dMass) const
{
	if(empty()) return end();

	const_iterator itNext=lower_bound(dMass);

	if(itNext==end()) return --itNext;
	if(itNext==begin()) return itNext;

	const_iterator itPrev=itNext; itPrev--;

	return itNext->first-dMass<dMass-itPrev->first?itNext:itPrev;
}

ostream& operator<<(ostream &out, const HrMassIntensMap &MIM)
{
	for(HrMassIntensMap::const_iterator it=MIM.begin();it!=MIM.end();it++)
		out << setw(10) << setiosflags(ios::fixed) << setprecision(5) << it->first << '\t'
			<< setw(10) << setiosflags(ios::fixed) << setprecision(8) << it->second << '\n';

	return out;
}


istream& operator>>(istream &in, HrMassIntensMap &MIM)
{
	MIM.clear();

	double dMass,dIntens;

	while(true)
	{
		in >> dMass >> dIntens;
		if(in) 
			MIM.AddPeak(dMass,dIntens);
		else
			break;
	}

	return in;
}

void ComputeDecomposition(const HrMassIntensMap& MIM,list<HrMassIntensMap>& listMIM,double dDiversity)
{
	HrMassIntensMap MIMPart;
	HrMassIntensMap::iterator itLast=MIMPart.begin();

	listMIM.clear();

	for(HrMassIntensMap::const_iterator itPeak=MIM.begin();itPeak!=MIM.end();itPeak++)
	{
		if(itPeak!=MIM.begin()&&itPeak->first-MIMPart.GetMaxMass()>dDiversity)
		{
			listMIM.push_back(MIMPart);
			MIMPart.clear();
		}

		itLast=MIMPart.insert(itLast,pair<double,double>(itPeak->first,itPeak->second));
	}

	if(MIMPart.size()) listMIM.push_back(MIMPart);
}

int GetMaxClusterSize(const HrMassIntensMap& MIM)
{
	list<HrMassIntensMap> listMIM;
	int iMaxSize=0;

	ComputeDecomposition(MIM,listMIM);
	for(list<HrMassIntensMap>::const_iterator itMIM=listMIM.begin();itMIM!=listMIM.end();itMIM++)
		iMaxSize=max<int>(iMaxSize,itMIM->size());

	return iMaxSize;
}

HrMassIntensMap& ComputeNeutralLossSpec(const HrMassIntensMap& mapSource,HrMassIntensMap& mapTarget,double dMass)
{
	mapTarget.clear();
	HrMassIntensMap::iterator itLast=mapTarget.begin();

	if(!dMass) dMass=mapSource.GetMaxMass();

	for(HrMassIntensMap::const_reverse_iterator itPeak=mapSource.rbegin();itPeak!=mapSource.rend();itPeak++)
		itLast=mapTarget.insert(itLast,pair<int,double>(dMass-itPeak->first,itPeak->second));
		
	return mapTarget;
}

double CalcDiffPPM(double dCalc, double dExp)
{
	return (dCalc-dExp)*1000000.0/dCalc;
}

double CalcPPM(double dCalc, double dDiff)
{
	return dDiff*1000000.0/dCalc;
}
