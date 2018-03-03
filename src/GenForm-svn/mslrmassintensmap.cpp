/*
 *  mslrmassintensmap.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include "mslrmassintensmap.h"
#include "msdef.h"
#include <iomanip>
#include <cmath>

using namespace std;

LrMassIntensMap::LrMassIntensMap() : 
	GenericMassIntensMap<int,double>()
{
}

LrMassIntensMap::LrMassIntensMap(const GenericMassIntensMap<int,double> &MIM) : 
	GenericMassIntensMap<int,double>(MIM)
{
}
	
const LrMassIntensMap& LrMassIntensMap::operator=(const std::map<int,double> &MIM)
{
	map<int,double>::operator=(MIM);

	return *this;
}

ostream& operator<<(ostream &out, const LrMassIntensMap &MIM)
{
	for(LrMassIntensMap::const_iterator it=MIM.begin();it!=MIM.end();it++)
		out << setw(4) << it->first << '\t'
			<< setw(10) << setiosflags(ios::fixed) << setprecision(8) << it->second << '\n';

	return out;
}
