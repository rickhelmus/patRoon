/*
 *  msstringconversion.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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
#include "msstringconversion.h"

using namespace std;

string IntToFormulaString(int iMult)
{
	ostringstream ost;
	if(iMult!=1) ost << iMult;
	return ost.str();
}

string IntToChargeString(int nCharge,bool bSignumFirst)
{
	ostringstream ost;

	if(nCharge>0)
	{
		if(bSignumFirst) ost << '+';
		if(nCharge>1) ost << nCharge;
		if(!bSignumFirst) ost << '+';
	}
	else
	if(nCharge<0)
	{
		if(bSignumFirst) ost << '-';
		if(nCharge<-1) ost << -nCharge;
		if(!bSignumFirst) ost << '-';
	}

	return ost.str();
}
