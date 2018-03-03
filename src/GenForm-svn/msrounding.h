/*
 *  msrounding.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef MS_ROUNDING_H
#define	MS_ROUNDING_H

#include <cmath>

inline double Round(double	dNumber, int nDigits=6)
{
	double dFactor=1;

	for(int iMult=0;iMult<nDigits;iMult++) dFactor*=10;

	return floor(dFactor*dNumber+.5)/dFactor;
}

inline int Round(int nNumber, int nDigits=5)
{
	return nNumber; 
}

#endif
