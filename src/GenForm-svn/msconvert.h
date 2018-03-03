/*
 *  msconvert.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_CONVERT_H
#define MS_CONVERT_H

#include "mslrmassintensmap.h"
#include "mshrmassintensmap.h"
	
HrMassIntensMap&	Init(HrMassIntensMap& HR,const LrMassIntensMap& LR);
LrMassIntensMap&	Init(LrMassIntensMap& LR,const HrMassIntensMap& HR, double dShift=0.22);
//dShift: the intensity at integer mass iMass results as sum of intensities for high res
//masses within [iMass-dShift,iMass+1-dShift] 
//a recommended value for this parameter is 0.22; this  value results from the intersection
//point for mass defect of C_n H_2n+2 and O_n

#endif
