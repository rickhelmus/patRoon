/*
 *  mslrmassintensmap.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_LR_MASS_INTENS_MAP_H
#define MS_LR_MASS_INTENS_MAP_H

#include <list>
#include <iostream>
#include "msgenericmassintensmap.h"

/******************************* DEFINES      ********************************/

/******************************* DECLARATIONS ********************************/

/******************************* CLASSES      ********************************/


class LrMassIntensMap : public GenericMassIntensMap<int,double>
{
public:
							LrMassIntensMap();
							LrMassIntensMap(const GenericMassIntensMap<int,double> &MIM);

	const LrMassIntensMap&	operator=(const std::map<int,double> &MIM);
};
	
std::ostream& operator<<(std::ostream &out, const LrMassIntensMap& MIM);
std::istream& operator>>(std::ostream &in, const LrMassIntensMap& MIM);

#endif
