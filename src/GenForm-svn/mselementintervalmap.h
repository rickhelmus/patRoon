/*
 *  mselementintervalmap.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef MS_ELEMENT_INTERVAL_MAP_H
#define MS_ELEMENT_INTERVAL_MAP_H

#include <map>
#include <iostream>
#include "msgenericinterval.h"
#include "mselementmultmap.h"

typedef GenericInterval<ElementMult> Interval;

class ElementIntervalMap : public std::map<Element, Interval>
{
public:
					ElementIntervalMap();
					ElementIntervalMap(const ElementIntervalMap& EIM);
					ElementIntervalMap(const ElementMultMap& EMM);

	inline
	Interval		GetInterval(const Element& E)							const;
	inline void		SetInterval(const Element& E, const Interval& I);
	void			SetInterval(const Interval& I);
	bool			RaiseLowerLimits(const ElementMultMap& EMM);
};

std::ostream& operator<<(std::ostream &out, const ElementIntervalMap& EIM);
std::istream& operator>>(std::istream &in, ElementIntervalMap& EIM);

Interval ElementIntervalMap::GetInterval(const Element& E) const
{
	const_iterator it=find(E);
	return it==end()?Interval():it->second;
}

void ElementIntervalMap::SetInterval(const Element& E, const Interval& I)
{
	if(I.IsValid()&&I!=Interval()) operator[](E)=I;
}

#endif
