/*
 *  mselementmultmap.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef MS_ELEMENT_MULT_MAP_H
#define MS_ELEMENT_MULT_MAP_H

#include <map>
#include <iostream>

typedef unsigned char Element;
typedef short ElementMult;

class ElementMultMap : public std::map<Element,unsigned short>
{
public:
					ElementMultMap();
					ElementMultMap(const ElementMultMap& EMM);

	bool			operator< (const ElementMultMap& EMM)					const;
	bool			operator<=(const ElementMultMap& EMM)					const;
	bool			operator> (const ElementMultMap& EMM)					const;
	bool			operator>=(const ElementMultMap& EMM)					const;

	bool			IsSubset(const ElementMultMap& EMM)						const;
	const 
	ElementMultMap& operator+=(const ElementMultMap& EMM);
	const 
	ElementMultMap& operator-=(const ElementMultMap& EMM);
	const
	ElementMultMap& operator|=(const ElementMultMap& EMM);			

	inline
	ElementMult		GetMult(const Element& E)								const;
	inline void		SetMult(const Element& E, ElementMult iMult);
	inline void		AddMult(const Element& E, ElementMult iMult);

	ElementMult		GetAtomCount()											const;
};

std::ostream& operator<<(std::ostream &out, const ElementMultMap& EMM);
std::istream& operator>>(std::istream &in, ElementMultMap& EMM);

ElementMultMap operator-(const ElementMultMap& first,const ElementMultMap& second);

/******************************* IMPLEMENTATIONS *****************************/

ElementMult ElementMultMap::GetMult(const Element& E) const
{
	const_iterator it=find(E);
	return it==end()?0:it->second;
}

void ElementMultMap::SetMult(const Element& E, ElementMult iMult)
{
	if(iMult) operator[](E)=iMult;
}

void ElementMultMap::AddMult(const Element& E, ElementMult iMult)
{
	if(iMult) operator[](E)+=iMult;
}

#endif
