/*
 *  mselementinfo.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef MS_ELEMENT_INFO_H
#define MS_ELEMENT_INFO_H

#include <map>
#include <vector>
#include <string>
#include "mslrmassintensmap.h"
#include "mshrmassintensmap.h"

class ElementInfo
{
public:
	std::string				m_strSymbol;
	int						m_iValency;

	int						m_iNominalMass;
	double					m_dNominalMass;

	LrMassIntensMap			m_mapLrIsotope;
	HrMassIntensMap			m_mapHrIsotope;

	std::vector<int>		m_vecStateVal;
	std::vector<int>		m_vecStateLon;
	std::vector<int>		m_vecStateRad;
	std::vector<int>		m_vecStateCha;

	int						GetNominalMass(int)	const
							{ return m_iNominalMass; }
	double					GetNominalMass(double) const
							{ return m_dNominalMass; }
};

#endif
