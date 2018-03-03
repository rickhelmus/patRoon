/*
 *  msaddion.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_CHEM_ELEMENT_H
#define	MS_CHEM_ELEMENT_H

#include "msknownelements.h"

class ChemElement
{
protected:
	static KnownElements s_KnownElements;

	unsigned char		m_cCode;

public:
						ChemElement(const std::string& strSymbol)
						{ m_cCode = s_KnownElements.GetCode(strSymbol); }

						ChemElement(unsigned char cCode)
						{ m_cCode = cCode;}
						
	unsigned char		GetCode()									const
						{ return m_cCode; }
	const ElementInfo&	GetElementInfo()							const
						{return s_KnownElements.GetElementInfo(m_cCode); }
};

#endif
