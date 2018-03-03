/*
 *  msknownelements.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_KNOWN_ELEMENTS_H
#define	MS_KNOWN_ELEMENTS_H

#include "mselementinfo.h"
#include <map>

// Chemical elements of periods 1 to 6
enum ElementCode{
		  e_H , e_He, e_Li, e_Be, e_B ,
          e_C , e_N , e_O , e_F , e_Ne,
          e_Na, e_Mg, e_Al, e_Si, e_P ,
          e_S , e_Cl, e_Ar, e_K , e_Ca,
          e_Sc, e_Ti, e_V , e_Cr, e_Mn,
          e_Fe, e_Co, e_Ni, e_Cu, e_Zn,
          e_Ga, e_Ge, e_As, e_Se, e_Br,
          e_Kr, e_Rb, e_Sr, e_Y , e_Zr,
          e_Nb, e_Mo, e_Tc, e_Ru, e_Rh,
          e_Pd, e_Ag, e_Cd, e_In, e_Sn,
          e_Sb, e_Te, e_I , e_Xe, e_Cs,
          e_Ba, e_La, e_Ce, e_Pr, e_Nd,
          e_Pm, e_Sm, e_Eu, e_Gd, e_Tb,
          e_Dy, e_Ho, e_Er, e_Tm, e_Yb,
          e_Lu, e_Hf, e_Ta, e_W , e_Re,
          e_Os, e_Ir, e_Pt, e_Au, e_Hg,
          e_Tl, e_Pb, e_Bi, e_Po, e_At,
          e_Rn, e_D, e_ElementCount};
// Rn completes the 6th period,
// D is just a workaround as long as no general isotopes are available as input

// Mass of an electron
const double c_dMassElectron=5.485799E-004; //5.485799110E-004 in german wiki

class KnownElements
{
private:
	std::map<std::string,int>   m_mapSymbolCode;
	std::vector<ElementInfo>	m_vecElementInfo;

public:
								KnownElements();

	int							GetCode(const std::string &strSymbol);

	const ElementInfo&			GetElementInfo(int iCode)		const
	                            { return m_vecElementInfo[iCode]; }
};

#endif




