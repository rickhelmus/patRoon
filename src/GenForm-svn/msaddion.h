/*
 *  msaddion.h, part of GenForm by M. Meringer, Copyright (C) 2015
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

#ifndef	MS_ADD_ION_H
#define MS_ADD_ION_H

#include "mselementmultmap.h"

// Class for representation of (additional) ions occuring in (ESI) MS
// like for instance the protonated molecular ion, dimeric ion, adduct ions, ...

class AddIon
{
private:
	ElementMultMap		m_emmPlus;
	ElementMultMap		m_emmMinus;
	int					m_iCharge;
	int					m_iMolMult;
	std::string			m_strName;

public:
						AddIon();
						AddIon(const AddIon& AI);
						AddIon(
							const ElementMultMap& emmPlus,
							const ElementMultMap& emmMinus,
							int iCharge=1,int iMolMult=1, 
							const std::string strName="");
						AddIon(const std::string strName);

	const	AddIon&		operator= (const AddIon& AI);
	bool				operator==(const AddIon& AI)			const;
	bool				operator!=(const AddIon& AI)			const;

	const
	ElementMultMap&		GetFormulaPlus()						const
						{ return m_emmPlus; }
	const
	ElementMultMap&		GetFormulaMinus()						const
						{ return m_emmMinus; }
	int					GetMolMult()							const
						{ return m_iMolMult; }
	int					GetCharge()								const
						{ return m_iCharge; }
	const std::string&	GetName()								const
						{ return m_strName; }

	double				CalcMass(double dMoleculeMass)			const;
	double				CalcIonMass(double dMoleculeMass)		const;
	double				CalcMolMass(double dRatioMz)			const;
	double				CalcRatioMZ(double dMoleculeMass)		const;
	double				CalcRatioMZ(const ElementMultMap& EMM)	const;
	const
	ElementMultMap&		CalcEMM(
							const ElementMultMap& emmSource,
							ElementMultMap& emmTarget)			const;
};

std::string ChemName(const AddIon& AI);

#endif
