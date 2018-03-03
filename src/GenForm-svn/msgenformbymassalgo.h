/*
 *  msgenformbymassalgo.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_GEN_FORM_BY_MASS_ALGO_H
#define MS_GEN_FORM_BY_MASS_ALGO_H

#include "msgenformalgo.h"

//only for the implementation needed:
#include "msrounding.h"
#include "mschemelement.h"
#include <iomanip>

template<class MASS>
class GenFormByMassAlgo : public GenFormAlgo<MASS>
{
private:
	
	std::vector<int>			m_vecMin, m_vecMax;
	std::vector<MASS>			m_vecMass;

	void						ComputeMass();

protected:

	void						InitBegin();
	void						ElementaryStep();
	MASS						m_iMass, m_iMinMass, m_iMaxMass;
	bool						m_bEnd;

public:							
								GenFormByMassAlgo(){}
								GenFormByMassAlgo(
									const ElementIntervalMap& EIM,
									MASS iMinMass, MASS iMaxMass);
								~GenFormByMassAlgo(){}

	virtual void				Init(
									const ElementIntervalMap& EIM,
									MASS iMinMass, MASS iMaxMass);


	virtual void				begin();
	bool						end()									const;
	virtual void				operator++();

	MASS						GetMass()								const;
	
	const std::vector<int>&		GetVecCount()							const
								{ return this->m_vecCount; }
	const std::vector<int>&		GetVecValency()							const
								{ return this->m_vecValency; }
	const std::vector<int>&		GetVecCode()							const
								{ return this->m_vecCode; }
												
	std::ostream&				operator<<(std::ostream& out)			const;
};

template<class MASS>
GenFormByMassAlgo<MASS>::GenFormByMassAlgo(
	const ElementIntervalMap& EIM, MASS iMinMass, MASS iMaxMass)
{
	Init(EIM,iMinMass,iMaxMass);
}

template<class MASS>
void GenFormByMassAlgo<MASS>::Init(
	const ElementIntervalMap& EIM, MASS iMinMass, MASS iMaxMass)
{
	this->m_EIM=EIM;
	m_iMinMass=iMinMass;
	m_iMaxMass=iMaxMass;

	m_vecMass.clear(); m_vecMass.reserve(EIM.size());
	m_vecMin.clear(); m_vecMin.reserve(EIM.size());
	m_vecMax.clear(); m_vecMax.reserve(EIM.size());
	this->m_vecValency.clear(); this->m_vecValency.reserve(EIM.size());
	this->m_vecCode.clear(); this->m_vecValency.reserve(EIM.size());

	for(ElementIntervalMap::const_iterator it=EIM.begin();it!=EIM.end();it++)
	{
		const ElementInfo& EI=ChemElement(it->first).GetElementInfo();
		this->m_vecCode.push_back(it->first);
		this->m_vecValency.push_back(EI.m_iValency);
		m_vecMin.push_back(it->second.first);
		m_vecMax.push_back(it->second.second);
		m_vecMass.push_back(EI.GetNominalMass(MASS()));
	}
}

template<class MASS>
void GenFormByMassAlgo<MASS>::ComputeMass()
{
	m_iMass=MASS();

	for(int iElem=0;iElem<m_vecMass.size();iElem++)
		m_iMass+=m_vecMass[iElem]*this->m_vecCount[iElem];
}

template<class MASS>
void GenFormByMassAlgo<MASS>::ElementaryStep()
{
	for(unsigned int iElem=0;iElem<this->m_vecCount.size();iElem++)
	{
		MASS iMass=Round(m_iMass+m_vecMass[iElem],6);
		if(this->m_vecCount[iElem]<m_vecMax[iElem]&&iMass<=m_iMaxMass)
		{
			m_iMass=iMass;
			this->m_vecCount[iElem]++;
			return;
		}
		else
		{	
			m_iMass=Round(m_iMass-m_vecMass[iElem]*(this->m_vecCount[iElem]-m_vecMin[iElem]),6);
			this->m_vecCount[iElem]=m_vecMin[iElem];
		}
	}

	m_bEnd=true;
}

template<class MASS>
void GenFormByMassAlgo<MASS>::InitBegin()
{
	this->m_vecCount=m_vecMin;
	m_iMass=MASS();
	m_bEnd=false;

	for(unsigned int iElem=0;iElem<m_vecMass.size();iElem++)
	{
		m_iMass+=m_vecMass[iElem]*this->m_vecCount[iElem];
		if(this->m_vecCount[iElem]>m_vecMax[iElem]) m_bEnd=true;
	}

	m_iMass=Round(m_iMass,6);
	if(m_iMass>m_iMaxMass) m_bEnd=true;
}

template<class MASS>
void GenFormByMassAlgo<MASS>::begin()
{
	InitBegin();
	if(m_iMass<m_iMinMass) operator++();
}

template<class MASS>
void GenFormByMassAlgo<MASS>::operator++()
{
	do ElementaryStep();
	while(m_iMass<m_iMinMass&&!m_bEnd);
}

template<class MASS>
bool GenFormByMassAlgo<MASS>::end() const
{
	return m_bEnd;
}

template<class MASS>
MASS GenFormByMassAlgo<MASS>::GetMass() const
{
	return m_iMass;
}

//template<class MASS>
//std::ostream& BrFormByMassAlgo<MASS>::operator<<(std::ostream& out) const
//{
//	int iElem,iWidth=12;
//
//	for(iElem=0;iElem<this->m_BFWI.size();iElem++)
//		out << std::setw(iWidth) << this->m_BFWI.Type(iElem);
//	out << "\n";
//	for(iElem=0;iElem<this->m_vecCount.size();iElem++)
//		out << std::setw(iWidth) << this->m_vecCount[iElem];
//	out << "\n";
//
//	return out;
//}


#endif
