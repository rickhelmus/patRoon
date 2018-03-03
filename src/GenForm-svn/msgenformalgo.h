/*
 *  msgenformalgo.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_GEN_FORM_ALGO_H
#define MS_GEN_FORM_ALGO_H

#include <vector>
#include "mselementmultmap.h"
#include "mselementintervalmap.h"
#include "msdef.h"

template<class MASS>
class GenFormAlgo
{
protected:
	std::vector<int>			m_vecCount, m_vecValency, m_vecCode;
	ElementIntervalMap			m_EIM;
public:
	virtual void				begin()											=0;
	virtual bool				end()									const	=0;
	virtual void				operator++()									=0;

	ElementMultMap&				GetElementMultMap(ElementMultMap& EMM)	const;
	virtual MASS				GetMass()								const	=0;

	//There are three criteria for existence of a connected multigraph
	// (GR1) valency sum = 0 mod 2
	// (GR2) valency sum - 2 * max valency >= 0
	// (CON) valency sum - 2 * atom count +2 >= 0
	//IsConnected() only tests (GR2) and (CON) and stores the parity in 
	//bOddElectronIon (true, if valency sum is even, false otherwise)
	//IsValid is more sophisticated, as you can specify the parity
	//for (GR1) and especially lower bounds iMinDiffMaxVal for (GR2)
	//and iMinDiffAtomCount for (GR3) instead of 0
	bool						IsValid(
									bool3 b3EvenValSum=BOOL3_TRUE,
									int iMinDiffMaxVal=0, 
									int iMinDiffAtomCount=0)			const;

	//Obsolete versions of IsValid()
	bool						IsConnected(bool& bOddElectronIon)		const;
	bool						IsConnected()							const
								{ bool b; return IsConnected(b); }

};

template<class MASS>
bool GenFormAlgo<MASS>::IsValid(bool3 b3EvenValSum,int iMinDiffMaxVal,int iMinDiffAtomCount) const
{
	int iAtomCount=0,iValencyCount=0,iMaxValency=0;

	for(unsigned int iElem=0;iElem<m_vecCount.size();iElem++)
	{
		iAtomCount+=m_vecCount[iElem];
		iValencyCount+=m_vecCount[iElem]*m_vecValency[iElem];
		if(m_vecCount[iElem]&&m_vecValency[iElem]>iMaxValency)
			iMaxValency=m_vecValency[iElem];
    }

	bool bParityOK=(b3EvenValSum==BOOL3_PERH)?true:iValencyCount%2!=b3EvenValSum;

	return bParityOK
		&& iValencyCount-iMaxValency*2>=iMinDiffMaxVal
		&& iValencyCount-iAtomCount*2+2>=iMinDiffAtomCount;
}

template<class MASS>
bool GenFormAlgo<MASS>::IsConnected(bool& bOddElectronIon) const
{
	int iAtomCount=0,iValencyCount=0,iMaxValency=0;

	for(int iElem=0;iElem<m_vecCount.size();iElem++)
	{
		iAtomCount+=m_vecCount[iElem];
		iValencyCount+=m_vecCount[iElem]*m_vecValency[iElem];
		if(m_vecCount[iElem]&&m_vecValency[iElem]>iMaxValency)
			iMaxValency=m_vecValency[iElem];
    }   

	bOddElectronIon=iValencyCount%2?false:true;

	return iValencyCount>=iMaxValency*2 && iValencyCount>=iAtomCount*2-2;
}


template<class MASS>
ElementMultMap& GenFormAlgo<MASS>::GetElementMultMap(ElementMultMap& EMM) const
{
	EMM.clear();
	ElementMultMap::iterator itEMM=EMM.begin();

	for(unsigned int iElem=0;iElem<m_vecCount.size();iElem++)
		if(m_vecCount[iElem])
			itEMM=EMM.insert(itEMM,std::pair<Element,unsigned short>(m_vecCode[iElem],m_vecCount[iElem]));

	return EMM;
}

#endif
