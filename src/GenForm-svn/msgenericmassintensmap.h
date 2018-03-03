/*
 *  msgenericmassintensmap.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_GENERIC_MASS_INTENS_MAP_H
#define MS_GENERIC_MASS_INTENS_MAP_H

#include <map>
#include <list>

template<class M,class I> class GenericMassIntensMap : public std::map<M,I>
{
public:
								GenericMassIntensMap();
								GenericMassIntensMap(
									const std::map<M,I> &MIM);

	const GenericMassIntensMap<M,I>&	operator=(const std::map<M,I> &MIM);
	const GenericMassIntensMap<M,I>&	operator+=(const std::map<M,I> &MIM);
	const GenericMassIntensMap<M,I>&	operator*=(const std::map<M,I> &MIM);
	const GenericMassIntensMap<M,I>&	operator*=(I dFactor);

	inline I					GetIntens(M nMass)					const;
	inline M					GetMinMass()						const;	
	inline M					GetMaxMass()						const;

	typename GenericMassIntensMap<M,I>::const_iterator
								GetBasePeak()						const;
	typename GenericMassIntensMap<M,I>::iterator
								GetBasePeak();

	typename GenericMassIntensMap<M,I>::const_iterator
								GetNearestPeak(M nMass)				const;

	I							GetMaxIntens()						const;
	I							GetTotalIntens()					const;
	I							GetTotalIntensQuad()				const;
	I							GetTotalIntensMassProd()			const;
	I							GetTotalIntensQuadMass()			const;

	inline void					AddPeak(M nMass,I dIntens);

	const GenericMassIntensMap&	Normalize(I dIntens=1.0,bool bMax=true);
	const GenericMassIntensMap&	Round(int nDigits=6);
	const GenericMassIntensMap&	CutNoise(I dCut=0.0);
	const GenericMassIntensMap&	Convolution(
								const std::map<M,I> &MIM1, 
								const std::map<M,I> &MIM2);
	const GenericMassIntensMap&	DelMasses(const std::map<M,I> &MIM);
	const GenericMassIntensMap&	Maximum(const std::map<M,I> &MIM);
	const GenericMassIntensMap&	ShiftMasses(const M &dShift);
};

template<class M,class I>
GenericMassIntensMap<M,I>::GenericMassIntensMap() : std::map<M,I>()
{
}

template<class M,class I>
GenericMassIntensMap<M,I>::GenericMassIntensMap(const std::map<M,I> &MIM) : std::map<M,I>(MIM)
{
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::operator=(const std::map<M,I> &MIM)
{
	std::map<M,I>::operator=(MIM);

	return *this;
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::operator+=(const std::map<M,I> &MIM)
{
	for(typename std::map<M,I>::const_iterator itMIM=MIM.begin();itMIM!=MIM.end();itMIM++)
		this->operator[](itMIM->first)+=itMIM->second;
	//could be improved by running through both sequences

	return *this;
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::operator*=(const std::map<M,I> &MIM)
{
	GenericMassIntensMap<M,I> MIM2(*this);

	return this==&MIM?Convolution(MIM2,MIM2):Convolution(MIM,MIM2);
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::operator*=(I dFactor)
{
	for(typename GenericMassIntensMap<M,I>::iterator it=this->begin();it!=this->end();it++)
		it->second*=dFactor;

	return *this;	
}

template<class M,class I>
I GenericMassIntensMap<M,I>::GetMaxIntens() const
{
	return this->empty()?0.0:this->GetBasePeak()->second;
}

template<class M,class I>
I GenericMassIntensMap<M,I>::GetTotalIntens() const
{
	I dSum=0.0;

	for(typename GenericMassIntensMap<M,I>::const_iterator it=this->begin();it!=this->end();it++)
		dSum+=it->second;

	return dSum;
}

template<class M,class I>
I GenericMassIntensMap<M,I>::GetTotalIntensQuad() const
{
	I dSum=0.0;

	for(typename GenericMassIntensMap<M,I>::const_iterator it=this->begin();it!=this->end();it++)
		dSum+=it->second*it->second;

	return dSum;
}

template<class M,class I>
I GenericMassIntensMap<M,I>::GetTotalIntensMassProd() const
{
	I dSum=0.0;

	for(typename GenericMassIntensMap<M,I>::const_iterator it=this->begin();it!=this->end();it++)
		dSum+=it->second*it->first;

	return dSum;
}

template<class M,class I>
I GenericMassIntensMap<M,I>::GetTotalIntensQuadMass() const
{
	I dSum=0.0;

	for(typename GenericMassIntensMap<M,I>::const_iterator it=this->begin();it!=this->end();it++)
		dSum+=it->second*it->second*it->first;

	return dSum;
}

template<class M,class I>
typename GenericMassIntensMap<M,I>::const_iterator GenericMassIntensMap<M,I>::GetBasePeak() const
{
	typename GenericMassIntensMap<M,I>::const_iterator itBasePeak=this->begin();

	for(typename GenericMassIntensMap<M,I>::const_iterator itPeak=this->begin();itPeak!=this->end();itPeak++)
		if(itBasePeak->second<itPeak->second)
			itBasePeak=itPeak;

	return itBasePeak;
}

template<class M,class I>
typename GenericMassIntensMap<M,I>::iterator GenericMassIntensMap<M,I>::GetBasePeak()
{
	typename GenericMassIntensMap<M,I>::iterator itBasePeak=this->begin();

	for(typename GenericMassIntensMap<M,I>::iterator itPeak=this->begin();itPeak!=this->end();itPeak++)
		if(itBasePeak->second<itPeak->second)
			itBasePeak=itPeak;

	return itBasePeak;
}

template<class M,class I>
typename GenericMassIntensMap<M,I>::const_iterator GenericMassIntensMap<M,I>::GetNearestPeak(M nMass) const
{
	typename GenericMassIntensMap<M,I>::const_iterator itUpperNeighb=lower_bound(nMass);
	typename GenericMassIntensMap<M,I>::const_iterator itLowerNeighb=itUpperNeighb;itLowerNeighb--;

	if(itUpperNeighb==this->end()) return itLowerNeighb;
	if(itLowerNeighb==this->end()) return itUpperNeighb;

	return nMass-itLowerNeighb->first<=itUpperNeighb->first-nMass?itLowerNeighb:itUpperNeighb;
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::Normalize(I dIntens,bool bMax)
//if bMax==false, the basepeak will have intensity 1, else the intensity sum
{
	I dIntensOld=bMax?GetMaxIntens():GetTotalIntens();
	if(!dIntensOld) return *this;

	return operator*=(dIntens/dIntensOld);
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::Round(int nDigits)
{
	bool bNullIntens=false;
	I dFactor=1;

	for(int iMult=0;iMult<nDigits;iMult++) dFactor*=10.0;

	for(typename GenericMassIntensMap<M,I>::iterator it=this->begin();it!=this->end();it++)
		if(!(it->second=floor(dFactor*it->second+.5)/dFactor)) 
			bNullIntens=true;

	if(bNullIntens) CutNoise(0.0);

	return *this;
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::CutNoise(I dCut)
{
	typename GenericMassIntensMap<M,I>::iterator itPeak=this->begin(),itDelete;

	while(itPeak!=this->end())
	{
		itDelete=itPeak++;
		if(itDelete->second<=dCut) erase(itDelete);
	}

	return *this;
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::Convolution(
	const std::map<M,I> &MIM1, const std::map<M,I> &MIM2)
{
	this->clear();
	for(typename GenericMassIntensMap<M,I>::const_iterator it1=MIM1.begin();it1!=MIM1.end();it1++)
		for(typename GenericMassIntensMap<M,I>::const_iterator it2=MIM2.begin();it2!=MIM2.end();it2++)
			this->operator[](it1->first+it2->first)+=it1->second*it2->second;

	return *this;
}

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::DelMasses(const std::map<M,I> &MIM)
{
	typename GenericMassIntensMap<M,I>::iterator itThis=this->begin();
	typename GenericMassIntensMap<M,I>::const_iterator itMIM=MIM.begin();

	while(itThis!=this->end()&&itMIM!=MIM.end())
	{
		if(itThis->first<itMIM->first)
			itThis++;
		else
		if(itThis->first>itMIM->first)
			itMIM++;
		else // itThis->first==itMIM->first
		{
			itMIM++;
			// erase(itThis++);
			// not sure if the operator++ still works after *itThis is erased;
			// better use the following code:
			typename GenericMassIntensMap<M,I>::iterator itDel=itThis++;
			this->erase(itDel);
		}
	}

	return *this;
} 

template<class M,class I>
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::ShiftMasses(const M &dShift)
{
	GenericMassIntensMap<M,I> MIM;
	typename GenericMassIntensMap<M,I>::iterator itMIM=MIM.begin();
	
	for(typename GenericMassIntensMap<M,I>::iterator it=this->begin();it!=this->end();it++)
		itMIM=MIM.insert(itMIM,std::pair<M,I>(it->first+dShift,it->second));

	return *this=MIM;
}

template<class M,class I>	
const GenericMassIntensMap<M,I>& GenericMassIntensMap<M,I>::Maximum(const std::map<M,I> &MIM)
{//not yet optimized
	for(typename GenericMassIntensMap<M,I>::const_iterator it=MIM.begin();it!=MIM.end();it++)
	{
		I& dIntens=operator[](it->first);
		dIntens=MAX(dIntens,it->second);
	}

	return *this;
}

template<class M,class I> 
GenericMassIntensMap<M,I> operator+(
	const GenericMassIntensMap<M,I> &MIM1,
	const GenericMassIntensMap<M,I> &MIM2)
{
	return GenericMassIntensMap<M,I>(MIM1)+=MIM2;
}

template<class M,class I> 
GenericMassIntensMap<M,I> operator*(
	const GenericMassIntensMap<M,I> &MIM1,
	const GenericMassIntensMap<M,I> &MIM2)
{
	return GenericMassIntensMap<M,I>().Convolution(MIM1,MIM2);
}

template<class M,class I> 
GenericMassIntensMap<M,I> operator*(
	const GenericMassIntensMap<M,I> &MIM1,I dFactor)
{
	return GenericMassIntensMap<M,I>(MIM1)*=dFactor;
}

template<class M,class I>
void GenericMassIntensMap<M,I>::AddPeak(M nMass,I dIntens)
{
	typename GenericMassIntensMap<M,I>::iterator it=this->empty()?this->begin():this->end()--;

	this->insert(it,std::pair<M,I>(nMass,dIntens));
}

template<class M,class I>
I GenericMassIntensMap<M,I>::GetIntens(M nMass) const
{
	typename GenericMassIntensMap<M,I>::const_iterator it=this->find(nMass);
	return it==this->end()?0.0:it->second;
}

template<class M,class I>
M	GenericMassIntensMap<M,I>::GetMinMass() const
{
	return (M) (this->empty()?0.0:this->begin()->first);
}

template<class M,class I>
M	GenericMassIntensMap<M,I>::GetMaxMass() const
{
	return (M) (this->empty()?0.0:this->rbegin()->first);
}

#endif
