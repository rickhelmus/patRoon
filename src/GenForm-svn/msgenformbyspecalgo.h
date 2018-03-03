/*
 *  msgenformbyspecalgo.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_GEN_FORM_BY_SPEC_ALGO_H
#define MS_GEN_FORM_BY_SPEC_ALGO_H

#include "msgenformbymassalgo.h"
#include "mslrmassintensmap.h"

class GenFormBySpecAlgo : public GenFormByMassAlgo<int>
{
private:
	
	std::vector<bool>			m_vecNullIntens;
								
public:							
								GenFormBySpecAlgo(){}
								GenFormBySpecAlgo(
									const ElementIntervalMap& EIM,
									const LrMassIntensMap& MIM);
								~GenFormBySpecAlgo(){}

	void						Init(
									const ElementIntervalMap& EIM,
									const LrMassIntensMap& MIM);

	void						Init(
									const ElementIntervalMap& EIM,
									int iMinMass, int iMaxMass);

	void						begin();
	void						operator++();
};

#endif
