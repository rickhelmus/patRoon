/*
 *  mshrmassintensmap.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_HR_MASS_INTENS_MAP_H
#define MS_HR_MASS_INTENS_MAP_H

#include <list>
#include <iostream>
#include "msgenericmassintensmap.h"

class HrMassIntensMap : public GenericMassIntensMap<double,double>
{
public:
							HrMassIntensMap();
							HrMassIntensMap(const GenericMassIntensMap<double,double> &MIM);

	const HrMassIntensMap&	operator=(const std::map<double,double> &MIM);
	const HrMassIntensMap&	operator=(const std::map<int,double> &MIM);

	bool					IsHighRes()										const;
							// TRUE, if at least one mass is not integer

	double					GetSumIntens(double dMinMass,double dMaxMass)	const;
	double					GetSumIntens(double dMass,int iPPM)				const;
							//iPPM: Precision given in parts per million

	double					GetMaxIntens(double dMinMass,double dMaxMass)	const;
	double					GetMaxIntens(double dMass, int iPPM)			const;

	const_iterator			GetMaxPeak(double dMinMass,double dMaxMass)		const;
	const_iterator			GetMaxPeak(double dMass, int iPPM)				const;

	double					GetNearestIntens(double dMass, double dMaxDiff)	const;
	double					GetNearestIntens(double dMass, int iPPM)		const;

	const_iterator			GetNearestPeak(double dMass)					const;
};
	
std::ostream& operator<<(std::ostream &out, const HrMassIntensMap& MIM);
std::istream& operator>>(std::istream &in, HrMassIntensMap& MIM);

void ComputeDecomposition(
	const HrMassIntensMap& MIM,
	std::list<HrMassIntensMap>& listMIM,
	double dDiversity=2.0);
	
int GetMaxClusterSize(const HrMassIntensMap& MIM);
// size in number of peaks

HrMassIntensMap& ComputeNeutralLossSpec(
	const HrMassIntensMap& MIMSource, 
	HrMassIntensMap& MIMTarget, 
	double sMolecularMass=0.0);


double CalcDiffPPM(double dCalc, double dExp);

double CalcPPM(double dCalc, double dDiff);

#endif
