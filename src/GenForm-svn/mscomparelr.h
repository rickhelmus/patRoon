/*
 *  mscomparelr.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef	MS_COMPARE_LR_H
#define MS_COMPARE_LR_H

#include "mslrmassintensmap.h"

typedef double (*CompareMassIntensMap)(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical);

double CompareEuklidDistance(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical);

double CompareNormalizedSumSquares(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical);

double CompareNormalizedDotProduct(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical);

double CompareShimadzuSimilarityIndex(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical);

double CompareLeastSquares(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical);

double SimilarityWeightedCosine(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical,
	double dWeightMass = 1.0,
	double dWeightIntens = 0.5);

double SimilarityComposite(
	const LrMassIntensMap& mapObserved,
	const LrMassIntensMap& mapTheoretical,
	double dWeightMass = 1.0,
	double dWeightIntens = 0.5);


#endif
