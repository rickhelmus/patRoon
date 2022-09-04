/*
 *  msgenformmatchisotopemsms.cpp, part of GenForm by M. Meringer Copyright (C) 2015
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>

#include "msdef.h"
#include "msstringconversion.h"
#include "mshrmassintensmaputil.h"
#include "mslrmassintensmap.h"
#include "msconvert.h"
#include "mselementmultmaputil.h"
#include "mselementmultmapfilter.h"
#include "mselementintervalmaputil.h"
#include "msexception.h"
#include "msaddion.h"
#include "msgenformbymassalgo.h"
#include "msisotope.h"
#include "mscomparelr.h"
#include "mshrmatchformula.h"

using namespace std;

void GenFormMatchIsotopeMsMsReference()
{
	cerr << "Formula calculation from MS and MS/MS data as described in\n"
		 <<	"Meringer et al (2011) MATCH Commun Math Comput Chem 65: 259-290\n";
//		 << "http://www.pmf.kg.ac.rs/match/electronic_versions/match65/n2/match65%282%29_259-290.pdf\n";
}

void GenFormMatchIsotopeMsMsUsage(const string& strProgName)
{
	GenFormMatchIsotopeMsMsReference();

	cerr << "Usage: " << strProgName << " ms=<filename> [msms=<filename>] [out=<filename>]\n";
	cerr << "\t[exist[=mv]] [m=<number>] [ion=-e|+e|-H|+H|+Na] [cha=<number>]\n";
	cerr << "\t[ppm=<number>] [msmv=ndp|nsse|nsae] [acc=<number>] [rej=<number>]\n";
	cerr << "\t[thms=<number>] [thmsms=<number>] [thcomb=<number>]\n";
	cerr << "\t[sort[=ppm|msmv|msmsmv|combmv]] [el=<elements> [oc]] [ff=<fuzzy formula>]\n";
	cerr << "\t[vsp[=<even|odd>]] [vsm2mv[=<value>]] [vsm2ap2[=<value>]] [hcf] [kfer[=ex]]\n";
	cerr << "\t[wm[=lin|sqrt|log]] [wi[=lin|sqrt|log]] [exp=<number>] [oei]\n";
	cerr << "\t[dbeexc=<number>] [ivsm2mv=<number>] [vsm2ap2=<number>]\n";
	cerr << "\t[oms[=<filename>]] [omsms[=<filename>]] [oclean[=<filename>]]\n";
	cerr << "\t[analyze [loss] [intens]] [dbe] [cm] [pc] [sc] [max=<number>]\n";
	cerr << "Explanation:\n";
//5.1
	cerr << "\tms\t: filename of MS data (*.txt)\n";
	cerr << "\tmsms\t: filename of MS/MS data (*.txt)\n";
//5.2
	cerr << "\tout\t: output generated formulas\n";
	cerr << "\texist\t: allow only molecular formulas for that at least one\n\t\t  structural formula exists;overrides vsp, vsm2mv, vsm2ap2;\n\t\t  argument mv enables multiple valencies for P and S\n";
	cerr << "\tm\t: experimental molecular mass (default: mass of MS basepeak)\n";
	cerr << "\tion\t: type of ion measured (default: M+H)\n";
	cerr << "\tppm\t: accuracy of measurement in parts per million (default: 5)\n";
	cerr << "\tmsmv\t: MS match value based on normalized dot product, normalized\n\t\t  sum of squared or absolute errors (default: nsae)\n";
	cerr << "\tacc\t: allowed deviation for full acceptance of MS/MS peak in ppm\n\t\t  (default: 2)\n";
	cerr << "\trej\t: allowed deviation for total rejection of MS/MS peak in ppm\n\t\t  (default: 4)\n";
	cerr << "\tthms\t: threshold for the MS match value\n";
	cerr << "\tthmsms\t: threshold for the MS/MS match value\n";
	cerr << "\tthcomb\t: threshold for the combined match value\n";
	cerr << "\tsort\t: sort generated formulas according to mass deviation in ppm,\n\t\t  MS match value, MS/MS match value or combined match value\n";
//6.1
	cerr << "\tel\t: used chemical elements (default: CHBrClFINOPSSi)\n";
	cerr << "\toc\t: only organic compounds, i.e. with at least one C atom\n";
	cerr << "\tff\t: overwrites el and oc and uses fuzzy formula for limits of\n\t\t  element multiplicities\n";
	cerr << "\thet\t: formulas must have at least one hetero atom\n";
//6.2
	cerr << "\tvsp\t: valency sum parity (even for graphical formulas)\n";
	cerr << "\tvsm2mv\t: lower bound for valency sum - 2 * maximum valency\n\t\t  (>=0 for graphical formulas)\n";
	cerr << "\tvsm2ap2\t: lower bound for valency sum - 2 * number of atoms + 2\n\t\t  (>=0 for graphical connected formulas)\n";
//6.3
	cerr << "\thcf\t: apply Heuerding-Clerc filter\n";
	cerr << "\tkfer\t: apply Kind-Fiehn element ratio (extended) ranges\n";
//6.4
	cerr << "\twm\t: m/z weighting for MS/MS match value\n";
	cerr << "\twi\t: intensity weighting for MS/MS match value\n";
	cerr << "\texp\t: exponent used, when wi is set to log\n";
	cerr << "\toei\t: allow odd electron ions for explaining MS/MS peaks\n";
	cerr << "\tdbeexc\t: excess of double bond equivalent for ions\n";
	cerr << "\tivsm2mv\t: lower bound for valency sum - 2 * maximum valency\n\t\t  for fragment ions\n";
	cerr << "\tivsm2ap2: lower bound for valency sum - 2 * number of atoms + 2\n\t\t  for fragment ions\n";
//6.5
	cerr << "\toms\t: write scaled MS peaks to output\n";
	cerr << "\tomsms\t: write weighted MS/MS peaks to output\n";
	cerr << "\toclean\t: write explained MS/MS peaks to output\n";
	cerr << "\tanalyze\t: write explanations for MS/MS peaks to output\n";
	cerr << "\tloss\t: for analyzing MS/MS peaks write losses instead of fragments\n";
	cerr << "\tintens\t: write intensities of MS/MS peaks to output\n";
	cerr << "\tdbe\t: write double bond equivalents to output\n";
	cerr << "\tcm\t: write calculated ion masses to output\n";
	cerr << "\tpc\t: output match values in percent\n";
//6.6
	cerr << "\tsc\t: strip calculated isotope distributions\n";
//new
	cerr << "\tnoref\t: hide the reference information\n";
	cerr << "\tmax\t: maximum number of final candidates (0 is no limit)\n";
}

int GenFormMatchIsotopeMsMs(const string& strProgName,map<string,string>& mapArgValue)
{
	map<string,string>::iterator it;
	string strEl,strEIM,strMs,strIon="M+H",strMsMs,strMsMv="nsae";
	bool3 b3EvenValSum=BOOL3_PERH;
	int iMinDiffMaxVal=numeric_limits<int>::min();
	int iMinDiffAtomCount=numeric_limits<int>::min();
	int iDbeExcess=3;
	int iFragMinDiffMaxVal=0;
	int iFragMinDiffAtomCount=0;
	bool bHCFilter=false,bAnalyze=false,bLoss=false;
	bool bKFElementRatios=false, bKFElementRatiosExtendedRange=false;
	bool bAllowRadicalIons=false,bTex=false;
	bool bHighMass=true,bStripCalc=false;
	bool bPercent=false,bWriteCalcMass=false,bWriteDBE=false;
	bool bObligatoryHeteroAtom=false,bMultipleValencies=false;
	bool bOrganicCompound=false;
	bool bIntens=false, bShowRef=true;
	double dMass=0.0,dPPM=5.0;
	double dMinMsMv=-numeric_limits<double>::max();
	double dMinMsMsMv=-numeric_limits<double>::max();
	double dMinCombMv=-numeric_limits<double>::max();
	int iPpmAccept=2,iPpmReject=4,iExponent=5;
	string strWeightMass="false",strWeightIntens="false";
	ostream *pOst=NULL,*pOutMs=NULL,*pOutMsMs=NULL,*pOutCleanMsMs=NULL;

	enum SortMethod { SortPpm, SortMsMv, SortMsMsMv, SortCombMv, SortUndefined };
	const char* pSortValue[SortUndefined]={"ppm","msmv","msmsmv","combmv"};
	unsigned int iSortMethod=SortUndefined;
	unsigned int iMaxFinal = 0; // Added by Rick Helmus
	 
	if((it=mapArgValue.find("m"))!=mapArgValue.end())
	{ dMass=atof(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("ppm"))!=mapArgValue.end())
	{ dPPM=atof(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("acc"))!=mapArgValue.end())
	{ iPpmAccept=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("rej"))!=mapArgValue.end())
	{ iPpmReject=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("ms"))!=mapArgValue.end())
	{ strMs=it->second; mapArgValue.erase(it); }
	else
	{ cerr << "Error: missing argument ms\n"; GenFormMatchIsotopeMsMsUsage(strProgName); exit(1); }

	if((it=mapArgValue.find("thms"))!=mapArgValue.end())
	{ dMinMsMv=atof(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("thmsms"))!=mapArgValue.end())
	{ dMinMsMsMv=atof(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("thcomb"))!=mapArgValue.end())
	{ dMinCombMv=atof(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("el"))!=mapArgValue.end())
	{ strEl=it->second; mapArgValue.erase(it); }

	if((it=mapArgValue.find("oc"))!=mapArgValue.end())
	{ bOrganicCompound=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("ff"))!=mapArgValue.end())
	{ strEIM=it->second; mapArgValue.erase(it); }

	if((it=mapArgValue.find("het"))!=mapArgValue.end())
	{ bObligatoryHeteroAtom=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("vsp"))!=mapArgValue.end())
	{ b3EvenValSum=it->second!="odd"?BOOL3_TRUE:BOOL3_FALSE; mapArgValue.erase(it); }

	if((it=mapArgValue.find("vsm2mv"))!=mapArgValue.end())
	{ iMinDiffMaxVal=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("vsm2ap2"))!=mapArgValue.end())
	{ iMinDiffAtomCount=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("exist"))!=mapArgValue.end())
	{ b3EvenValSum=BOOL3_TRUE;iMinDiffMaxVal=iMinDiffAtomCount=0;
	  bMultipleValencies=(it->second=="mv")?true:false; mapArgValue.erase(it); }

	if((it=mapArgValue.find("out"))!=mapArgValue.end())
	{ pOst=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("oms"))!=mapArgValue.end())
	{ pOutMs=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("omsms"))!=mapArgValue.end())
	{ pOutMsMs=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("oclean"))!=mapArgValue.end())
	{ pOutCleanMsMs=it->second.empty()?&cout:new ofstream(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("hcf"))!=mapArgValue.end())
	{ bHCFilter=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("kfer"))!=mapArgValue.end())
	{ bKFElementRatios=true;
	  bKFElementRatiosExtendedRange=(it->second=="ex")?true:false; mapArgValue.erase(it); }

	if((it=mapArgValue.find("oei"))!=mapArgValue.end())
	{ bAllowRadicalIons=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("analyze"))!=mapArgValue.end())
	{ bAnalyze=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("loss"))!=mapArgValue.end())
	{ bLoss=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("sc"))!=mapArgValue.end())
	{ bStripCalc=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("wi"))!=mapArgValue.end())
	{ strWeightIntens=it->second; mapArgValue.erase(it); }

	if((it=mapArgValue.find("wm"))!=mapArgValue.end())
	{ strWeightMass=it->second; mapArgValue.erase(it); }

	if((it=mapArgValue.find("ion"))!=mapArgValue.end())
	{ strIon=it->second; mapArgValue.erase(it); }

	if((it=mapArgValue.find("exp"))!=mapArgValue.end())
	{ iExponent=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("dbeexc"))!=mapArgValue.end())
	{ iDbeExcess=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("ivsm2mv"))!=mapArgValue.end())
	{ iFragMinDiffMaxVal=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("ivsm2ap2"))!=mapArgValue.end())
	{ iFragMinDiffAtomCount=atoi(it->second.c_str()); mapArgValue.erase(it); }

	if((it=mapArgValue.find("msmv"))!=mapArgValue.end())
	{ strMsMv=it->second; mapArgValue.erase(it); }

	if((it=mapArgValue.find("msms"))!=mapArgValue.end())
	{ strMsMs=it->second; mapArgValue.erase(it); }

	if((it=mapArgValue.find("tex"))!=mapArgValue.end())
	{ bTex=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("pc"))!=mapArgValue.end())
	{ bPercent=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("dbe"))!=mapArgValue.end())
	{ bWriteDBE=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("cm"))!=mapArgValue.end())
	{ bWriteCalcMass=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("intens"))!=mapArgValue.end())
	{ bIntens=true; mapArgValue.erase(it); }

	if((it=mapArgValue.find("noref"))!=mapArgValue.end())
	{ bShowRef=false; mapArgValue.erase(it); }

	if((it=mapArgValue.find("sort"))!=mapArgValue.end())
	{
		if(it->second.empty()) iSortMethod=SortCombMv;
		else
		for(iSortMethod=0;iSortMethod!=SortUndefined;iSortMethod++)
			if(it->second==pSortValue[iSortMethod]) break;
		
		if(iSortMethod==SortUndefined)
		{
			cerr << "Error:\tinvalid value for key 'sort'\n\t " << it->second;
			cerr << "\n"; GenFormMatchIsotopeMsMsUsage(strProgName); exit(1);
		}

		mapArgValue.erase(it); 
	}
	
	// Added by Rick Helmus
	if((it=mapArgValue.find("max"))!=mapArgValue.end())
	{ iMaxFinal=atoi(it->second.c_str()); mapArgValue.erase(it); }
	
	if(!mapArgValue.empty())
	{
		cerr << "Error:\tunknown key(s)\n\t ";
		for(it=mapArgValue.begin();it!=mapArgValue.end();it++)
			cerr << it->first;
		cerr << "\n"; GenFormMatchIsotopeMsMsUsage(strProgName); exit(1);
	}

	if(bShowRef)
		GenFormMatchIsotopeMsMsReference();

	LrMassIntensMap lrmCalc,lrmExp;
	HrMassIntensMap hrmMs,hrmMsMs,hrmMsMsNorm,hrmMsMsClean;

	ifstream in(strMs.c_str());
	if(!in)
	{ cerr << "Error opening " << strMs << "\n"; exit(1); }
	while(in)
	{
		double dMass,dIntens;
		in >> dMass >> dIntens;

		if(in)
			hrmMs[dMass]=dIntens;
	}

	Init(lrmExp,hrmMs);

	cout << "finished reading MS file (" << lrmExp.size() << " peaks)\n";

	//Normalize to intensity sum = 1; 
	//this is important for the ShimadzuSimilarityIndex
	lrmExp.Normalize(1.0,false);
	hrmMs.Normalize(1.0,false);

	if(pOutMs) *pOutMs << hrmMs;

	if(dMass==0.0)
	{  
		if(hrmMs.GetBasePeak()!=hrmMs.end())
			dMass=hrmMs.GetBasePeak()->first;
		else
		{ cerr << "Error searching ms base peak\n"; exit(1); }
	}

	if(!strMsMs.empty())
	{
		ifstream inMsMs(strMsMs.c_str());
		if(!inMsMs)
		{ cerr << "Error opening " << strMsMs << "\n"; exit(1); }
		while(inMsMs)
		{
			double dMass,dIntens;
			inMsMs >> dMass >> dIntens;

			if(inMsMs)
				hrmMsMs[dMass]=dIntens;
		}

		cout << "finished reading MS/MS file (" << hrmMsMs.size() << " peaks)\n" << flush;
	}

	hrmMsMs.Normalize();
	hrmMsMsNorm=hrmMsMs;
	WeightSpectrum(hrmMsMs,strWeightMass,strWeightIntens,iExponent);
	hrmMsMs.Normalize();
	if(pOutMsMs) *pOutMsMs << hrmMsMs;

	AddIon AI;
	try{AI=AddIon(strIon);}
	catch(MsError& ME)
	{
		cerr << ME.what() << "\n";
		return 1;
	}
	dMass=AI.CalcMolMass(dMass);

	double dMinMass=dMass*(1.0-dPPM/1000000.0),dMaxMass=dMass*(1.0+dPPM/1000000.0);
	unsigned long iTotalCount=0,iValidCount=0,iValidElementRatioCount=0;
	unsigned long iIsotopeCount=0,iMsMsCount=0,iCombCount=0;
	ElementIntervalMap EIM;

	if(!strEIM.empty())
	{
		try{Init(EIM,strEIM);}
		catch(MsError& ME)
		{
			cerr << ME.what() << "\n";
			return 1;
		}
	}
	else
	{
		Interval I(0,numeric_limits<ElementMult>::max());
		if(strEl.empty()) strEl="CHBrClFINOPSSi";

		try{Init(EIM,strEl);}
		catch(MsError& ME)
		{
			cerr << ME.what() << "\n";
			return 1;
		}

		EIM.SetInterval(I);

		if(bOrganicCompound)
			EIM.SetInterval(e_C,Interval(1,numeric_limits<ElementMult>::max()));
	}

	// raise lower limits
	if(!EIM.RaiseLowerLimits(AI.GetFormulaMinus()))
	{
		cerr << "Invalid  element bounds when considering ionization type: " << EIM << endl;
		return 1;
	}

	// in order to write explained peaks and explaining formula to output
	map<double,list<ElementMultMap> > mapMassEMM;

	// in order to write candidate formulas sorted by MV
	multimap<double,string> mapSort;

//	CpuTime	time; time.Start();

	GenFormByMassAlgo<double> BFBMA(EIM,dMinMass,dMaxMass);

	// Modified by Rick Helmus
	//for(BFBMA.begin();!BFBMA.end();BFBMA.operator++())
	for(BFBMA.begin();!BFBMA.end()&&(iMaxFinal==0 || iCombCount<iMaxFinal);BFBMA.operator++())
	{
		ElementMultMap emmIon,emmMol;

		if(bObligatoryHeteroAtom)
		{
			BFBMA.GetElementMultMap(emmMol);
			if(!HasHeteroAtom(emmMol))
				continue;
		}

		iTotalCount++;

		if(emmMol.empty()&&bMultipleValencies)
			BFBMA.GetElementMultMap(emmMol);

		bool bValid=bMultipleValencies?
			IsGraphicalMultVal(emmMol):
			BFBMA.IsValid(b3EvenValSum,iMinDiffMaxVal,iMinDiffAtomCount);

		if(bValid) 
		{
			iValidCount++;

			if(emmMol.empty())
				BFBMA.GetElementMultMap(emmMol);

			if(bHCFilter)
			{ if(HeuerdingClercFilter(emmMol)==false) continue; }

			if(bKFElementRatios)
			{ if(KindFiehnElementRatios(emmMol,bKFElementRatiosExtendedRange)==false) continue; }

			iValidElementRatioCount++;

			AI.CalcEMM(emmMol,emmIon);
			Init(lrmCalc,emmIon);

			if(bHighMass)
				Init(lrmExp,hrmMs,GetNominalMass<int>(emmIon)-GetNominalMass<double>(emmIon)+0.5);

//cout << setprecision(12) << GetNominalMass<double>(emmIon)-(c_dMassElectron*iCharge) << "\n";

			if(bStripCalc)
			{
				//remove Peaks in calculated spectrum that are not present in Exp spectrum
				lrmCalc.erase(lrmCalc.begin(),lrmCalc.lower_bound(lrmExp.begin()->first));
				lrmCalc.erase(lrmCalc.upper_bound(lrmExp.rbegin()->first),lrmCalc.end());
			}

//cout << lrmCalc << '\n' << lrmExp;

			double dMs=1.0,dMsMs=0.0;

			if(strMsMv=="ndp")
				dMs-=CompareNormalizedDotProduct(lrmExp,lrmCalc);
			else
			if(strMsMv=="nsse")
				dMs-=CompareNormalizedSumSquares(lrmExp,lrmCalc);
			else
			{
				//LrMassIntensMap hrmMsLr;

				//Init(hrmMsLr,ShiftMass(hrmMs,HrMassIntensMap(),GetNominalMass<int>(emmIon)-GetNominalMass<double>(emmIon)));
				dMs-=CompareShimadzuSimilarityIndex(lrmExp,lrmCalc);
								
			//	dMs-=CompareShimadzuSimilarityIndexSpecial(lrmExp,lrmCalc);
			}
			
			if(dMs<dMinMsMv) continue;
			iIsotopeCount++;

			if(hrmMsMs.size())
				dMsMs=ComputeHrMatchvalue(hrmMsMs,emmIon,iPpmAccept,iPpmReject,
					iDbeExcess,iFragMinDiffMaxVal,iFragMinDiffAtomCount,
					bAllowRadicalIons,AI.GetCharge()>0,NULL,(bAnalyze||pOutCleanMsMs)?&mapMassEMM:NULL);

			if(dMsMs<dMinMsMsMv) continue;
			iMsMsCount++;

			double dMvComb=dMsMs*dMs;

			if(dMsMs*dMs<dMinCombMv) continue;
			iCombCount++;
			
			if(bPercent)
			{
				dMs*=100.0;dMsMs*=100.0;dMvComb*=100.0;
			}

			if(pOst)
			{
//	Use this code, if all three MS matchvalues shall be printed
//				pOst->setf(ios::fixed,ios::floatfield);
//				double dMassCalc=GetNominalMass<double>(emmMol);
//				*pOst << setw(15) << setiosflags(ios_base::left) << (bTex?TexName(emmMol):ChemName(emmMol)) << (bTex?"&\t":"\t")
//					  << setiosflags(ios_base::right) << setw(6) << setprecision(1) << CalcDiffPPM(dMassCalc,dMass)
//					  << (bTex?"&\t":"\t") << setprecision(6) << 1.0-CompareNormalizedDotProduct(lrmExp,lrmCalc)
//					  << (bTex?"&\t":"\t") << setprecision(6) << 1.0-CompareShimadzuSimilarityIndex(lrmExp,lrmCalc)
//					  << (bTex?"&\t":"\t") << setprecision(6) << 1.0-CompareNormalizedSumSquares(lrmExp,lrmCalc) << "\n" << resetiosflags(ios::adjustfield);

				int iWritePrec=bPercent?3:5;
				stringstream sStream;

				pOst->setf(ios::fixed,ios::floatfield);
				sStream.setf(ios::fixed,ios::floatfield);
				double dMassCalc=GetNominalMass<double>(emmMol);
				double dPpmCalc=CalcDiffPPM(dMassCalc,dMass);
				sStream << setw(15) << setiosflags(ios_base::left) << (bTex?TexName(emmMol):ChemName(emmMol)) << (bTex?"&\t":"\t");
				if(bWriteDBE)
					sStream << setiosflags(ios_base::right) << setw(5) << setprecision(1) << GetDBE(emmMol) << (bTex?"&\t":"\t");
				if(bWriteCalcMass)
					sStream << setiosflags(ios_base::right) << setw(10) << setprecision(5) << AI.CalcIonMass(dMassCalc) << (bTex?"&\t":"\t");
				sStream << setiosflags(ios_base::right) << setw(6) << setprecision(1) << dPpmCalc
					  << (bTex?"&\t":"\t") << setiosflags(ios_base::right) << setw(7) << setprecision(iWritePrec) << dMs;
				if(hrmMsMs.size())
					sStream << (bTex?"&\t":"\t") << setiosflags(ios_base::right) << setw(7) << setprecision(iWritePrec) << dMsMs
						  << (bTex?"&\t":"\t") << setiosflags(ios_base::right) << setw(7) << setprecision(iWritePrec) << dMvComb;
				sStream << (bTex?"\\\\ \\hline\n":"\n") << resetiosflags(ios::adjustfield);
				*pOst << sStream.str();

				if(iSortMethod!=SortUndefined)
					switch(iSortMethod)
					{
						case SortPpm   : mapSort.insert(make_pair(-fabs(dPpmCalc),sStream.str())); break;
						case SortMsMv  : mapSort.insert(make_pair(dMs,sStream.str())); break;
						case SortMsMsMv: mapSort.insert(make_pair(dMsMs,sStream.str())); break;
						case SortCombMv: mapSort.insert(make_pair(dMvComb,sStream.str())); break;
					}

				if(bAnalyze)
					for(map<double,list<ElementMultMap> >::const_iterator itMass=mapMassEMM.begin();itMass!=mapMassEMM.end();itMass++)
					{
						const list<ElementMultMap>& listEMM=itMass->second;
						*pOst << setw(10) << setiosflags(ios::fixed) << setprecision(5) << itMass->first << '\t';
						if(bIntens)
							*pOst << setw(10) << setiosflags(ios::fixed) << setprecision(8) << hrmMsMsNorm.GetIntens(itMass->first) << '\t';
						for(list<ElementMultMap>::const_iterator itEMM=listEMM.begin();itEMM!=listEMM.end();itEMM++)
						{
							double dIonMassCalc=GetNominalMass<double>(*itEMM)-(c_dMassElectron*AI.GetCharge());
							ElementMultMap emmLoss;
							string strLoss;
							if(bLoss)
							{
								emmLoss=emmIon-*itEMM;
								strLoss=emmLoss.empty()?"no loss":"-"+ChemName(emmLoss);
							}
							if(itEMM!=listEMM.begin()) 
								*pOst << "          \t" << (bIntens?"          \t":"");
							*pOst << setw(15) << setiosflags(ios_base::left) << (bLoss?strLoss:ChemName(*itEMM))  << '\t';
							if(bWriteDBE)
								*pOst << setiosflags(ios_base::right) << setw(5) << setprecision(1) << CalcDBE(*itEMM) << (bTex?"&\t":"\t");
							if(bWriteCalcMass) 
								*pOst << setw(10) << setiosflags(ios::fixed) << setprecision(5) << dIonMassCalc << '\t';
							*pOst << setw(6) << setprecision(1) << CalcDiffPPM(dIonMassCalc,itMass->first) << '\n' << resetiosflags(ios::adjustfield);
						}
					}

				if(pOutCleanMsMs)
					for(map<double,list<ElementMultMap> >::const_iterator itMass=mapMassEMM.begin();itMass!=mapMassEMM.end();itMass++)
						hrmMsMsClean.AddPeak(itMass->first,hrmMsMsNorm.GetIntens(itMass->first));
			}
		}
	}

	if(pOutCleanMsMs) *pOutCleanMsMs << hrmMsMsClean;

	if(bHCFilter||bKFElementRatios)
		cout << iCombCount << "/" << iMsMsCount << "/" << iIsotopeCount << "/" << iValidElementRatioCount << "/" << iValidCount << "/" << iTotalCount
		      << " (final/MSMS-/MS-/element-ratio-filter/valid/total) formula(s)\n";// in " << setprecision(1) << time.GetSecs() << "s\n";
	else
		cout << iCombCount << "/" << iMsMsCount << "/" << iIsotopeCount << "/" <<  iValidCount << "/" << iTotalCount 
		      << " (final/MSMS-/MS-filter/valid/total) formula(s)\n";// in " << setiosflags(ios::fixed) << setprecision(1) << time.GetSecs() << "s\n";

	if(iSortMethod!=SortUndefined)
		for(multimap<double,string>::reverse_iterator itSort=mapSort.rbegin();itSort!=mapSort.rend();itSort++)
			*pOst << itSort->second;

	if(pOst!=NULL&&pOst!=&cout) delete pOst;
	if(pOutMs!=NULL&&pOutMs!=&cout) delete pOutMs;
	if(pOutMsMs!=NULL&&pOutMsMs!=&cout) delete pOutMsMs;
	if(pOutCleanMsMs!=NULL&&pOutCleanMsMs!=&cout) delete pOutCleanMsMs;

	return 0;

}
