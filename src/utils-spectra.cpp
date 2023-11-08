#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h>

#include <Rcpp.h>

#include "utils.h"
#include "utils-spectra.h"

namespace {

void specApply(Rcpp::List spectra, double rtMin, double rtMax, double mzMin, double mzMax,
               std::function<void(double, const Rcpp::NumericVector &, const Rcpp::NumericVector &)> func)
{
    const Rcpp::DataFrame specHeader = Rcpp::as<Rcpp::DataFrame>(spectra["header"]);
    const Rcpp::NumericVector hdRetTimes = specHeader["retentionTime"];
    const Rcpp::NumericVector hdMSLevels = specHeader["msLevel"];
    const Rcpp::NumericVector hdSeqNums = specHeader["seqNum"];
    const Rcpp::List specList = spectra["spectra"];
    const int specCount = hdRetTimes.size();
    const double rttol = 1E-4;
    
    // NOTE: assume header is RT ordered (low to high)
    for (int spi=0; spi<specCount; ++spi)
    {
        const double specRt = hdRetTimes[spi];

        if (hdMSLevels[spi] == 1 && numberGTE(specRt, rtMin, rttol))
        {
            if (numberLTE(specRt, rtMax, rttol))
            {
                const Rcpp::NumericMatrix peaklist = Rcpp::as<Rcpp::NumericMatrix>(specList[hdSeqNums[spi] - 1]); // -1: R 1 based index
                func(specRt, peaklist(Rcpp::_, 0), peaklist(Rcpp::_, 1));
            }
            else // surpassed RT range, abort
                return;
        }
    }
}
    
double getTotMZIntFromSpec(const Rcpp::NumericVector &peakMZs, const Rcpp::NumericVector &peakInts,
                           double mzMin, double mzMax)
{
    // NOTE: assume spectrum is mz ordered (low to high)
    
    const double mztol = 1E-8;
    const int pCount = peakMZs.length();
    double totInt = 0;
    for (int pi=0; pi<pCount; ++pi)
    {
        if (numberGTE(peakMZs[pi], mzMin, mztol))
        {
            if (numberLTE(peakMZs[pi], mzMax, mztol))
                totInt += peakInts[pi];
            else
                break; // surpassed mz range, abort
        }
    }
    
    return totInt;
}

struct BinnedSpectrum
{
    std::vector<double> mzs, intsLeft, intsRight;
    std::vector<int> IDsLeft, IDsRight;
};
    

}


// [[Rcpp::export]]
Rcpp::NumericVector loadEICIntensities(Rcpp::List spectra, Rcpp::DataFrame featList, Rcpp::NumericVector rtWindow)
{
    const Rcpp::NumericVector rets = featList["ret"];
    const Rcpp::NumericVector mzmins = featList["mzmin"];
    const Rcpp::NumericVector mzmaxs = featList["mzmax"];
    
    const double rtWin = Rcpp::as<double>(rtWindow);
    const int ftCount = rets.length();
    
    Rcpp::NumericVector intensities(ftCount);
    for (int fti=0; fti<ftCount; ++fti)
    {
        const double featRt = rets[fti];
        const double mzMin = mzmins[fti], mzMax = mzmaxs[fti];
        double closestRTDiff = -1, closestInt = 0;
        
        specApply(spectra, featRt - rtWin, featRt + rtWin, mzMin, mzMax,
                  [&](double specRt, const Rcpp::NumericVector &peakMZs, const Rcpp::NumericVector &peakInts)
                  {
                      const double rtdiff = fabs(featRt - specRt);
                      if (closestRTDiff == -1)
                          closestRTDiff = rtdiff;
                      else if (closestRTDiff < rtdiff)
                          return;
                      
                      closestRTDiff = rtdiff;
                      closestInt = getTotMZIntFromSpec(peakMZs, peakInts, mzMin, mzMax);
                  });
        
        intensities[fti] = closestInt;
    }
    
    return intensities;
}

// [[Rcpp::export]]
Rcpp::List loadEICs(Rcpp::List spectra, Rcpp::NumericVector rtMins, Rcpp::NumericVector rtMaxs,
                    Rcpp::NumericVector mzMins, Rcpp::NumericVector mzMaxs)
{
    const int EICCount = rtMins.length();
    Rcpp::List ret(EICCount);
    for (int eici=0; eici<EICCount; ++eici)
    {
        const double mzMin = mzMins[eici], mzMax = mzMaxs[eici];
        std::vector<double> EICTimes;
        std::vector<double> EICIntensities;
        
        specApply(spectra,  rtMins[eici], rtMaxs[eici], mzMin, mzMax,
                  [&](double specRt, const Rcpp::NumericVector &peakMZs, const Rcpp::NumericVector &peakInts)
                  {
                      EICTimes.push_back(specRt);
                      EICIntensities.push_back(getTotMZIntFromSpec(peakMZs, peakInts, mzMin, mzMax));
                  });

        if (EICTimes.size() > 2)
        {
            // compress data by removing any zero intensity datapoints that are inbetween two other zero intensity points.
            for (size_t ind=0; ind<(EICTimes.size()-2); )
            {
                if (EICIntensities[ind+2] != 0)
                    ind += 3;
                else if (EICIntensities[ind+1] != 0)
                    ind += 2;
                else if (EICIntensities[ind] != 0)
                    ++ind;
                else // all zero
                {
                    EICTimes.erase(EICTimes.begin() + (ind + 1));
                    EICIntensities.erase(EICIntensities.begin() + (ind + 1));
                }
            }
        }
        
        ret[eici] = Rcpp::DataFrame::create(Rcpp::Named("time") = EICTimes,
                                            Rcpp::Named("intensity") = EICIntensities);
    }
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::List makeSAFDInput(Rcpp::List spectra, Rcpp::NumericVector mzRange)
{
    const Rcpp::DataFrame specHeader = Rcpp::as<Rcpp::DataFrame>(spectra["header"]);
    const Rcpp::NumericVector hdMSLevels = specHeader["msLevel"];
    const Rcpp::NumericVector hdSeqNums = specHeader["seqNum"];
    const Rcpp::List specList = spectra["spectra"];
    const double minMZ = mzRange[0], maxMZ = mzRange[1];
    
    int cols = 0;
    std::vector<int> specInds;
    for (int i=0; i<hdMSLevels.size(); ++i)
    {
        if (hdMSLevels[i] != 1)
            continue;
        
        specInds.push_back(hdSeqNums[i]);
        
        int c = 0;
        const Rcpp::NumericMatrix spec = Rcpp::as<Rcpp::NumericMatrix>(specList[i]);
        for (int j=0; j<spec.nrow(); ++j)
        {
            if (numberGTE(spec(j, 0), minMZ, 1E-8) && numberLTE(spec(j, 0), maxMZ, 1E-8) &&
                spec(j, 1) != 0)
                ++c;
        }
        cols = std::max(cols, c);
    }
    
    Rcpp::NumericMatrix mzM(specInds.size(), cols), intM(specInds.size(), cols);
    Rcpp::Rcout << "creating " << specInds.size() << "x" << cols << " matrix\n";
    for (int i=0; i<specInds.size(); ++i)
    {
        Rcpp::NumericMatrix spec = Rcpp::as<Rcpp::NumericMatrix>(specList[specInds[i] - 1]);
        Rcpp::NumericVector mzs = spec(Rcpp::_, 0), ints = spec(Rcpp::_, 1);
        
        int ind = 0;
        for (int j=0; j<spec.nrow(); ++j)
        {
            if (numberGTE(mzs[j], minMZ, 1E-8) && numberLTE(mzs[j], maxMZ, 1E-8) &&
                ints[j] != 0)
            {
                mzM(i, ind) = mzs[j];
                intM(i, ind) = ints[j];
                ++ind;
            }
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("mzM") = mzM, Rcpp::Named("intM") = intM);
}

BinnedSpectrum doBinSpectra(const Spectrum &specLeft, Spectrum specRight,
                            const std::string &shift, double precDiff, double mzWindow)
{
    // assumptions: specs are ordered on mz
    
    if (shift == "precursor")
    {
        for (double &m : specRight.mzs)
            m -= precDiff;
        // NOTE: negative m/z values will be skipped below
    }
    else if (shift == "both") // UNDONE: other name for "both"?
    {
        // first bin as normal (recursive call)
        const BinnedSpectrum binNone = doBinSpectra(specLeft, specRight, "none", precDiff, mzWindow);
        Spectrum specLeftUn, specRightUn;
        BinnedSpectrum binOverlap;
        for (size_t i=0; i<binNone.mzs.size(); ++i)
        {
            const double m = binNone.mzs[i];
            if (binNone.intsLeft[i] == 0)
            {
                specRightUn.IDs.push_back(binNone.IDsRight[i]);
                specRightUn.mzs.push_back(m);
                specRightUn.intensities.push_back(binNone.intsRight[i]);
            }
            else if (binNone.intsRight[i] == 0)
            {
                specLeftUn.IDs.push_back(binNone.IDsLeft[i]);
                specLeftUn.mzs.push_back(m);
                specLeftUn.intensities.push_back(binNone.intsLeft[i]);
            }
            else
            {
                binOverlap.mzs.push_back(m);
                binOverlap.intsLeft.push_back(binNone.intsLeft[i]);
                binOverlap.intsRight.push_back(binNone.intsRight[i]);
                binOverlap.IDsLeft.push_back(binNone.IDsLeft[i]);
                binOverlap.IDsRight.push_back(binNone.IDsRight[i]);
            }
        }
        
        // bin missing with shift
        BinnedSpectrum binShift = doBinSpectra(specLeftUn, specRightUn, "precursor", precDiff, mzWindow);
        
        // merge both: add missing from binNone
        binShift.mzs.insert(binShift.mzs.end(), binOverlap.mzs.begin(), binOverlap.mzs.end());
        binShift.intsLeft.insert(binShift.intsLeft.end(), binOverlap.intsLeft.begin(), binOverlap.intsLeft.end());
        binShift.intsRight.insert(binShift.intsRight.end(), binOverlap.intsRight.begin(), binOverlap.intsRight.end());
        binShift.IDsLeft.insert(binShift.IDsLeft.end(), binOverlap.IDsLeft.begin(), binOverlap.IDsLeft.end());
        binShift.IDsRight.insert(binShift.IDsRight.end(), binOverlap.IDsRight.begin(), binOverlap.IDsRight.end());
        
        // UNDONE: sort?
        
        return binShift;
    }
    
    BinnedSpectrum ret;
    std::vector<size_t> usedRightInds;
    size_t lastRightInd = 0;

    for (size_t i=0; i<specLeft.mzs.size(); ++i)
    {
        const double leftMZ = specLeft.mzs[i];
        double rightMZ = 0, rightInt;
        int rightID;
        bool foundRight = false;
        
        while (lastRightInd < specRight.mzs.size())
        {
            const double rmz = specRight.mzs[lastRightInd], rmzmin = rmz - mzWindow, rmzmax = rmz + mzWindow;
            
            if (rmz <= 0) // may be (below) zero due to precursor shift
            {
                ++lastRightInd;
                continue;
            }
            
            if (leftMZ < rmzmin)
                break; // surpassed range for left
            
            if (leftMZ >= rmzmin && leftMZ <= rmzmax)
            {
                // overlap
                rightMZ = rmz; rightInt = specRight.intensities[lastRightInd]; rightID = specRight.IDs[lastRightInd];
                foundRight = true;
                usedRightInds.push_back(lastRightInd);
            }
            
            ++lastRightInd;
            
            if (foundRight)
                break; // done or surpassed range
        }
        
        if (foundRight)
        {
            ret.mzs.push_back((leftMZ + rightMZ) / 2.0);
            ret.intsRight.push_back(rightInt);
            ret.IDsRight.push_back(rightID);
        }
        else
        {
            ret.mzs.push_back(leftMZ);
            ret.intsRight.push_back(0);
            ret.IDsRight.push_back(0);
        }
        ret.intsLeft.push_back(specLeft.intensities[i]);
        ret.IDsLeft.push_back(specLeft.IDs[i]);
    }
    
    // add missing from right
    for (size_t j=0; j<specRight.mzs.size(); ++j)
    {
        if (std::find(usedRightInds.begin(), usedRightInds.end(), j) == usedRightInds.end() &&
            specRight.mzs[j] > 0)
        {
            ret.mzs.push_back(specRight.mzs[j]);
            ret.intsLeft.push_back(0);
            ret.IDsLeft.push_back(0);
            ret.intsRight.push_back(specRight.intensities[j]);
            ret.IDsRight.push_back(specRight.IDs[j]);
        }
    }
    
    // UNDONE: sort?
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::DataFrame binSpectra(Rcpp::DataFrame sp1, Rcpp::DataFrame sp2, Rcpp::CharacterVector shift,
                           Rcpp::NumericVector precDiff, Rcpp::NumericVector mzWindow)
{
    Spectrum specLeft{ sp1["ID"], sp1["mz"], sp1["intensity"] };
    Spectrum specRight{ sp2["ID"], sp2["mz"], sp2["intensity"] };
    
    normalizeNums(specLeft.intensities); normalizeNums(specRight.intensities);

    BinnedSpectrum binnedSpec = doBinSpectra(specLeft, specRight, Rcpp::as<std::string>(shift),
                                             Rcpp::as<double>(precDiff), Rcpp::as<double>(mzWindow));
    return Rcpp::DataFrame::create(Rcpp::Named("mz") = binnedSpec.mzs,
                                   Rcpp::Named("intensity_1") = binnedSpec.intsLeft,
                                   Rcpp::Named("intensity_2") = binnedSpec.intsRight,
                                   Rcpp::Named("ID_1") = binnedSpec.IDsLeft,
                                   Rcpp::Named("ID_2") = binnedSpec.IDsRight);
}

double doCalcSpecSimilarity(const BinnedSpectrum &binnedSpec, const std::string &method,
                            const std::string &shift, double precDiff,
                            double mzWeight, double intWeight, double mzWindow)
{
    // UNDONE: pearsons/spearman? needs sorting?
    if (method == "cosine")
    {
        std::vector<double> u, v;
        for (size_t i=0; i<binnedSpec.mzs.size(); ++i)
        {
            const double m = std::pow(binnedSpec.mzs[i], mzWeight);
            u.push_back(m * std::pow(binnedSpec.intsLeft[i], intWeight));
            v.push_back(m * std::pow(binnedSpec.intsRight[i], intWeight));
        }
        
        const double dp = std::inner_product(u.begin(), u.end(), v.begin(), 0.0);
        // sqrt(sum(u^2)) * sqrt(sum(v^2)))
        double divu = 0.0, divv = 0.0;
        for (size_t i=0; i<u.size(); ++i)
        {
            divu += (u[i] * u[i]);
            divv += (v[i] * v[i]);
        }
        
        if (divu == 0.0 || divv == 0.0) // lack of any overlap
            return 0.0;
        
        const double div = std::sqrt(divu) * std::sqrt(divv);
        
        return dp / div;
    }
    else if (method == "jaccard")
    {
        // binnedPL[intensity_1 != 0 & intensity_2 != 0, .N] / nrow(binnedPL)
        int both = 0;
        for (size_t i=0; i<binnedSpec.mzs.size(); ++i)
        {
            if (binnedSpec.intsLeft[i] != 0 && binnedSpec.intsRight[i] != 0)
                ++both;
        }
        return static_cast<double>(both) / static_cast<double>(binnedSpec.mzs.size());
    }
    
    return NA_REAL; // shouldn't be here
}

double doCalcSpecSimilarity(Spectrum sp1, Spectrum sp2, const std::string &method,
                            const std::string &shift, double precDiff,
                            double mzWeight, double intWeight, double mzWindow)
{
    normalizeNums(sp1.intensities); normalizeNums(sp2.intensities);
    const BinnedSpectrum binnedSpec = doBinSpectra(sp1, sp2, shift, precDiff, mzWindow);
    return doCalcSpecSimilarity(binnedSpec, method, shift, precDiff, mzWeight, intWeight, mzWindow);
}

// [[Rcpp::export]]
Rcpp::NumericVector calcSpecSimilarity(Rcpp::DataFrame sp1, Rcpp::DataFrame sp2, Rcpp::CharacterVector method,
                                       Rcpp::CharacterVector shift, Rcpp::NumericVector precDiff,
                                       Rcpp::NumericVector mzWeight, Rcpp::NumericVector intWeight, Rcpp::NumericVector mzWindow)
{
    const Spectrum specLeft{ sp1["ID"], sp1["mz"], sp1["intensity"] };
    const Spectrum specRight{ sp2["ID"], sp2["mz"], sp2["intensity"] };
    return Rcpp::NumericVector::create(doCalcSpecSimilarity(specLeft, specRight, Rcpp::as<std::string>(method),
                                                            Rcpp::as<std::string>(shift), Rcpp::as<double>(precDiff),
                                                            Rcpp::as<double>(mzWeight), Rcpp::as<double>(intWeight),
                                                            Rcpp::as<double>(mzWindow)));
}

// [[Rcpp::export]]
std::vector<double> calcAnnSims(Rcpp::DataFrame spectrum, Rcpp::List annotatedInds, const std::string &method,
                                double mzWeight, double intWeight, double mzWindow)
{
    // Create a dummy bin spectrum: the reference spectrum is compared to itself, minus the unannotated peaks. The right
    // spectrum is the one with just the annotated peaks. The intensities of left/right are the same, but will only be
    // set below if the peak was annotated.
    BinnedSpectrum binsp;
    binsp.mzs = Rcpp::as<std::vector<double>>(spectrum["mz"]);
    binsp.intsLeft = Rcpp::as<std::vector<double>>(spectrum["intensity"]);
    binsp.IDsLeft = binsp.IDsRight = Rcpp::as<std::vector<int>>(spectrum["ID"]);
    
    std::vector<double> ret(annotatedInds.size());
    for (int i=0; i<annotatedInds.size(); ++i)
    {
        const Rcpp::IntegerVector ai = Rcpp::as<Rcpp::IntegerVector>(annotatedInds[i]);
        const std::set<int> annInds(ai.begin(), ai.end());

        binsp.intsRight.clear();        
        for (size_t i=0; i<binsp.IDsLeft.size(); ++i)
        {
            if (annInds.find(binsp.IDsLeft[i]) != annInds.end())
                binsp.intsRight.push_back(binsp.intsLeft[i]);
            else
                binsp.intsRight.push_back(0.0);
        }
        
        ret[i] = doCalcSpecSimilarity(binsp, method, "none", 0.0, mzWeight, intWeight, mzWindow);
    }
    
    return ret;
}
