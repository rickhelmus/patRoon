#include <string>
#include <vector>

#include <Rcpp.h>

#include "utils.h"

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
Rcpp::List loadEICs(Rcpp::List spectra, Rcpp::List rtRanges, Rcpp::List mzRanges)
{
    const int EICCount = rtRanges.length();
    Rcpp::List ret(EICCount);
    for (int eici=0; eici<EICCount; ++eici)
    {
        const Rcpp::NumericVector rtr = rtRanges[eici], mzr = mzRanges[eici];
        const double mzMin = mzr[0], mzMax = mzr[1];
        
        std::vector<double> EICTimes;
        std::vector<double> EICIntensities;
        specApply(spectra, rtr[0], rtr[1], mzMin, mzMax,
                  [&](double specRt, const Rcpp::NumericVector &peakMZs, const Rcpp::NumericVector &peakInts)
                  {
                      EICTimes.push_back(specRt);
                      EICIntensities.push_back(getTotMZIntFromSpec(peakMZs, peakInts, mzMin, mzMax));
                  });

        // compress data by removing any zero intensity datapoints that are inbetween two other zero intensitiy points.
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
