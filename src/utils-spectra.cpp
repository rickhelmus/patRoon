#include <string>
#include <vector>

#include <Rcpp.h>

#include "utils.h"

// [[Rcpp::export]]
Rcpp::NumericVector loadEICIntensities(Rcpp::List spectra, Rcpp::DataFrame featList, Rcpp::NumericVector rtWindow)
{
    const Rcpp::DataFrame specHeader = spectra["header"];
    const Rcpp::NumericVector hdRetTimes = specHeader["retentionTime"];
    const Rcpp::NumericVector hdMSLevels = specHeader["msLevel"];
    const Rcpp::NumericVector hdSeqNums = specHeader["seqNum"];
    
    const Rcpp::List specList = spectra["spectra"];
    
    const Rcpp::NumericVector rets = featList["ret"];
    const Rcpp::NumericVector mzmins = featList["mzmin"];
    const Rcpp::NumericVector mzmaxs = featList["mzmax"];
    
    const double rtWin = Rcpp::as<double>(rtWindow);
    const int specCount = hdRetTimes.size();
    const int ftCount = rets.length();
    
    Rcpp::NumericVector intensities(ftCount);
    for (int fti=0; fti<ftCount; ++fti)
    {
        const double featRt = rets[fti];
        const double rtMin = featRt - rtWin;
        const double rtMax = featRt + rtWin;
        const double mzMin = mzmins[fti], mzMax = mzmaxs[fti];
        
        //std::vector<double> EICInts;
        double closestRTDiff = -1, closestInt = 0;
        for (int spi=0; spi<specCount; ++spi)
        {
            const double specRt = hdRetTimes[spi];
            
            if (!numberWithin(specRt, rtMin, rtMax, 1E-4) || hdMSLevels[spi] != 1)
                continue;
            
            const double rtdiff = fabs(featRt - specRt);
            if (closestRTDiff == -1)
                closestRTDiff = rtdiff;
            else if (closestRTDiff < rtdiff)
                continue;
            
            const Rcpp::DataFrame peaklist = specList[hdSeqNums[spi] - 1]; // -1: R 1 based index
            const Rcpp::NumericVector peakMZs = peaklist["mz"], peakInts = peaklist["intensity"];
            const int pCount = peakMZs.length();
            
            double totInt = 0;
            for (int pi=0; pi<pCount; ++pi)
            {
                const double mz = peakMZs[pi];
                
                if (numberWithin(mz, mzMin, mzMax, 1E-8))
                    totInt += peakInts[pi];
            }
            
            closestRTDiff = rtdiff;
            closestInt = totInt;
        }
        
        intensities[fti] = closestInt;
    }
    
    return intensities;
}
