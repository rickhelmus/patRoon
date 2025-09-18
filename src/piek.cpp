#include <Rcpp.h>

#include "spectrum-raw.h"
#include "msdata.h"

namespace {

/*Supporting Information*/
/*Open source feature detection for non-target LC-MS analytics*/
/*Christian Dietrich, Arne Wick, Thomas A. Ternes*/

/*licensed under the GNU GPL version 3 or later*/

/*This algorithm is designed to process centroided high-resolution MS data provided as mzML or mzXML*/

/* The function was modified by Rick Helmus:
 * Small formatting (eg indenting)
 * Return struct instead of Rcpp objects
 * Convert vector args to const refs
 * Removed unnecessary RcppArmadillo dependencies
 * Converted RT min/max args from scans to times and disable check if its zero
 * Disabled noise removal for reported intensities/areas
 * Converted peakwidth_min and max from ints to doubles
 * Use SpectrumRawTypes for intensity and time vectors
 * Few small fixes, marked in the code below
 * 
 * UNDONE: convert int args to doubles?
 */ 

struct PeakPickingResults
{
    std::vector<double> times, timesLeft, timesRight, XICIntensities, noiseDeviations, areas, FWHMsLeft, FWHMsRight,
                        baselines;
    std::vector<int> scans, scansLeft, scansRight, noiseScans;
    PeakPickingResults(void) = default;
    PeakPickingResults(size_t s) : times(s), timesLeft(s), timesRight(s), XICIntensities(s), noiseDeviations(s),
                                   areas(s), FWHMsLeft(s), FWHMsRight(s), baselines(s), scans(s), scansLeft(s),
                                   scansRight(s), noiseScans(s) { }
};

PeakPickingResults peakPicking_cpp(const std::vector<SpectrumRawTypes::Intensity> &intensity, const std::vector<SpectrumRawTypes::Time> &scantime,
                                   double min_intensity, int sn, double peakwidth_min, double peakwidth_max, double rt_min,
                                   double rt_max, int maxPeaksPerSignal)
{
    
    std::vector<double> derivative((intensity.size()-1));
    std::vector<int> maxima(intensity.size(),0);
    int anzahlmaxima=0;
    int j = 0;
    double min_intensity2 = 0.1;
    
    /*calculate the lower 10-quantile of all XIC intensities and intermediately set the minimum intensity to this value (unless you did not put it even lower) */
    /*makes sure to not overlook small peaks later on and erroneously consider them as noise*/  
    std::vector<double> intensities(intensity.size(),0);
    intensities.assign(intensity.begin(),intensity.end());
    std::sort(intensities.begin(),intensities.end());
    min_intensity2 = intensities[intensities.size()*0.9];
    if ((min_intensity2 == 0) || (min_intensity2 > min_intensity)) {
        min_intensity2 = min_intensity;
    }
    
    /*calculate the first derivative of the XIC intensities and check for local maxima in each scan*/  
    derivative[0] = intensity[1]-intensity[0];
    for(int i = 1; i < (intensity.size()-1); ++i) {
        derivative[i] = intensity[i+1]-intensity[i];
        if ((derivative[(i-1)] > 0) && (derivative[i] <= 0) /*this is basically checking 1.derivative = 0 and 2.derivative < 0 in one step*/
    && (intensity[i] >= min_intensity2)) {
            maxima[anzahlmaxima] = i;
            anzahlmaxima++;
        }
    }
    
    maxima.resize(anzahlmaxima);
    std::vector<int> left_end(anzahlmaxima,0);	 
    std::vector<int> right_end(anzahlmaxima,0);
    std::vector<double> noiselevel(anzahlmaxima,0);
    std::vector<int> amountofpeaks(anzahlmaxima,0);
    std::vector<int> noisescans(anzahlmaxima,0);
    double noiseintensity=0;
    int noisecounter=0;
    
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        /* find start of peak (left_end) */
        j = maxima[i]-1;
        while ((derivative[j] > 0) && (j > 0)) {
            j--;
        }
        left_end[i] = j+1; 
        /* find end of peak (right_end) */
        j = maxima[i];
        while ((derivative[j] < 0) && (j < intensity.size())) {
            j++;
        }
        right_end[i] = j; 
        
        /* 1st noiselevel calculation based on mean intensity around the peak */
        noiseintensity = 0;
        noisecounter = 0;  
        
        noisescans[i] = (right_end[i]-left_end[i])*10;	    	
        j = left_end[i];
        while ((j > (left_end[i]-noisescans[i])) && (j > 0)) {
            j--;
            noisecounter++;
            noiseintensity += intensity[j];
        }
        
        
        noiselevel[i] = noiseintensity/noisecounter;
        
        amountofpeaks[i] = 1;
        /* check if two peaks belong together */
        if ((i > 0) && (maxima[i] > 0)) {
            if ((right_end[i-1] >= left_end[i]) && (maxima[i-1] > 0)) {
                if ((intensity[left_end[i]]-noiselevel[i-1]) > (std::min(intensity[maxima[i]],intensity[maxima[i-1]])-noiselevel[i-1])/2) {
                    if (std::min(intensity[maxima[i-1]],intensity[maxima[i]])-noiselevel[i-1] > (std::max(intensity[maxima[i-1]],intensity[maxima[i]])-noiselevel[i-1])/2) {
                        amountofpeaks[i] = amountofpeaks[i]+amountofpeaks[i-1];
                    } else {
                        if (amountofpeaks[i-1] > maxPeaksPerSignal) {
                            noiselevel[i-1] = intensity[maxima[i-1]];
                            left_end[i-1] = left_end[i];
                            maxima[i-1] = maxima[i]; // ADDED: FIX: don't change maxima of this peak in the if line below as it is considered noise  
                        }
                        if ((intensity[maxima[i-1]] > intensity[maxima[i]]) && (anzahlmaxima > i)) { // CHANGED: fixed single '&' usage
                            if (std::min(intensity[maxima[i]],intensity[maxima[i+1]])-noiselevel[i] > (std::max(intensity[maxima[i]],intensity[maxima[i+1]])-noiselevel[i])/2) {
                                right_end[i] = right_end[i-1];
                                amountofpeaks[i] = amountofpeaks[i-1];
                            }
                        }
                    }
                    if (intensity[maxima[i-1]] > intensity[maxima[i]]) maxima[i] = maxima[i-1];
                    noiselevel[i] = noiselevel[i-1];
                    left_end[i] = left_end[i-1];
                    noisescans[i] = (right_end[i]-left_end[i])*10;
                    noisescans[i-1] = 0;
                    maxima[i-1] = 0;
                    left_end[i-1] = 0;
                    right_end[i-1] = 0;
                    noiselevel[i-1] = 0;
                    amountofpeaks[i-1] = 0;
                } else { 
                    if (amountofpeaks[i-1] < maxPeaksPerSignal) noiselevel[i] = noiselevel[i-1];
                }
            }
        }
        
    }
    
    
    /* delete maxima that have been set to 0, those below intensity threshold and too narrow ones*/
    j = 0;
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        if ((maxima[i] > 0) && (intensity[maxima[i]] > min_intensity) && (scantime[right_end[i]]-scantime[left_end[i]] > peakwidth_min)) {
            maxima[j] = maxima[i];
            left_end[j] = left_end[i];
            right_end[j] = right_end[i];
            noiselevel[j] = noiselevel[i];
            noisescans[j] = noisescans[i];
            amountofpeaks[j] = amountofpeaks[i];
            j++;
        }
    }
    anzahlmaxima = j;
    maxima.resize(anzahlmaxima); 
    
    
    
    int peaksectionstartid = 0;
    int amountofpeaks_sum = 0;
    
    /* check if there are peak clusters */
    for(int i = 1; i < (anzahlmaxima); ++i) { 
        amountofpeaks_sum += amountofpeaks[i-1];
        if ((right_end[i-1] < left_end[i]) || (std::min(intensity[maxima[i]],intensity[maxima[i-1]])-noiselevel[i-1] < (std::max(intensity[maxima[i-1]],intensity[maxima[i]])-noiselevel[i-1])/2) || ((maxima[i]-maxima[i-1]) > noisescans[i]) || (i==anzahlmaxima-1)) {
            if ((amountofpeaks_sum > maxPeaksPerSignal) && (peaksectionstartid != i)) {
                for(int n = peaksectionstartid; n < i; n++) {
                    maxima[n] = 0;
                }
            }
            peaksectionstartid = i;
            amountofpeaks_sum = 0;	
        }
        
    }
    
    /* delete maxima that have been set to 0 */
    j = 0;
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        if (maxima[i] > 0)  {
            maxima[j] = maxima[i];
            left_end[j] = left_end[i];
            right_end[j] = right_end[i];
            noiselevel[j] = noiselevel[i];
            noisescans[j] = noisescans[i];
            amountofpeaks[j] = amountofpeaks[i];
            j++;
        }
    }
    anzahlmaxima = j;
    maxima.resize(anzahlmaxima); 
    
    
    /* accurate noise calculation */
    std::vector<double> noisedeviation(anzahlmaxima,0);
    std::vector<double> noiseRegion(intensity.size(),0); 
    int otherMaximum;
    
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        noiseintensity = 0;
        noisecounter = 0;  
        std::fill(noiseRegion.begin(),noiseRegion.end(),0);
        
        /* check if there are peaks in front of this one, otherwise add the respective scan to the noiseRegion*/
        j = left_end[i];
        otherMaximum = i-1;
        while ((j > (left_end[i]-noisescans[i])) && (j > 0)) {
            if ((otherMaximum >= 0) && (right_end[otherMaximum] >= j)) {
                /* if there is another peak, jump to the beginning of that peak*/
                j = left_end[otherMaximum];
                otherMaximum -= 1;
            } else {
                j--;
                noiseRegion[noisecounter] = intensity[j];
                noisecounter++;
                noiseintensity += intensity[j];
            }
        }
        
        /* check if there are peaks after this one, otherwise add the respective scan to the noiseRegion*/
        j = right_end[i];
        otherMaximum = i+1;
        while ((j < (right_end[i]+noisescans[i])) && (j < intensity.size()-1)) {
            if ((otherMaximum < anzahlmaxima) && (left_end[otherMaximum] <= j)) {
                /* if there is another peak, jump to the beginning of that peak*/
                j = right_end[otherMaximum];
                otherMaximum += 1;
            } else {
                j++;
                noiseRegion[noisecounter] = intensity[j];
                noisecounter++;
                noiseintensity += intensity[j];
            }
        }
        std::sort(noiseRegion.begin(),noiseRegion.begin()+noisecounter);
        noisedeviation[i] = noiseRegion[noisecounter*0.9]-noiseRegion[noisecounter*0.1];
        noiselevel[i] = noiseintensity/noisecounter;
    } 
    
    /* delete maxima that are outside predefined RT range and those below S/N treshold */
    j = 0;
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        if ((scantime[maxima[i]] >= rt_min) && (rt_max == 0.0 || scantime[maxima[i]] <= rt_max) && ((intensity[maxima[i]]-noiselevel[i])*2 >= noisedeviation[i]*sn))  {	
            maxima[j] = maxima[i];
            left_end[j] = left_end[i];
            right_end[j] = right_end[i];
            noiselevel[j] = noiselevel[i];
            noisescans[j] = noisescans[i];
            noisedeviation[j] = noisedeviation[i];
            amountofpeaks[j] = amountofpeaks[i];
            j++;
        }
    }
    anzahlmaxima = j;
    maxima.resize(anzahlmaxima); 
    
    
    /*calculate FWHM*/
    std::vector<double> FWHM_left(anzahlmaxima,0);
    std::vector<double> FWHM_right(anzahlmaxima,0);
    double slope = 0;
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        j = left_end[i];
        while (intensity[j]-noiselevel[i] < (intensity[maxima[i]]-noiselevel[i])/2) {
            j++;
        }
        slope = (intensity[j]-intensity[j-1])/(scantime[j]-scantime[j-1]);
        if (slope == 0) {
            FWHM_left[i] = scantime[j];
        } else {
            FWHM_left[i] = scantime[j-1]+(intensity[maxima[i]]/2-intensity[j-1]+noiselevel[i])/slope;
        }
        
        
        j = right_end[i];
        while (intensity[j]-noiselevel[i] < (intensity[maxima[i]]-noiselevel[i])/2) {
            j--;
        }
        slope = (intensity[j+1]-intensity[j])/(scantime[j+1]-scantime[j]);
        if (slope == 0) {
            FWHM_right[i] = scantime[j];
        } else {
            FWHM_right[i] = scantime[j]+(intensity[maxima[i]]/2-intensity[j]+noiselevel[i])/slope;
        }
    }
    
    
    
    /* delete too broad peaks (2 x FWHM > peakwidth_max) and calculate area*/
    j = 0;
    std::vector<double> area(anzahlmaxima,0);
    
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        if ((FWHM_right[i]-FWHM_left[i])*2 <= peakwidth_max)  {
            maxima[j] = maxima[i];
            left_end[j] = left_end[i];
            right_end[j] = right_end[i];
            noiselevel[j] = noiselevel[i];
            noisedeviation[j] = noisedeviation[i];
            FWHM_left[j]= FWHM_left[i];
            FWHM_right[j] = FWHM_right[i];
            amountofpeaks[j] = amountofpeaks[i];
            for (int ii = left_end[i]; ii < right_end[i]; ++ii) {
                // CHANGED: don't subtract noise
                //if ((intensity[ii] >= noiselevel[i]) && (intensity[ii+1] >= noiselevel[i])) area[j] += ((intensity[ii]-noiselevel[i])+(intensity[ii+1]-noiselevel[i]))*(scantime[ii+1]-scantime[ii])/2;
                area[j] += (intensity[ii]+intensity[ii+1])*(scantime[ii+1]-scantime[ii])/2;
            }
            j++;
        }
    }
    anzahlmaxima = j;
    maxima.resize(anzahlmaxima); 
    
    PeakPickingResults ret(anzahlmaxima);
    for(int i = 0; i < (anzahlmaxima); ++i)
    {
        ret.times[i] = scantime[maxima[i]];
        ret.XICIntensities[i] = intensity[maxima[i]]; // -noiselevel[i]; --> CHANGED: report raw intensities
        ret.scans[i] = maxima[i]+1;
        ret.timesLeft[i] = scantime[left_end[i]];
        ret.timesRight[i] = scantime[right_end[i]];
        ret.scansLeft[i] = left_end[i]+1;
        ret.scansRight[i] = right_end[i]+1;
        ret.noiseScans[i] = noisescans[i];
        ret.noiseDeviations[i] = noisedeviation[i];
        ret.areas[i] = area[i];
        ret.FWHMsLeft[i] = FWHM_left[i];
        ret.FWHMsRight[i] = FWHM_right[i];
        ret.baselines[i] = noiselevel[i];
    }
    
    return ret;
    
#if 0
    NumericMatrix ergebnis(anzahlmaxima,17);
    
    for(int i = 0; i < (anzahlmaxima); ++i) { 
        ergebnis(i,0) = 0; /*will be filled in R*/
    ergebnis(i,1) = scantime[maxima[i]];
    ergebnis(i,2) = 0; /*will be filled in R*/
    ergebnis(i,3) = intensity[maxima[i]]-noiselevel[i];
    ergebnis(i,4) = maxima[i]+1;
    ergebnis(i,5) = scantime[left_end[i]];
    ergebnis(i,6) = scantime[right_end[i]];
    ergebnis(i,7) = left_end[i]+1;
    ergebnis(i,8) = right_end[i]+1;
    ergebnis(i,9) = noisescans[i];
    ergebnis(i,10) = noisedeviation[i];
    ergebnis(i,11) = area[i];
    ergebnis(i,12) = FWHM_left[i];
    ergebnis(i,13) = FWHM_right[i];
    ergebnis(i,14) = noiselevel[i];
    ergebnis(i,15) = 0; /*will be filled in R*/
    ergebnis(i,16) = 0; /*will be filled in R*/
    
    c("mz","RT","Intensity","XIC_Intensity","Scan","LeftendRT","RightendRT","Leftendscan","Rightendscan","NoiseScans", "NoiseDeviation","Area","FWHM_left","FWHM_right","Baseline","XICStartMass","MS2scan")
    }
    
    
    return(ergebnis);
#endif
}

}

// [[Rcpp::export]]
Rcpp::List doFindPeaksPiek(Rcpp::List EICs, bool fillEICs, double minIntensity, int SN, double peakWidthMin,
                           double peakWidthMax, double RTMin, double RTMax, int maxPeaksPerSignal,
                           bool verbose = true)
{
    // UNDONE: set default args
    
    const size_t entries = EICs.size();
    
    if (entries == 0)
        return Rcpp::List();

    const std::vector<SpectrumRawTypes::Time> fullEICTimes = (fillEICs) ? EICs.attr("allXValues") : std::vector<SpectrumRawTypes::Time>();
    std::vector<std::vector<SpectrumRawTypes::Time>> allTimes;
    std::vector<std::vector<SpectrumRawTypes::Intensity>> allIntensities;
    for (size_t i=0; i<entries; ++i)
    {
        Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(EICs[i]);
        allTimes.emplace_back(Rcpp::as<std::vector<SpectrumRawTypes::Time>>(df["time"]));
        allIntensities.emplace_back(Rcpp::as<std::vector<SpectrumRawTypes::Intensity>>(df["intensity"]));
    }
    
    std::vector<PeakPickingResults> ppresults(entries);
    
    if (verbose) // UNDONE: progress bar?
        Rcpp::Rcout << "Finding peaks for " << entries << " traces...";
    
    if (fillEICs)
    {
        #pragma omp parallel for
        for (size_t i=0; i<entries; ++i)
        {
            const auto intens = fillEIXIntensities(fullEICTimes, allTimes[i], allIntensities[i]);
            ppresults[i] = peakPicking_cpp(intens, fullEICTimes, minIntensity, SN, peakWidthMin, peakWidthMax,
                                           RTMin, RTMax, maxPeaksPerSignal);
        }
    }
    else
    {
        #pragma omp parallel for
        for (size_t i=0; i<entries; ++i)
        {
            ppresults[i] = peakPicking_cpp(allIntensities[i], allTimes[i], minIntensity, SN, peakWidthMin, peakWidthMax,
                                           RTMin, RTMax, maxPeaksPerSignal);
        }
    }

    if (verbose)
        Rcpp::Rcout << " Done!\n";
    
    Rcpp::List ret(entries);
    for (size_t i=0; i<entries; ++i)
    {
        ret[i] = Rcpp::List::create(Rcpp::Named("ret") = ppresults[i].times,
                                    Rcpp::Named("retmin") = ppresults[i].timesLeft,
                                    Rcpp::Named("retmax") = ppresults[i].timesRight,
                                    Rcpp::Named("area") = ppresults[i].areas,
                                    Rcpp::Named("intensity") = ppresults[i].XICIntensities,
                                    Rcpp::Named("noiseDeviation") = ppresults[i].noiseDeviations,
                                    Rcpp::Named("FWHMLeft") = ppresults[i].FWHMsLeft,
                                    Rcpp::Named("FWHMRight") = ppresults[i].FWHMsRight,
                                    Rcpp::Named("baseline") = ppresults[i].baselines,
                                    Rcpp::Named("scan") = ppresults[i].scans,
                                    Rcpp::Named("scanLeft") = ppresults[i].scansLeft,
                                    Rcpp::Named("scanRight") = ppresults[i].scansRight,
                                    Rcpp::Named("noiseScans") = ppresults[i].noiseScans);
    }
    
    return ret;
}


// utils for findFeaturesPiek()

// [[Rcpp::export]]
Rcpp::LogicalVector findFeatTableDups(const Rcpp::NumericVector &rts, const Rcpp::NumericVector &mzs,
                                      const Rcpp::NumericVector &mobs, const Rcpp::NumericVector &ints,
                                      double tolRT, double tolMZ, double tolMob, const Rcpp::LogicalVector &mzCentered,
                                      const Rcpp::LogicalVector &mobCentered, double prefIntensityRatio)
{
    // NOTE: rts, mobs and ints are optional
    
    Rcpp::LogicalVector ret(mzs.size(), false);
    
    if (mzs.size() < 2)
        return ret; // no duplicates possible
    
    const auto isSame = [](double a, double b, double tol)
    {
        // handling of NAs
        // both NA: ignore dimension in comparison
        // one NA: don't compare, skip
        // both not NA: compare values
        
        if (Rcpp::NumericVector::is_na(a) && Rcpp::NumericVector::is_na(b))
            return true;
        if (Rcpp::NumericVector::is_na(a) || Rcpp::NumericVector::is_na(b))
            return false;
        return numberLTE(std::abs(a - b), tol);
    };

    const bool hasMobs = mobs.size() > 0;
    const auto intSortedInds = getSortedInds(ints, [&ints](size_t i, size_t j) { return ints[i] > ints[j]; });
    
    for (size_t i=0; i<intSortedInds.size()-1; ++i)
    {
        const auto indi = intSortedInds[i];
        const bool mzCenti = mzCentered[indi];
        const bool mobCenti = (hasMobs) ? mobCentered[indi] : true;
        const double inteni = ints[indi];
        bool needCentered = !mzCenti || !mobCenti;
        size_t partialCentInd = (mzCenti || mobCenti) ? indi : std::numeric_limits<size_t>::max();
        for (size_t j=i+1; j<intSortedInds.size(); ++j)
        {
            const auto indj = intSortedInds[j];
            
            if (ret[indj]) // already marked as duplicate
                continue;
            
            if (!isSame(rts[indi], rts[indj], tolRT))
                continue;
            if (!isSame(mzs[indi], mzs[indj], tolMZ))
                continue;
            if (hasMobs && !isSame(mobs[indi], mobs[indj], tolMob))
                continue;
            
            const bool mzCentj = mzCentered[indj];
            const bool mobCentj = (hasMobs) ? mobCentered[indj] : true;
            const double intenj = ints[indj];
            
            if (!needCentered || (inteni * prefIntensityRatio) > intenj)
                ret[indj] = true; // feat i is centered or considerably more intense, so just mark all others as duplicate
            else if (mzCentj && mobCentj) // i is not centered, but this is --> take this instead
            {
                if (partialCentInd != std::numeric_limits<size_t>::max())
                    ret[partialCentInd] = true; // mark previous partial centered as duplicate
                ret[indi] = true;
                needCentered = false;
            }
            else if (partialCentInd == std::numeric_limits<size_t>::max() && (mzCentj || mobCentj))
            {
                partialCentInd = indj;
                ret[indi] = true; // prefer the partially centered one
            }
            else
                ret[indj] = true; // we already have a more intense partially intense centered, so mark this as duplicate
        }
    }
    
    return ret;
}

// [[Rcpp::export]]
Rcpp::LogicalVector filterEICBins(const Rcpp::NumericVector &binStartMZs, const double mzBinWidth,
                                  const Rcpp::NumericVector &binStartMobs, const double mobBinWidth,
                                  const Rcpp::NumericVector &checkStartMZs, const Rcpp::NumericVector &checkEndMZs,
                                  const Rcpp::NumericVector &checkStartMobs, const Rcpp::NumericVector &checkEndMobs)
    
{
    // NOTE: bins should be sorted by start MZs
    // NOTE: if bins contain mobility information, the mzs are repeated and mobilities are sorted
    
    Rcpp::LogicalVector keep(binStartMZs.size(), false);
    if (binStartMZs.size() == 0 || checkStartMZs.size() == 0)
        return keep;
    
    const bool hasMobs = binStartMobs.size() > 0 && checkStartMobs.size() > 0;
    const auto sortedCheckStartMZInds = getSortedInds(checkStartMZs);
    auto prevBinIt = binStartMZs.begin();
    
    for (size_t i=0; i<sortedCheckStartMZInds.size(); ++i)
    {
        const size_t j = sortedCheckStartMZInds[i];
        auto binIt = prevBinIt;
        if (i == 0 || *binIt > checkEndMZs[j] || (*binIt + mzBinWidth) < checkStartMZs[j])
            binIt = std::lower_bound(prevBinIt, binStartMZs.end(), checkStartMZs[j] - mzBinWidth);
        
        if (binIt == binStartMZs.end())
            break; // no more bins to check
        
        if (hasMobs)
        {
            // For IMS checks: start here again for next check item as it may have different mobilities, otherwise we
            // continue after all bins that will be checked below
            prevBinIt = binIt;
        }
        
        while (binIt!=binStartMZs.end() && numberLTE(*binIt, checkEndMZs[j]))
        {
            const size_t binIndex = std::distance(binStartMZs.begin(), binIt);
            
            if (!numberLTE(checkStartMZs[j], binStartMZs[binIndex] + mzBinWidth) ||
                !numberGTE(checkEndMZs[j], binStartMZs[binIndex]))
            {
                // didn't reach proper mz bin yet
                ++binIt;
                continue;
            }
            
            if (hasMobs)
            {
                const auto binMZ = binStartMZs[binIndex];
                const auto endMZBinIt = std::upper_bound(binIt, binStartMZs.end(), binMZ);
                
                if (Rcpp::NumericVector::is_na(checkStartMobs[j]))
                {
                    // no mobility check, just keep all bins in this MZ bin
                    for (auto ind = binIndex; binIt!=endMZBinIt; ++binIt, ++ind)
                        keep[ind] = true;
                }
                else
                {
                    const auto endMZBinInd = std::distance(binStartMZs.begin(), endMZBinIt);
                    const auto mobEndIt = std::next(binStartMobs.begin(), endMZBinInd);
                    
                    auto mobIt = std::lower_bound(std::next(binStartMobs.begin(), binIndex), mobEndIt,
                                                  checkStartMobs[j] - mobBinWidth);
                    for (; mobIt!=mobEndIt && numberLTE(*mobIt, checkEndMobs[j]); ++mobIt)
                    {
                        if (!numberLTE(checkStartMobs[j], *mobIt + mobBinWidth) || !numberGTE(checkEndMobs[j], *mobIt))
                            continue;
                        const size_t mobIndex = std::distance(binStartMobs.begin(), mobIt);
                        // Rcpp::Rcout << "Found it!\n";
                        keep[mobIndex] = true;
                    }
                    binIt = endMZBinIt; // skip to end of MZ bin
                }
            }
            else
            {
                keep[binIndex] = true;
                ++binIt;
            }
        }
        
        if (!hasMobs)
            prevBinIt = binIt;
    }

    return keep;
}

// [[Rcpp::export]]
Rcpp::LogicalVector filterPiekResults(const Rcpp::NumericVector &resultRTs, const Rcpp::NumericVector &resultMZs,
                                      const Rcpp::NumericVector &resultMobs, const Rcpp::NumericVector &checkRTs,
                                      const Rcpp::NumericVector &checkStartMZs, const Rcpp::NumericVector &checkEndMZs,
                                      const Rcpp::NumericVector &checkStartMobs,
                                      const Rcpp::NumericVector &checkEndMobs, double tolRT)
{
    // NOTE: mz range may be isolation ranges with unknown centers, so these are given as start/end range
    // NOTE: RT/mob checking is only done if check vectors are given
    
    Rcpp::LogicalVector keep(resultMZs.size(), false);
    const auto sortedStartMZInds = getSortedInds(checkStartMZs);
    const bool doCheckRTs = checkRTs.size() > 0, doCheckMobs = checkStartMobs.size() > 0;
    
    double maxMZWidth = 0.0;
    for (size_t i=0; i<checkStartMZs.size(); ++i)
        maxMZWidth = std::max(maxMZWidth, checkEndMZs[i] - checkStartMZs[i]);
    
    for (int i=0; i<resultMZs.size(); ++i)
    {
        const auto mz = resultMZs[i];
        auto it = std::lower_bound(sortedStartMZInds.begin(), sortedStartMZInds.end(), mz - maxMZWidth,
                                   [&checkStartMZs](size_t ind, double m) { return checkStartMZs[ind] < m; });
        for (; it!=sortedStartMZInds.end() && numberLTE(checkStartMZs[*it], mz); ++it)
        {
            const size_t j = *it;
            if (!numberGTE(checkEndMZs[j], mz))
                continue;
            if (doCheckRTs && !numberLTE(std::abs(checkRTs[j] - resultRTs[i]), tolRT))
                continue;
            if (doCheckMobs && (!numberLTE(checkStartMobs[j], resultMobs[i]) ||
                !numberGTE(checkEndMobs[j], resultMobs[i])))
                continue;
            
            keep[i] = true;
            break; // found a match
        }
    }
    
    return keep;
}
