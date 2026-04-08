// qalgorithms_datatypes.h
#ifndef QALGORITHMS_DATATYPE_PEAK_H
#define QALGORITHMS_DATATYPE_PEAK_H

#include <vector>
#define _USE_MATH_DEFINES
#include <math.h> // INFINITY and other number macros
#include <string>
#include "qalgorithms_utils.h"

/*  This file includes the structs used for data management in qAlgorithms.
    Anything required by multiple parts of the full program should be listed here.
*/

namespace qAlgorithms
{
    // handle polarity switching
    enum Polarities
    {
        unknown_polarity,
        positive,
        negative,
        mixed,
    };

    struct CompoundFilter
    {
        // this struct is used to limit the amount of operations performed by the program to
        // a selection of masses and RTs that is relevant to the analysis at hand.
        double mz_expected[16] = {0}; // @todo reasonable amount of points?
        double mz_tolerance_ppm = 0;  // @todo is ppm the only relevant choice?

        double RT = -1;
        double RT_tol = -1; // tolerance to either side, assumes symmetrical peaks

        std::string compoundName = ""; // @todo better solution?
        std::string methodName = "";
        Polarities polarity = Polarities::unknown_polarity;
    };

    struct slice_d
    {
        double *const start;
        size_t length;
    };

    struct ProfileBlock
    {
        const float *intensity;
        const float *mz;
        size_t startPos;
        size_t length;
    };

    struct RegCoeffs
    {
        double b0 = 0, b1 = 0, b2 = 0, b3 = 0;
        size_t scale = 0, x0 = 0;
    };

    struct RegressionGauss
    {
        RegCoeffs coeffs = {0};   // regression coefficients
        Range_i regSpan = {0, 0}; // limits of the peak regression window
        int df = 0;               // degrees of freedom, interpolated data points will not be considered
        float apex_position = 0;  // position of the apex of the peak
        float mse = 0;            // mean squared error
        float area = 0;           // area of the peak (in evenly spaced x dimension, scaled later)
        float area_uncert = 0, position_uncert = 0, height_uncert = 0;
        int numCompetitors = 0; // number of points that were discarded in favour of this regression
        float jaccard = 0;
        bool isValid = false; // flag to indicate if the regression is valid
    };

    // The distinction between centroid and feature is not really sensible as a core part of the project
    struct PeakFit
    {
        RegCoeffs coeffs = {0};
        float position = 0;
        float position_uncert = 0;
        float height = 0;
        float height_uncert = 0;
        float fwhm = 0;
        float area = 0;
        float area_uncert = 0;
        float dqs = 0;
        float jaccard = 0;
    };

    struct ProfilePos // gives the range of points covered by a centroid and the access index for streamfind
    {
        unsigned int ID_spectrum = 0;
        unsigned int start = 0, end = 0; // start and end into the original MS1 spectrum
        int start_rel = 0, end_rel = 0;  // start and end in the abstract dimension where a peak is centered on 0
    };

    struct CentroidPeak
    {
        double mz = 0;
        float RT = 0;
        float height = 0;
        float area = 0;
        float width = 0;
        float heightUncertainty = 0, areaUncertainty = 0, mzUncertainty = 0;
        float DQSC = 0;
        // the binning tolerates at most three non-occurrences of a mass in order, but should not include interpolated spectra for this.
        // for conversion, number_MS1 is also the index into a vector that stores the "corrected" scan numbers after interpolation
        unsigned int number_MS1 = 0;
        // unsigned int df = 0; // degrees of freedom
        ProfilePos trace = {0};
        unsigned int numCompetitors = 0;
        unsigned int scale = 0;
        unsigned int ID = 0;
        // unsigned int interpolations;
    };

    struct RT_Grouping
    {
        size_t originalIndex = -1;
        size_t interpolatedIndex = -1;
        float trueRT = -1;
        bool interpolated = true;
    };

    struct RT_Converter
    {
        std::vector<RT_Grouping> groups;
        // index into the groups vector. The "originalIndex" field ind the RT_Grouping struct is the index into this vector
        std::vector<size_t> indexOfOriginalInInterpolated = {0};
    };

    struct EIC // Extracted Ion Chromatogram
    {
        std::vector<unsigned int> scanNumbers = {0};
        std::vector<float> mz = {0};
        std::vector<float> predInterval{0};
        std::vector<float> ints_area{0};
        std::vector<float> ints_height{0};
        std::vector<unsigned int> df{0}; // this is required for dealing with interpolations, but should be moved into qPeaks eventually @todo
        std::vector<float> DQSB{0};
        std::vector<float> DQSC{0};
        std::vector<unsigned int> cenID{0};
        std::vector<unsigned int> interpolatedScans{0};
        std::vector<float> RT{0};
        // std::vector<float> interpolatedDQSB;
        unsigned int componentID = 0; // this is only set during componentisation
        // bool interpolations;          // are there interpolated values?
    };

    struct FeaturePeak
    {
        RegCoeffs coefficients{0};
        float height = 0;
        float area = 0;
        // float width;
        float heightUncertainty = 0;
        float areaUncertainty = 0;
        float DQSF = 0, DQSB = 0, DQSC = 0;
        float retentionTime = 0;
        float mz = 0;
        float RT_Uncertainty = 0;
        float mzUncertainty = 0;
        unsigned int componentID = 0; // this is only set during execution of qPattern / qComponent / whatever better name i think of. Zero means uninitialised -> components start at 1!
        unsigned int idxBin = 0;
        // these refer to the interpolated EIC!
        unsigned int idxPeakStart = 0, idxPeakEnd = 0, idxCenter_offset = 0;
        // relates to abstracted MS1 scan counts, starts at 2 for real points
        unsigned int scanPeakStart = 0, scanPeakEnd = 0;
        // indices into the non-interpolated bin; degrees of freedom = idxBinEnd - idxBinStart + 1
        unsigned int idxBinStart = 0, idxBinEnd = 0;
        // temporary values, @todo remove?
        unsigned int interpolationCount = 0;
        unsigned int competitorCount = 0;
        unsigned int scale = 0;
        float mse_base = 0;
        float lowerRT = 0;
        float upperRT = 0; // @todo set these during feature construction
    };
}

#endif // QALGORITHMS_DATATYPE_PEAK_H
