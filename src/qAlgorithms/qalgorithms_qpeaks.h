// qalgorithms_qpeaks.h
#ifndef QALGORITHMS_QPEAKS_H
#define QALGORITHMS_QPEAKS_H

// internal
#include "./qalgorithms_datatypes.h"
// #include "qalgorithms_read_file.h" // @todo get rid of this coupling!

#include <vector>

namespace qAlgorithms
{
    // @todo rework this module so all centroiding / feature construction parts are in a separate block and qPeaks can be used as a library header

    // ### General Implementation ### //

    /// @brief This function serves as the central interface for performing peak detection on an arbitrary set of x and y data.
    /// @param y_values Values used to determine height of the peak
    /// @param x_values Values used as width of the peak, same length as y_values. It is a requirement that x values are evenly
    ///                 spaced when applying this function, since it depends on that assumption for its precalculated matric transform.
    /// @param degreesOfFreedom_cum cumulative degrees of freedom to account for interpolations and mean values. This is the cumulative
    ///                             sum of the number of measured points for all y. It has the same length as x and y.
    ///                             This can be set to nullptr if all df are one.
    /// @param length length of y, x and df if supplied
    /// @param maxScale_in up to which scale regressions should be attempted. The maximum peak width (in points) is maxScale * 2 + 1. maxScale must be > 1.
    ///                 If maxScale would exceed the length of x or y, it is set to (the length of either -1) / 2.
    /// @param result All found peaks are appended to this vector.
    /// @return the total number of appended peaks. If this is < 0, the function did not execute correctly.
    /// @details Errorcodes:
    ///  0 = no valid peaks were found (this is not an error per se, when supplying ex. a constant y this is the expected result)
    /// -1 = one of y_values, x_values or detectedPeaks was nullptr
    /// -2 = y_values, x_values or degreesOfFreedom_cum have unequal lengths
    /// -3 maxScale is < 2
    /// ### Assumptions: ###  these are not tested @todo
    /// x and degreesOfFreedom_cum increase monotonically.
    /// y has equal variance at every point
    /// there is enough space to write the results
    int qpeaks_find(
        const float *y_values,
        const float *x_values,
        const unsigned int *degreesOfFreedom_cum,
        const size_t length,
        const size_t maxScale_in,
        std::vector<PeakFit> *result);

    /// @param intensity_log logarithmy of the intensity values. This is the y axis of the fit. The x axis is required to be equidistant.
    /// @param maxscale maximum scale of a peak that should be attempted to fit.
    /// @param coeffs Sets of coefficients for all possible regressions. They are written out in the order scale = 2, scale = 3, ... , scale = maxscale
    ///               Note that not all scales have the same number of regressions
    void findCoefficients(
        const float *intensity_log,
        const size_t length,
        size_t maxscale,
        std::vector<RegCoeffs> *coeffs);

    void findBestScales(std::vector<RegressionGauss> *validRegressions,
                        std::vector<RegressionGauss> *validRegsTmp,
                        const float *intensities,
                        const unsigned int *const degreesOfFreedom_cum);

    void runningRegression(
        const float *intensities,
        std::vector<float> *intensities_log,
        const unsigned int *const degreesOfFreedom_cum,
        const size_t length,
        const size_t maxScale,
        std::vector<RegressionGauss> *validRegressions);

    // mutate b0 so that it is optimal for the exponential case if b1, b2 and b3 are identical

    /// @brief adjust the height of a regression to better fit the exponential data
    /// @param intensities non-logarithmic intensity values the regression was fitted to
    /// @param r range of the regression
    /// @param predicted empty vector that the predicted values for intensity (AFTER correction) are written to
    /// @param coeff coefficients that should be updated
    /// @return used correction factor
    double correctB0(const float *const intensities,
                     const Range_i *r,
                     std::vector<float> *predicted,
                     RegCoeffs *coeff);

    enum invalid
    {
        ok,
        no_apex,
        invalid_apex,
        no_df,
        invalid_apexToEdge, // this will probably be removed
        f_test_fail,
        invalid_quadratic,
        invalid_area,
        invalid_height,
        invalid_chisq,
        none = -1
    };

    /// @brief perform various statistical tests to see if a regression describes a valid peak
    /// @param degreesOfFreedom_cum cumulative degrees of freedom (only relevant for interpolated data)
    /// @param intensities measured intensities
    /// @param intensities_log log of measured intensities - must have same length as intensities
    /// @param mutateReg regression that should be mutated by this function
    /// @return 0 if the regression is valid, otherwise the filter step which kicked it out
    invalid makeValidRegression(
        const float *intensities,
        const std::vector<float> *intensities_log,
        const std::vector<float> *predict,
        const size_t df_sum,
        const size_t length,
        RegressionGauss *mutateReg);

    invalid validRegWidth(const RegCoeffs *coeffs, Range_i *range);

    // ### Centroiding-specific Code ### //

    // int findCentroids(XML_File &data, // @todo get rid of the direct coupling to pugixml
    //                   const std::vector<unsigned int> *selectedIndices,
    //                   std::vector<CentroidPeak> *centroids);

    void findCentroidPeaks(std::vector<CentroidPeak> *retPeaks, // results are appended to this vector
                           const std::vector<ProfileBlock> *subprofiles,
                           const size_t scanNumber,
                           const size_t ID_spectrum,
                           const size_t maxWindowSize);

    void createCentroidPeaks(
        std::vector<CentroidPeak> *peaks,
        const std::vector<RegressionGauss> *validRegressionsVec,
        const ProfileBlock *block,
        const size_t scanNumber,
        const size_t ID_spectrum);

    // ### Feature-specific Code ### //

    std::vector<FeaturePeak> findFeatures(std::vector<EIC> &data,
                                          const RT_Converter *convertRT);

    // ### Retention Time Conversion @todo remove ### //

    // RT_Converter interpolateScanNumbers(const std::vector<float> *retentionTimes);

    void createFeaturePeaks(
        std::vector<FeaturePeak> *peaks,
        const std::vector<RegressionGauss> *validRegressionsVec,
        const RT_Converter *convertRT,
        const std::vector<float> *RTs);

    /// @brief calculate the residual sum of squares for the log regression / data
    /// @param mutateReg relevant regression
    /// @param y_start log data
    /// @return RSS value
    double calcRSS_log(const RegressionGauss *mutateReg, const std::vector<float> *y_start);

    /// @brief performs two F-tests against the log data. First H0 is the mean, second y = mx + b
    /// @param observed log data (or normal data, depends on the use case)
    /// @param RSS_reg previously calculated residual sum of squares of the complex model. Hard assumpion of four coefficients.
    /// @param range range of the regression.
    /// @return true: Regression is significant; false: Regression is not better than either alternative.
    bool f_testRegression(const float *observed, double RSS_reg, const Range_i *range);

    double calcMSE_exp(const RegCoeffs *coeff,
                       const float *observed,
                       const Range_i *regSpan,
                       const double df);

    double calcSSE_chisqared(const RegressionGauss *mutateReg,
                             const float *observed,
                             const std::vector<float> *predict);

    struct RegPair
    {
        unsigned int idx;
        double mse;
    };

    RegPair findBestRegression(
        const float *intensities,
        const std::vector<RegressionGauss> *regressions,
        const unsigned int *const degreesOfFreedom_cum,
        const Range_i regSpan);

    /**
     * @brief updates apex position field of mutateReg and the supplied valley_position
     * @return 0 if positions could be calculated, otherwise the point of failure in the
     * function is returned.
     */
    int calcApexAndValleyPos(
        RegressionGauss *mutateReg,
        double *valley_position);

    // take a jacobian matrix as input and return the transpose at scale
    double matProductReg(const double J[4], const size_t scale);

    bool isValidQuadraticTerm(const RegCoeffs *coeffs, const double mse, const size_t df_sum);

    bool isValidPeakHeight(
        const RegressionGauss *mutateReg,
        const double valley_position,
        const size_t df_sum,
        const double apexToEdge);

    void calcPeakHeightUncert(RegressionGauss *mutateReg);

    /**
     * @brief Check if the peak area and the covered peak area are valid using t-test.
     * @details The function calculates the peak area and the covered peak area using the
     * regression coefficients. The peak area is the integral of the regression model
     * from -infinity to +infinity. The covered peak area is the integral of the regression
     * model from the left limit of the regression window to the right limit of the
     * regression window. Moreover, the trapzoid under the peak is considered as
     * not covered peak area.
     *
     * @param coeff : Matrix of regression coefficients
     * @param mse : mean squared error
     * @param scale : Window size scale, e.g., 5 means the window size is 11 (2*5+1)
     * @param df_sum : sum of the degree of freedom of the regression model
     * @param area : area of the peak
     * @param area_uncert : uncertainty of the area
     * @return true : if the peak area is valid
     * @return false : if the peak area is not valid
     */

    void calcPeakAreaUncert(RegressionGauss *mutateReg);

    bool isValidPeakArea(const RegCoeffs *coeffs, const double mse, const size_t df_sum);

    void calcUncertaintyPos(RegressionGauss *mutateReg);

    struct MeanVar
    {
        float mean;
        float var;
    };

    MeanVar weightedMeanAndVariance_EIC(const std::vector<float> *weight,
                                        const std::vector<float> *values,
                                        const Range_i regSpan);

    // utility functions for calculating regression values
    double regAt(const RegCoeffs *coeff, const double x);

    double fullWidthHalfMax(const RegCoeffs *coeff, const double height, const double delta_x);
    double FWHM_to_sdev(const double fwhm);

    // this one does not include b0
    double regExp_fac(const RegCoeffs *coeff, const double x);

    double medianVec(const std::vector<float> *vec);

    // double peakAreaFull(const RegCoeffs *coeff, const double delta_x);

    double peakArea(const double b0, const double b1, const double b2, const double b3, const double delta_x);

    // ### pre-calculate the regression matrix ### //
#define MAXSCALE 63
#define GLOBAL_MAXSCALE_CENTROID 8 // @todo this is a critical part of the algorithm and should not be hard-coded
#define GLOBAL_MINSCALE 2

#include "./qalgorithms_matinverse.h"

    ///     This function performs a convolution with the kernel (xTx)^-1 xT and the data array intensity_log.

    ///     (xTx)^-1 is pre-calculated and stored in the vector INV_ARRAY (calculated in the "initialise" function).
    ///     Only six values of the final matrix are required for the simple case, see below:

    ///     xT is the transpose of the design matrix X.
    ///     for scale = 2:
    ///     xT = | 1  1  1  1  1 |    : all ones
    ///          |-2 -1  0  1  2 |    : from -scale to scale
    ///          | 4  1  0  0  0 |    : x^2 values for x < 0
    ///          | 0  0  0  1  4 |    : x^2 values for x > 0

    ///     It contains one additional row of all ones for every additional peak that is added into the model

    ///     When adding multiple peaks to the regression model, we need to adjust the inverse values.
    ///     This will change the number of unique values in the inv_values array from 6 to 7.
    ///     Here we use the inv_array[1] position and shift all values from that point onwards to the right.
    ///     example for num_peaks = 2:
    ///     original matrix with the unique values [a, b, c, d, e, f] (six unique values)
    ///     | a  0  b  b |
    ///     | 0  c  d -d |
    ///     | b  d  e  f |
    ///     | b -d  f  e |

    ///     new matrix with the unique values [A1, A2, B, C, D, E, F] (seven unique values)
    ///     | A1  A2  0  B  B |
    ///     | A2  A1  0  B  B |
    ///     | 0   0   C  D -D |
    ///     | B   B   D  E  F |
    ///     | B   B  -D  F  E |

    ///     for num_peaks = 3:
    ///     new matrix with the unique values [A1, A2, B, C, D, E, F] (the same seven unique values)
    ///     | A1  A2  A2  0  B  B |
    ///     | A2  A1  A2  0  B  B |
    ///     | A2  A2  A1  0  B  B |
    ///     | 0   0   0   C  D -D |
    ///     | B   B   B   D  E  F |
    ///     | B   B   B  -D  F  E |

    ///     Note that no more than seven different values are needed per scale, even for a multidimensional approach.
}

#endif