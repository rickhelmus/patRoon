#include "qalgorithms_qpeaks.h"
#include "qalgorithms_utils.h"
#include "qalgorithms_datatypes.h"
// #include "qalgorithms_read_file.h" // @todo remove coupling

#include "liberfc_reduced.h"

#include <cassert>
#define _USE_MATH_DEFINES // relevant for windows to have math constants
#include <math.h>
#include <vector>

#include <algorithm> // sorting @todo get rid of this
// #include <stdint.h> // max numeric vals

namespace qAlgorithms
{
    struct LoggingParams
    {
        bool features = false;
        int dataID = 0;
        int checkedRegionCount = 0;
        int totalPeakCount = 0;
        int failData[invalid::invalid_chisq + 1] = {0};
    };

    void loggerPrint(LoggingParams *log)
    {
        // print log to hardcoded logfile path. Example output:
        // C, 1000000, 150000, 1234567, 1234567, 1234567, 1234567, 1234567, 1234567, 1234567, 1234567, 1234567, 1234567 \n
        // -> 300 chars space will be more than sufficient

        // header: type, regions, peaks, passed, no_apex, invalid_apex, no_df, invalid_apexToEdge, f_test_fail, invalid_quadratic, invalid_area, invalid_height, invalid_chisq
        char buffer[300];
        sprintf(buffer, "%c, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
                log->features ? 'F' : 'C', log->checkedRegionCount, log->totalPeakCount, log->failData[invalid::ok],
                log->failData[invalid::no_apex], log->failData[invalid::invalid_apex], log->failData[invalid::no_df],
                log->failData[invalid::invalid_apexToEdge], log->failData[invalid::f_test_fail], log->failData[invalid::invalid_quadratic],
                log->failData[invalid::invalid_area], log->failData[invalid::invalid_height], log->failData[invalid::invalid_chisq]);

        FILE *file = fopen("./logdata.csv", "a");
        assert(file);
        fprintf(file, "%s", buffer);
        fclose(file);
    }

    void resetLogger(LoggingParams *log)
    {
        log->checkedRegionCount = 0;
        log->totalPeakCount = 0;

        size_t len = invalid::invalid_chisq + 1;
        for (size_t i = 0; i < len; i++)
        {
            log->failData[i] = 0;
        }
    }

    LoggingParams globalLogStruct;

    void retransformPeaks(
        const std::vector<RegressionGauss> *peaks,
        const float *x_values,
        const size_t length,
        std::vector<PeakFit> *result)
    {
        // retransforming is probably best if the range of x values is as small as necessary
        for (size_t i = 0; i < peaks->size(); i++)
        {
            const RegressionGauss *regression = peaks->data() + i;
            assert(regression->isValid);
            const RegCoeffs coeff = regression->coeffs;

            double apex = regression->apex_position;
            size_t start = regression->regSpan.startIdx;
            size_t end = regression->regSpan.endIdx;

            PeakFit peak;

            // @todo ensure that this is a good way to estimate delta_x for real data
            // we could also use a strategy such as taking the distance closest to the apex.

            // no +1 for the length is used here because there are n-1 distances for n points
            double delta_x = (x_values[end] - x_values[start]) / (end - start);

            // position is determined relative to the point left of the apex
            size_t leftOfApex = size_t(apex);
            double apexFraction = delta_x * (apex - trunc(apex));
            peak.position = x_values[leftOfApex] + apexFraction;
            peak.position_uncert = regression->position_uncert * delta_x;

            // check if the chosen delta_x is acceptable
            double delta_x_2 = x_values[leftOfApex + 1] - x_values[leftOfApex];
            assert(abs(delta_x - delta_x_2) < delta_x * 0.05);

            // peak height (exp(b0 - b1^2/4/b2)) with position being -b1/2/b2
            peak.height = exp(coeff.b0 + (apex - coeff.x0) * coeff.b1 * 0.5);
            peak.height_uncert = regression->height_uncert * peak.height;

            double factorArea = exp(coeff.b0) * delta_x;
            peak.area = regression->area * factorArea;
            peak.area_uncert = regression->area_uncert * factorArea;

            // @todo, also check for negative width where appropriate
            // the empirical peak width is generally estimated at half maximum. Our peak
            // model only has a standard deviation for the apex peak
            peak.fwhm = fullWidthHalfMax(&coeff, peak.height, delta_x);

            peak.dqs = 1 - erf_approx_f(regression->area_uncert / regression->area);
            peak.jaccard = regression->jaccard;

            peak.coeffs = regression->coeffs;

            result->push_back(peak);
        }
    }

    int qpeaks_find(
        // SOME_IMPLEMENTATION_OF_LINEAR_ALLOCATOR_HERE
        const float *y_values,
        const float *x_values,
        const unsigned int *degreesOfFreedom_cum,
        const size_t length,
        const size_t maxScale_in,
        std::vector<PeakFit> *result)
    {
        // control input for nullpointers, mismatching x and y, and fitting maxScale
        // @todo use one enum for errors uniformly throughout the code
        if (y_values == nullptr || x_values == nullptr || result == nullptr)
        {
            return -1;
        }
        // if (y_values->size() != x_values->size())
        // {
        //     return -2;
        // }
        if (length < 5)
        {
            return -3;
        }
        if (maxScale_in < GLOBAL_MINSCALE)
        {
            return -4;
        }

        /*
        The first, implicit operation is norming x so that every value of x is at a distance
        of exactly 1 to both of its neighbouring elements. This allows us to use a pre-calculated
        design matrix for every operation, which in turn massively reduces computation time.

        Assuming that the points are equidistant in x with delta(x), the axis used for regression
        is x' = x / delta(x). Additionally, x0 is different for every individual regression. Both
        of these transformations must be reversed after finding the best fit, even if at no point
        in the program there is a function that actually performs either transform explicitly.

        At least for HRMS data, the majority of data regions does not contain any peaks (both
        for centroids and features). To avoid unnecessary computations of delta_x, it is only
        determined for cases where at least one regression was found.
        */

        // @todo the assumtion that all values of x are equidistant is conrolled. If gaps are
        // found, the missing values are interpolated assuming an exponential rate of change.
        // this should happen before calling this function (?)

        // core operation: identify best-fit regressions for the input data

        size_t maxScale = (length - 1) / 2; // apex is not included in scale -> -1
        maxScale = maxScale_in > maxScale ? maxScale : maxScale_in;

        /*
        The fitting routine assumes that all present peaks have a modified gaussian base function.
        This means that no baseline exists. Baseline substraction, if appropriate, has to be performed
        before calling the qpeaks_find.
        */
        std::vector<float> y_log;

        std::vector<RegressionGauss> validRegressions;
        runningRegression(
            y_values,
            &y_log,
            degreesOfFreedom_cum,
            length,
            maxScale,
            &validRegressions);

        // validregressions contains all relevant peak candidates

        // reverse the not-calculated transform of the x axis only for region where peaks exist
        retransformPeaks(&validRegressions, x_values, length, result);

        return -42;
    }

#pragma region "find peaks"

    void findCentroidPeaks(std::vector<CentroidPeak> *retPeaks, // results are appended to this vector
                           const std::vector<ProfileBlock> *subprofiles,
                           const size_t scanNumber,
                           const size_t ID_spectrum,
                           const size_t maxWindowSize)

    {
        assert(!subprofiles->empty());

        std::vector<float> logIntensity(maxWindowSize + 2 * GLOBAL_MAXSCALE_CENTROID, NAN);

        static_assert(GLOBAL_MAXSCALE_CENTROID <= MAXSCALE);

        std::vector<RegressionGauss> validRegressions;
        validRegressions.reserve(subprofiles->size() / 2); // probably too large, shouldn't matter
        for (size_t i = 0; i < subprofiles->size(); i++)
        {
            const ProfileBlock *const block = subprofiles->data() + i;

            // this is now filled inside the function, the vector only reserves space. We do not
            // perform this step in the function so that it is explicitly empty. This should be
            // replaced by a non-malloc calling scratch space eventually
            logIntensity.clear();
            runningRegression(
                block->intensity,
                &logIntensity,
                nullptr,
                block->length,
                GLOBAL_MAXSCALE_CENTROID,
                &validRegressions);

            volatile bool purge = false;
            if (purge)
                validRegressions.clear();

            if (validRegressions.empty())
            {
                continue; // no valid peaks
            }
            createCentroidPeaks(retPeaks, &validRegressions, block, scanNumber, ID_spectrum);
            validRegressions.clear();
        }
    }

#pragma endregion "find peaks"

#pragma region "running regression"

    volatile bool debug = false;
    std::vector<invalid> failCodes; // @todo: make a logging function / struct thing that can be used to check failure points later

    int validateRegressions( // @todo this should be specific to centroids, features or components
        const float *intensities,
        const std::vector<float> *intensities_log,
        const unsigned int *const degreesOfFreedom_cum,
        const std::vector<RegCoeffs> *coefficients,
        const size_t length,
        std::vector<RegressionGauss> *validRegressions)
    {
        size_t validCount = 0;
        for (size_t i = 0; i < coefficients->size(); i++)
        {
            Range_i range;
            invalid failpoint = validRegWidth(&(coefficients->at(i)), &range);

            if (failpoint != ok)
            {
                globalLogStruct.failData[failpoint] += 1;
                continue;
            }

            size_t df_sum = sumOfCumulative(degreesOfFreedom_cum, &range);
            if (df_sum < 5)
            {
                globalLogStruct.failData[invalid::no_df] += 1;
                continue;
            }
            df_sum -= 4; // four coefficients, adjust for components

            RegressionGauss reg;
            reg.coeffs = coefficients->at(i);
            reg.regSpan = range;

            /*
           Prediction for coefficients that are not b0. Since any "true" prediction
           value is produced by multiplying with exp(b0), this only needs to be calculated
           once, at least until regressions are compared. @todo consider having this at
           higher scope of validation
       */
            std::vector<float> predict(length, 0);

            /*
                Adjustment of b0 coefficient:
                When working with log-transformed data, the coefficients are suboptimal for the exponential case.
                Since we must work with a log system to perform a linear regression, there is a bias in the
                results which is somewhat corrected here. While correction could occur before validation, the
                initial tests filter out a lot of bad regressions which reduces processing time. The tests are
                presumed to be better when using the transformed coefficients in terms of applicability of the results.
                This function also modifies the "predict" vector supplied as its argument!
            */
            correctB0(intensities, &range, &predict, &reg.coeffs);

            failpoint = makeValidRegression(intensities,
                                            intensities_log,
                                            &predict,
                                            df_sum,
                                            length,
                                            &reg);
            failCodes.push_back(failpoint);
            globalLogStruct.failData[failpoint] += 1;

            validCount += failpoint == ok ? 1 : 0;
            if (failpoint == ok)
            {
                validRegressions->push_back(reg);
            }
        }
        assert(validCount == validRegressions->size());

        failCodes.push_back(invalid::none);

        return validCount;
    }

    int pruneConflictingRegs(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const unsigned int *const df_cum);

    int resolveScaleConflicts(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const unsigned int *const df_cum);

    int pruneRegsByApex(const float *intensities,
                        const unsigned int *const df_cum,
                        std::vector<RegressionGauss> *validRegressions);

    // ----------- old functions -------------//

    void mergeRegsInScale(
        const float *intensities,
        const unsigned int *const df_cum,
        std::vector<RegressionGauss> *validRegsTmp,
        std::vector<RegressionGauss> *validRegressions);

    void mergeRegressionsOverScales(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities);

    // -------------------------------------- //

    void runningRegression(
        const float *intensities,
        std::vector<float> *intensities_log,
        const unsigned int *const degreesOfFreedom_cum,
        const size_t length,
        const size_t maxScale,
        std::vector<RegressionGauss> *validRegressions)
    {
        assert(validRegressions->empty());

        intensities_log->clear();
        intensities_log->reserve(length);
        for (size_t i = 0; i < length; i++)
        {
            // assert(intensities[i] > 1); // @todo this could be relevant, currently ensured by the region pre-filter
            intensities_log->push_back(log(intensities[i]));
        }

        // coefficients for single-b0 peaks, spans all regressions over a peak window
        // all entries in coeff are sorted by scale and position in ascending order - this is not checked!
        std::vector<RegCoeffs> coefficients;
        findCoefficients(intensities_log->data(), length, maxScale, &coefficients);

        std::vector<RegressionGauss> validRegsTmp; // all independently valid regressions regressions
        validRegsTmp.reserve(coefficients.size() / 2);
        int validCount = validateRegressions(intensities,
                                             intensities_log,
                                             degreesOfFreedom_cum,
                                             &coefficients,
                                             length,
                                             &validRegsTmp);
        assert(validCount <= int(coefficients.size()));
        if (validCount == 0)
            return;

        int invalidCount = pruneConflictingRegs(&validRegsTmp, intensities, degreesOfFreedom_cum);

        assert(invalidCount < validCount);
        // if (validCount - invalidCount == 1) // terminate early if only one regression remains
        //     return; // @todo this vastly changes results, check for correctness!

        resolveScaleConflicts(&validRegsTmp, intensities, degreesOfFreedom_cum);

        int apexCount = pruneRegsByApex(intensities, degreesOfFreedom_cum, &validRegsTmp);

        // temp storage vector, eventually delete everything below this comment @todo
        std::vector<RegressionGauss> validRegsTmp2;
        for (size_t i = 0; i < validRegsTmp.size(); i++)
        {
            RegressionGauss *reg = validRegsTmp.data() + i;
            if (reg->isValid)
                validRegsTmp2.push_back(*reg);
            else
                reg->isValid = true; // required to force compliance with the old system for now
        }

        mergeRegsInScale(
            intensities,
            degreesOfFreedom_cum,
            &validRegsTmp,
            validRegressions);

        // there can be 0, 1 or more than one regressions in validRegressions
        mergeRegressionsOverScales(validRegressions, intensities);

        if (validRegsTmp2.size() != validRegressions->size())
        {
            for (size_t i = 0; i < length; i++)
            {
                printf("%f, ", intensities[i]);
            }
            printf("\n");

            // exit(1);
        }

        // assert(validRegressions->size() == apexCount);
        // assert(validRegsTmp2.size() == validRegressions->size());
        return;
    }

#pragma region "eliminate conflicting regs"

    int pruneRegsByApex(const float *intensities,
                        const unsigned int *const df_cum,
                        std::vector<RegressionGauss> *validRegressions)
    {
        // yet another attempt to solve this problem gracefully. This time, the
        // strategy is to pick the first apex described by a regression and, regardless
        // of scale, to select all regressions which describe the same peak. Then, the
        // best among them is chosen iteratively by testing from the largest scale downward.

        // 0: count the number of probable apexes
        std::vector<float> apexes;
        for (size_t i = 0; i < validRegressions->size(); i++)
            apexes.push_back(validRegressions->at(i).apex_position);

        std::sort(apexes.begin(), apexes.end());

        int apexcount = 1;
        for (size_t i = 1; i < apexes.size(); i++)
        {
            if (apexes[i] - apexes[i - 1] > 2)
                apexcount += 1;
        }

        return apexcount;

        // step one: for the first regression, select all other regressions that describe the same apex
        std::vector<size_t> conflictIdx;
        conflictIdx.reserve(validRegressions->size());

        RegressionGauss *startReg = validRegressions->data();
        double startApex_L = startReg->apex_position;
        double startApex_R = startReg->apex_position;
        for (size_t i = 0; i < validRegressions->size(); i++)
        {
            RegressionGauss *reg = validRegressions->data() + i;
            double apex = reg->apex_position;
            double diff = apex - startApex_L;

            if (abs(diff) > 2 * GLOBAL_MINSCALE)
                continue;
            diff = apex - startApex_R;
            if (abs(diff) > 2 * GLOBAL_MINSCALE)
                continue;

            // this should generally cover all cases due to the initial smaller scale being further left and right than later regressions
            startApex_L = min(apex, startApex_L);
            startApex_R = max(apex, startApex_R);

            // the regressions are in conflict
            conflictIdx.push_back(i);
        }

        // the main issue with multiple regressions is estimating if the greater scale ones cover too many points.
        // if not, should overlap between regressions be permitted?
        return 0;
    }

    double calcMSE_exp(const RegCoeffs *coeff,
                       const float *observed,
                       const Range_i *regSpan,
                       const double df)
    {
        double idxCenter = double(coeff->x0);
        double x = double(regSpan->startIdx) - idxCenter;
        double result = 0.0;
        for (size_t i = regSpan->startIdx; i < regSpan->endIdx + 1; i++)
        { // this loop is vectorised at level O3
            double pred = exp(regAt(coeff, x));
            double obs = observed[i];
            double newdiff = (obs - pred) * (obs - pred);
            result += newdiff;
            x += 1.0;
        }
        double mse = result / df;
        assert(mse > 0);
        return mse;
    }

    int pruneConflictingRegs(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const unsigned int *const df_cum)
    {
        // regressions start out valid and sorted by scale, and position within scales
        // S1R1 S1R2 S1 R3 S1 Rn | S2R1 S2R2 S2Rn | SnRn

        // case 0: SnRn is the only element with scale Sn
        //    -->  it stays valid, continue

        // case 1: SnR1 and SnR2 do not conflict
        //    -->  SnR1 stays valid and testing continues with SnR2 and SnR3

        // case 2: SnRn is the last regression in its scale
        //    -->  it stays valid / invalid and we move on to Sn+1R1

        // case 3: SnRn is invalid
        //    -->  continue

        // case 4: SnR1, SnR2 and SnR2, SnR3 conflict, but SnR1, SnR3 does not
        //    -->  this is assumed to be extremely unlikely at one scale and also not covered for by the existing solution
        // (this is only relevant if SnR1 and SnR3 are better than SnR2 when combined, which might not even be possible)
        // the main relaistic problem case is all three being interesting but conflicting due to baseline overlap

        size_t eliminations = 0;

        if (validRegressions->size() < 2) // only run function if there can be conflict
            return eliminations;

        RegressionGauss *compReg = validRegressions->data();
        assert(compReg->numCompetitors == 0);
        for (size_t i = 1; i < validRegressions->size(); i++)
        {
            RegressionGauss *reg = validRegressions->data() + i;
            assert(reg->isValid);
            size_t currentScale = reg->coeffs.scale;

            if (currentScale != compReg->coeffs.scale)
            {
                compReg = reg;
                continue;
            }
            if (compReg->regSpan.endIdx <= reg->regSpan.startIdx) // the regressions can switch at one point
            {
                // no conflict despite same scale
                compReg = reg;
                continue;
            }

            // both regressions conflict, one is invalid. The better regression is determined
            // by comparing MSEs for measured values over the region of both regressions
            // @todo could it be better to decide based on something else, since the range
            // is not adjusted afterward?
            Range_i newRange = {compReg->regSpan.startIdx, reg->regSpan.endIdx};
            size_t df = sumOfCumulative(df_cum, &newRange) - 4;

            double mse_reg = calcMSE_exp(&reg->coeffs,
                                         intensities,
                                         &newRange,
                                         df);

            double mse_prev = calcMSE_exp(&compReg->coeffs,
                                          intensities,
                                          &newRange,
                                          df);

            if (mse_reg < mse_prev)
            {
                // the current regression is better
                compReg->isValid = false;
                reg->numCompetitors += compReg->numCompetitors + 1;
                compReg = reg;
            }
            else
            {
                // the old regression is better
                compReg->numCompetitors += reg->numCompetitors + 1;
                reg->isValid = false;
            }
            assert(compReg->numCompetitors < int(validRegressions->size())); // at most better than all other regressions
            eliminations += 1;
        }

        // sanity check: the number of competitors in valid regressions is the number of invalid regressions
        size_t invalidCount = 0;
        size_t competitors = 0;
        for (size_t i = 0; i < validRegressions->size(); i++)
        {
            RegressionGauss *reg = validRegressions->data() + i;
            if (reg->isValid)
                competitors += reg->numCompetitors;
            else
                invalidCount += 1;
        }
        assert(invalidCount == eliminations);
        assert(invalidCount == competitors);

        return eliminations; // if only one regression remains, the next function can be skipped!
    }

    // these are helper functions for the big resolveScaleConflicts part of the program

    int findNextReg(
        const std::vector<RegressionGauss> *validRegressions,
        const size_t scale,
        const size_t startIdx)
    {
        // find the first valid regression that has the required scale at idx >= startidx or -1 if search fails
        int pos = int(startIdx);
        int len = validRegressions->size();
        bool found = false;
        for (int dpos = pos; dpos < len; dpos++)
        {
            const RegressionGauss *reg = validRegressions->data() + dpos;
            if (reg->coeffs.scale < scale)
                continue;
            if (!reg->isValid)
                continue;

            found = true;
            pos = dpos;
            break;
        }
        if (!found)
            return -1;

        return pos;
    }

    Range_i findComparisonRegs(
        const std::vector<RegressionGauss> *validRegressions,
        const Range_i *range,
        const size_t scale)
    {
        // since all regressions are in order, only the outermost two must be passed along
        // if regs = 1, 0: no valid regs found at this scale / in this range
        Range_i regs = {1, 0};
        bool first = true;

        for (size_t i = 0; i < validRegressions->size(); i++)
        {
            const RegressionGauss *reg = validRegressions->data() + i;
            // only consider regressions in the current scale
            if (reg->coeffs.scale < scale)
                continue;
            if (reg->coeffs.scale > scale)
                break;

            // only apex overlap is checked, good idea? @todo
            // this gave bad results where an intermediate regression was left despite conflicts
            // size_t apex = size_t(reg->apex_position) + 1;
            // if (apex < range->startIdx)
            //     continue;
            // if (apex > range->endIdx)
            //     break;
            // at most one point overlap is tolerated since the regression "ends" at that point
            if (reg->regSpan.endIdx <= range->startIdx)
                continue;
            if (reg->regSpan.startIdx >= range->endIdx)
                break;

            // validity check here due to early termination from previous checks
            if (!reg->isValid)
                continue;

            if (first)
            {
                first = false;
                regs.startIdx = i;
            }
            regs.endIdx = i;
        }

        return regs;
    }

    size_t invalidateRange(const Range_i *r, std::vector<RegressionGauss> *validRegressions)
    {
        // returns the total competitors eliminated
        size_t competitors = 0;
        for (size_t i = r->startIdx; i < r->endIdx; i++)
        {
            RegressionGauss *reg = validRegressions->data() + i;
            competitors += reg->isValid ? reg->numCompetitors + 1 : 0;
            reg->isValid = false;
        }
        return competitors;
    }

    void addCompetitors(const Range_i *r, std::vector<RegressionGauss> *validRegressions, size_t comp)
    {
        // distribute the number in comp evenly. Bias towards lower end regressions is accepted
        while (comp != 0)
        {
            for (size_t i = r->startIdx; i < r->endIdx; i++)
            {
                RegressionGauss *reg = validRegressions->data() + i;
                if (reg->isValid)
                {
                    reg->numCompetitors += 1;
                    comp -= 1;
                }
            }
        }
    }

    double multiMSE(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const Range_i *selectRegs,
        const Range_i *commonRange,
        size_t df)
    {
        // the df must be without coefficients! (10 = at least two scale 2 regressions)
        assert(df >= 10);

        double summedMSE = 0;
        size_t regCount = 1; // first one is not included in the loop
        size_t currentStart = commonRange->startIdx;
        size_t currentEnd = 0;

        RegressionGauss *regFront = validRegressions->data() + selectRegs->startIdx;

        for (size_t i = selectRegs->startIdx + 1; i < selectRegs->endIdx + 1; i++)
        {
            RegressionGauss *regCurr = validRegressions->data() + i;
            if (!regCurr->isValid)
                continue;

            // determine MSE for the left of both regressions
            // the start is known from the previous regression, the end is the middle between both regs
            currentEnd = (regFront->regSpan.endIdx + regCurr->regSpan.startIdx) / 2;
            Range_i rangeNow = {currentStart, currentEnd};

            // the mse is composed of individual regressions that do not share responsibility
            // for the same point
            summedMSE += calcMSE_exp(&regFront->coeffs, intensities, &rangeNow, df);

            // set current regression as the active element of the mse
            regFront = regCurr;
            currentStart = currentEnd + 1;
            regCount += 1;
        }

        // there is still one unprocessed regression
        currentEnd = commonRange->endIdx;
        Range_i rangeNow = {currentStart, currentEnd};
        summedMSE += calcMSE_exp(&regFront->coeffs, intensities, &rangeNow, df);

        return summedMSE / (df - regCount * 4); // important: two individual regs lose degrees of freedom
    }

    int resolveScaleConflicts(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities,
        const unsigned int *const df_cum)
    {
        // this function takes in the output of the previous function and then checks which regressions
        // should be preferred. It is important that comparisons are made in ascending scale order so that
        // overfitting due to comparing two too broad regressions is minimised

        // case 0: only one valid regression
        //    -->  do nothing

        // case 1: one valid regression at lower scale conflicts with a valid regression at a greater scale
        //    -->  pick the better one for the full range of both

        // case 2: Multiple regressions at lower scale conflict with the same regression at a greater scale
        //    -->  switch contribution to the mse at the mean of both limits

        // this function assumes that a conflict exists
        // there is never conflict within a scale! (previous function execution)

        // algorithm:
        // 0: start at scale = 3
        // is there a regression at this scale? -> function
        // F: scale += 1, go to 0
        // T: set this scale to the origin scale

        // compare at scale = 2
        // 1: Iterate from i = 0 over regressions
        // did the scale change from

        size_t minCompScale = validRegressions->front().coeffs.scale;
        size_t compScale = minCompScale;
        size_t scale = compScale + 1;
        size_t startIdx = 0;

        RegressionGauss *upperReg;
        bool getNextReg = true;

        while (true)
        {
            if (getNextReg)
            {
                int idx = findNextReg(validRegressions, scale, startIdx);
                if (idx == -1)
                    break; // no scale to compare with exists

                startIdx = idx + 1; // necessary since otherwise the regression will be reset to itself
                upperReg = validRegressions->data() + idx;
                assert(upperReg->isValid);

                // the new regression is not necessarily at the searched input scale
                scale = upperReg->coeffs.scale;
                getNextReg = false;

                // part of the reset is considering all possible regressions again. It is
                // important to always do a full reset in case more than one valid peak
                // exists within a selected block.
                compScale = minCompScale;
            }
            assert(scale > compScale);

            Range_i regRange = findComparisonRegs(validRegressions, &upperReg->regSpan, compScale);

            // check the special case of only one regression at the lower scale
            // first since this should be the most common one
            if (regRange.startIdx == 1 && regRange.endIdx == 0)
            {
                // invalid state, do nothing @todo adjust minimum checked scale
            }
            else if (regRange.startIdx == regRange.endIdx)
            {
                RegressionGauss *lowerReg = validRegressions->data() + regRange.endIdx;
                assert(lowerReg->isValid);
                assert(upperReg->isValid);
                Range_i commonRange = {
                    min(lowerReg->regSpan.startIdx, upperReg->regSpan.startIdx),
                    max(lowerReg->regSpan.endIdx, upperReg->regSpan.endIdx)};
                size_t df = sumOfCumulative(df_cum, &commonRange) - 4; // this also applies for features

                double mseUpper = calcMSE_exp(&upperReg->coeffs, intensities, &commonRange, df);
                double mseLower = calcMSE_exp(&lowerReg->coeffs, intensities, &commonRange, df);

                if (mseLower < mseUpper)
                {
                    getNextReg = true;
                    upperReg->isValid = false;
                    lowerReg->numCompetitors += upperReg->numCompetitors + 1;
                }
                else
                {
                    lowerReg->isValid = false;
                    upperReg->numCompetitors += lowerReg->numCompetitors + 1;
                }
                assert(lowerReg->numCompetitors < int(validRegressions->size()));
                assert(upperReg->numCompetitors < int(validRegressions->size()));
            }
            else // this branch was never taken for a test dataset!
            {
                // there are multiple regressions at the lower scale which need to be combined into one mse
                // the leftmost relevant point is the first point of the first regression, the rightmost the
                // last point of the last regression. This is because regressions are always sorted by x0 in a scale.
                RegressionGauss *lowerRegS = validRegressions->data() + regRange.startIdx;
                RegressionGauss *lowerRegE = validRegressions->data() + regRange.endIdx;
                Range_i commonRange = {
                    min(lowerRegS->regSpan.startIdx, upperReg->regSpan.startIdx),
                    max(lowerRegE->regSpan.endIdx, upperReg->regSpan.endIdx)};
                size_t df = sumOfCumulative(df_cum, &commonRange);

                double mseUpper = calcMSE_exp(&upperReg->coeffs, intensities, &commonRange, df - 4);

                // slightly convoluted method for merging MSEs.
                double mseLower = multiMSE(validRegressions, intensities, &regRange, &commonRange, df);

                if (mseLower < mseUpper)
                {
                    addCompetitors(&regRange, validRegressions, upperReg->numCompetitors + 1);
                    upperReg->isValid = false;
                    getNextReg = true;
                }
                else
                {
                    // the greater-span regression persists and is checked against other contenders
                    size_t competitors = invalidateRange(&regRange, validRegressions);
                    upperReg->numCompetitors += competitors;
                }
            }

            // the lower scale is advanced until the higher scale is matched
            compScale += 1;
            if (compScale == scale)
            {
                // We have reached the highest possible scale to compare against. This
                // means that the checked regression is fully validated (or invalidated).
                getNextReg = true;
            }

            // sanity check: the number of competitors in valid regressions is the number of invalid regressions
            size_t invalidCount = 0;
            size_t competitors = 0;
            for (size_t i = 0; i < validRegressions->size(); i++)
            {
                RegressionGauss *reg = validRegressions->data() + i;
                if (reg->isValid)
                    competitors += reg->numCompetitors;
                else
                    invalidCount += 1;
            }
            assert(invalidCount == competitors);
        }

        return 0;
    }

    // --------------- old functions ---------------- //
    void mergeRegsInScale(
        const float *intensities,
        const unsigned int *const df_cum,
        std::vector<RegressionGauss> *validRegsTmp,
        std::vector<RegressionGauss> *validRegressions)
    {
        std::vector<RegressionGauss> validRegsAtScale;
        size_t currentScale = GLOBAL_MINSCALE;
        validRegsTmp->push_back({0}); // doing this avoids a second check for the last scale group
        validRegsTmp->back().coeffs.scale = 0;
        RegressionGauss *currentReg = validRegsTmp->data();

        while (currentReg->coeffs.scale != 0)
        {
        repeat:
            if (currentReg->coeffs.scale == currentScale)
            {
                validRegsAtScale.push_back(*currentReg);
                currentReg += 1;
                goto repeat;
            }

            // nothing happens if the per-scale vector is empty
            if (validRegsAtScale.size() == 1)
            {
                // only one valid regression at scale
                validRegressions->push_back(validRegsAtScale.front());
            }
            else if (validRegsAtScale.size() > 1)
            {
                // resolve conflicting regressions
                findBestScales(validRegressions, &validRegsAtScale, intensities, df_cum);
            }
            // regression is not incremented because a toggle was triggered
            currentScale += 1;
            validRegsAtScale.clear();
        }
        validRegsTmp->pop_back();
    }

    void findBestScales(std::vector<RegressionGauss> *validRegressions,
                        std::vector<RegressionGauss> *validRegsTmp,
                        const float *intensities,
                        const unsigned int *const degreesOfFreedom_cum)
    {
        /*
            Grouping:
            This block of code implements the grouping. It groups the valid peaks based
            on the apex positions. Peaks are defined as similar, i.e., members of the
            same group, if they fullfill at least one of the following conditions:
            - The difference between two peak apexes is less than 4. (Nyquist Shannon
            Sampling Theorem, separation of two maxima)
            - At least one apex of a pair of peaks is within the window of the other peak.
            (Overlap of two maxima)
        */
        // @todo could this part be combined with merge over scales?
        // vector with the access pattern [2*i] for start and [2*i + 1] for end point of a regression group
        std::vector<Range_i> groups;
        groups.reserve(validRegsTmp->size());

        size_t prev_i = 0;

        for (size_t i = 0; i < validRegsTmp->size() - 1; i++)
        {
            // check if the difference between two peak apexes is less than 4 (Nyquist Shannon
            // Sampling Theorem, separation of two maxima), or if the apex of a peak is within
            // the window of the other peak (Overlap of two maxima)
            RegressionGauss *reg1 = &(validRegsTmp->at(i));
            RegressionGauss *reg2 = &(validRegsTmp->at(i + 1));

            if (std::abs(reg1->apex_position - reg2->apex_position) < 4)
                continue;

            if (reg1->apex_position > reg2->regSpan.startIdx)
                continue;

            if (reg1->regSpan.endIdx > reg2->apex_position)
                continue;

            // the two regressions differ, i.e. create a new group
            groups.push_back({prev_i, i});
            prev_i = i + 1;
        }
        groups.push_back({prev_i, validRegsTmp->size() - 1}); // last group ends with index of the last element

        /*
          Survival of the Fittest Filter:
          This block of code implements the survival of the fittest filter. It selects the peak with
          the lowest mean squared error (MSE) as the representative of the group. If the group contains
          only one peak, the peak is directly pushed to the valid regressions. If the group contains
          multiple peaks, the peak with the lowest MSE is selected as the representative of the group
          and pushed to the valid regressions.
        */
        for (size_t groupIdx = 0; groupIdx < groups.size(); groupIdx++)
        {
            if (groups[groupIdx].startIdx == groups[groupIdx].endIdx)
            { // already isolated peak => push to valid regressions
                size_t regIdx = groups[groupIdx].startIdx;
                auto onlyReg = validRegsTmp->at(regIdx);
                assert(onlyReg.isValid);
                validRegressions->push_back(onlyReg);
            }
            else
            { // survival of the fittest based on mse between original data and reconstructed (exp transform of regression)
                RegPair bestRegIdx = findBestRegression(intensities, validRegsTmp, degreesOfFreedom_cum,
                                                        groups[groupIdx]);

                RegressionGauss bestReg = validRegsTmp->at(bestRegIdx.idx);
                bestReg.mse = bestRegIdx.mse;
                assert(bestReg.isValid);
                validRegressions->push_back(bestReg);
            }
        }
    }

    RegPair findBestRegression(
        const float *intensities,
        const std::vector<RegressionGauss> *regressions,
        const unsigned int *const degreesOfFreedom_cum,
        const Range_i regSpan)
    {
        double best_mse = INFINITY;
        unsigned int bestRegIdx = 0;

        // identify left (smallest) and right (largest) limit of the grouped regression windows
        size_t left_limit = -1;
        size_t right_limit = 0;
        for (size_t i = regSpan.startIdx; i < regSpan.endIdx + 1; i++)
        {
            const RegressionGauss *reg = &(*regressions)[i];
            left_limit = min(left_limit, reg->regSpan.startIdx);
            right_limit = max(right_limit, reg->regSpan.endIdx);
        }

        Range_i newRange = {left_limit, right_limit};
        double df_sum = sumOfCumulative(degreesOfFreedom_cum, &newRange);
        df_sum -= 4; // four coefficients
        assert(df_sum > 0);

        for (size_t i = regSpan.startIdx; i < regSpan.endIdx + 1; i++)
        {
            // step 2: calculate the mean squared error (MSE) between the predicted and actual values
            const RegressionGauss *reg = &(*regressions)[i];
            const Range_i range = {left_limit, right_limit}; // @todo this should update with the regressions
            double mse = calcMSE_exp(&reg->coeffs,
                                     intensities,
                                     &range,
                                     df_sum);
            if (mse < best_mse)
            {
                best_mse = mse;
                bestRegIdx = i;
            }
        }
        assert(regressions->at(bestRegIdx).isValid);
        return {bestRegIdx, float(best_mse)};
    }

    void mergeRegressionsOverScales(
        std::vector<RegressionGauss> *validRegressions,
        const float *intensities)
    {
        if (validRegressions->size() < 2)
        {
            return;
        }

        /*
          Grouping Over Scales:
          This block of code implements the grouping over scales. It groups the valid
          peaks based on the apex positions. Peaks are defined as similar, i.e.,
          members of the same group, if they fullfill at least one of the following conditions:
          - The difference between two peak apexes is less than 4. (Nyquist Shannon Sampling Theorem, separation of two maxima)
          - At least one apex of a pair of peaks is within the window of the other peak. (Overlap of two maxima)
        */

        // iterate over the validRegressions vector
        RegressionGauss *firstReg = validRegressions->data();
        for (size_t i = 0; i < validRegressions->size(); i++)
        {
            RegressionGauss *activeReg = firstReg + i;
            assert(activeReg->isValid);
            double MSE_group = 0;
            int DF_group = 0;
            // only calculate required MSEs since this is one of the performance-critical steps
            std::vector<float> exponentialMSE(validRegressions->size(), 0);
            std::vector<unsigned int> validRegressionsInGroup; // vector of indices to validRegressions
            validRegressionsInGroup.reserve(64);
            size_t competitors = 0; // a competitor is a mutually exclusive alternative regression

            for (size_t j = 0; j < i; j++)
            {
                RegressionGauss *secondReg = firstReg + j;
                if (!secondReg->isValid) // check is needed because regressions are set to invalid in the outer loop
                    continue;

                if (activeReg->apex_position < secondReg->regSpan.startIdx)
                    continue;

                if (activeReg->apex_position > secondReg->regSpan.endIdx)
                    continue;

                if (secondReg->apex_position < activeReg->regSpan.startIdx)
                    continue;

                if (secondReg->apex_position > activeReg->regSpan.endIdx)
                    continue;

                if (exponentialMSE[j] == 0.0)
                {
                    exponentialMSE[j] = calcMSE_exp(
                        &secondReg->coeffs,
                        intensities,
                        &secondReg->regSpan,
                        secondReg->df);
                }
                DF_group += secondReg->df;                      // add the degree of freedom
                MSE_group += exponentialMSE[j] * secondReg->df; // add the sum of squared errors
                // add the iterator of the ref peak to a vector of iterators
                validRegressionsInGroup.push_back(j);
                competitors += secondReg->numCompetitors + 1; // a regression can have beaten a previous one

            } // after this loop, validRegressionsInGroup contains all regressions that are still valid and contend with the regression at position i

            if (validRegressionsInGroup.empty()) // no competing regressions exist
            {
                assert(DF_group < 1);
                continue;
            }

            MSE_group /= DF_group;

            if (exponentialMSE[i] == 0.0)
            { // calculate the mse of the current peak
                exponentialMSE[i] = calcMSE_exp(
                    &activeReg->coeffs,
                    intensities,
                    &activeReg->regSpan,
                    activeReg->df);
            }
            if (exponentialMSE[i] < MSE_group)
            {
                // Set isValid to false for the candidates from the group
                for (size_t it_ref_peak : validRegressionsInGroup)
                {
                    firstReg[it_ref_peak].isValid = false;
                }
                // only advance competitor count if regression is actually better
                activeReg->numCompetitors = competitors;
            }
            else
            { // Set isValid to false for the current peak
                activeReg->isValid = false;
            }
        }

        // remove invalid regressions
        size_t accessID = 0;
        for (size_t i = 0; i < validRegressions->size(); i++)
        {
            if (validRegressions->at(i).isValid)
            {
                validRegressions->at(accessID) = validRegressions->at(i);
                accessID += 1;
            }
        }
        validRegressions->resize(accessID);
    }

#pragma endregion "eliminate conflicting regs"

    void findCoefficients(
        const float *intensity_log,
        const size_t length,
        size_t maxScale,
        std::vector<RegCoeffs> *coeffs)
    {
        assert(maxScale >= GLOBAL_MINSCALE);
        assert(maxScale <= MAXSCALE);

        assert(length > 4);
        maxScale = min(maxScale, (length - 1) / 2);

        size_t totalRegs = 0;
        for (size_t scale = GLOBAL_MINSCALE; scale <= maxScale; scale++)
        {
            totalRegs += length - 2 * scale;
        }
        coeffs->resize(totalRegs);

        // the first and last MINSCALE elements of the data do not need to be checked for x0, since they are invalid by definition
        const size_t limit = length - GLOBAL_MINSCALE;
        for (size_t x0 = GLOBAL_MINSCALE; x0 < limit; x0++)
        {
            const float *cen = intensity_log + x0; // this is initially the third real point

            // calculate the convolution with the kernel of the lowest scale - 1 (= 1), i.e. xT * intensity_log[cen - 1 : cen + 1]
            // the product sums are the rows of the design matrix (xT) * intensity_log[i:i+4] (dot product)
            // they are set to scale = 1 so the first value written is at scale = 2
            // b0 is 1, 1, 1,
            double product_sum_b0 = cen[-1] + cen[0] + cen[1];
            // b1 is -1, 0, 1
            double product_sum_b1 = -cen[-1] + cen[1];
            // b2 is 1, 0, 0
            double product_sum_b2 = cen[-1];
            // b3 is 0, 0, 1
            double product_sum_b3 = cen[1];

            size_t maxScale_absolute = 0;
            { // the largest valid scale depends on x0
                size_t maxScale_left = x0;
                size_t maxScale_right = length - x0 - 1;
                size_t maxScale_limits = min(maxScale_left, maxScale_right);
                maxScale_absolute = min(maxScale_limits, maxScale);
            }

            // var for access in inner loop
            size_t offset_prev = 0;
            for (size_t scale = GLOBAL_MINSCALE; scale <= maxScale_absolute; scale++)
            {
                // random access is difficult to vectorise
                double leftVal = cen[-scale];
                double rightVal = cen[scale];

                // expand the kernel to the left and right of the intensity_log.
                // b0 is expanded by the two outer points * 1
                product_sum_b0 += leftVal + rightVal;
                // b1 is expanded by the points * scale, negative to the left
                // product_sum_b1 += -double(scale) * leftVal + double(scale) * rightVal;
                product_sum_b1 += double(scale) * (rightVal - leftVal);
                // b2 and b3 are expanded by scale^2 the outermost point to the left or right
                double scale_sqr = double(scale * scale);
                product_sum_b2 += scale_sqr * leftVal;
                product_sum_b3 += scale_sqr * rightVal;

                const MatInverse inv = qalgo_matInverse[scale];

                const double inv_B_b0 = inv.B * product_sum_b0;
                const double inv_D_b1 = inv.D * product_sum_b1;

                // access is determined by scale and x0.
                // scale 2: idx is x0 - scale
                // scale 3: length - ((scale - 1) * 2) + x0 - scale
                // scale 4: 2 * length - (scale - 1) * 2 - (scale - 2) * 2 + x0 - scale
                //          2 * length - scale * 4 + 6 + x0 - scale
                const size_t offset_front = x0 - scale;
                const size_t access = offset_front + offset_prev;
                assert(access < totalRegs);
#define current (*coeffs)[access]
                current.b0 = inv.A * product_sum_b0 + inv.B * (product_sum_b2 + product_sum_b3);
                current.b1 = inv.C * product_sum_b1 + inv.D * (product_sum_b2 - product_sum_b3);
                current.b2 = inv_B_b0 + inv_D_b1 + inv.E * product_sum_b2 + inv.F * product_sum_b3;
                current.b3 = inv_B_b0 - inv_D_b1 + inv.F * product_sum_b2 + inv.E * product_sum_b3;
                current.scale = scale;
                current.x0 = x0;
#undef current
                // next scale would access front of vector
                offset_prev += length - 2 * scale;
            }
        }

        return;
    }

#pragma endregion "running regression"

#pragma region "validate Regression"

    double correctB0(const float *const intensities,
                     const Range_i *range,
                     std::vector<float> *predicted,
                     RegCoeffs *coeff)
    {
        // problem: after the log transform, regression residuals are not directly transferable to
        // the retransformed model. This is corrected by adjusting b0 so that the MSE in the
        // exponential case is minimal. We perform another regression taking e^(x b1 + x^2 b23) as a constant
        // and adjusting e^b0 so that (intensities - e^b0 * constant)^2 is minimal.
        // This is the same as scaling the regression by a constant c.

        // The (XtX)^-1 term collapses into a scalar, 1 / sum(predict^2). Xt * y is also scalar, sum(predict * intensities).
        // The corrective factor is then sum(predict * intensities) / sum(predict^2).

        // predict intensity only within range to prevent unnecessary exp operations.
        // prediction is incomplete, so it has to be multiplied with exp(b0)

        // uncertainties of this prediction are derived from the general equation for every coefficient.
        // the derivative is the following general form where bx' depends on the derived coefficient (1 for b0, x for b1, x^2 for b23)
        // b0' + ln(sum(exp(b0 + b1 x + b23 x^2) * y_i)) * sum(exp(b0 + b1 x + b23 x^2) * y_i * bx')
        //     - ln(sum(exp(b0 + b1 x + b23 x^2)^2)) * sum(exp(b0 + b1 x + b23 x^2)^2 * bx'^2)
        // note that the value of bx' depends on x and cannot be removed from the sum. Precalculation is also not
        // possible due to the value of y being required. This is not a large performance concern since we can combine it
        // with the calculations that are already necessary for obntaining the modified b0.

        double b0_old = coeff->b0;

        int start = int(range->startIdx);
        assert(0 <= start);
        int end = int(range->endIdx) + 1;
        assert(start < end);
        double x0 = double(coeff->x0);
        for (int i = start; i < end; i++)
        {
            double x = double(i) - x0;
            predicted->at(i) = regExp_fac(coeff, x);
        }

        double sum_predictSq = 0;
        double sum_predictReal = 0;
        double b0_exp = exp(coeff->b0);

        // Regression correction is only calculated from the range in which the regression is relevant initially.
        for (int i = start; i < end; i++)
        {
            double pred = predicted->at(i) * b0_exp;
            assert(pred > 0);
            sum_predictSq += pred * pred;
            sum_predictReal += pred * intensities[i];
        }
        double correction = sum_predictReal / sum_predictSq;

        // exp(a) * exp(b) == exp(a + b), so b0 + log(correction) is the same as predict * correction
        coeff->b0 += log(correction);
        // assert(abs(coeff->b0 - b0_old) < 0.001);

        // adjust the now incorrect values for predict. Remember that the previous prediciton was incomplete!
        double factor = exp(coeff->b0);
        for (int i = start; i < end; i++)
        {
            predicted->at(i) *= factor;
        }

        return factor;
    }

    invalid validRegWidth(const RegCoeffs *coeffs, Range_i *range)
    {
        // test regression validity without depending on b0 or the degrees of freedom
        const bool valley_left = coeffs->b2 >= 0;
        const bool valley_right = coeffs->b3 >= 0;
        if (valley_left && valley_right)
        {
            // there is no peak since both halves have a minimum
            return no_apex;
        }

        const bool apexLeft = coeffs->b1 < 0;
        // the apex cannot be on the same side as the valley point
        if (valley_left && apexLeft)
            return no_apex;
        if (valley_right && (!apexLeft))
            return no_apex;

        // position maximum / minimum of b2 or b3. This is just the frst derivative of the peak half equation
        // (d of y = b0 + b1 x + b23 x^2 => y b1 + 2 * b23 x) solved for y = 0
        double position_b2 = -coeffs->b1 / (2 * coeffs->b2);
        double position_b3 = -coeffs->b1 / (2 * coeffs->b3);
        double x0d = double(coeffs->x0);

        double apexd = apexLeft ? position_b2 : position_b3;

        // there must be at least two full points to either side of the apex @todo consider degrees of freedom here?
        if (abs(apexd) >= double(coeffs->scale - 1))
            return invalid_apex; // apex outside of regression window

        size_t apex = size_t(x0d + apexd);
        size_t lim_l = coeffs->x0 - coeffs->scale;
        size_t lim_r = coeffs->x0 + coeffs->scale;
        *range = {lim_l, lim_r};

        if (!(valley_left || valley_right))
            return ok; // all ok by definition

        // it is already established that the apex side of the regression is acceptable.
        // now, we only need to check that the valley is at least two points away from the apex
        // since only the difference is relevant, this is independent of the side on which the apex is
        assert(position_b2 < position_b3);
        if (position_b3 - position_b2 <= 2)
            return invalid_apex;

        if (valley_left)
        {
            lim_l = max(lim_l, coeffs->x0 - size_t(-position_b2));
        }
        else // if (valley_right)
        {
            lim_r = min(lim_r, coeffs->x0 + size_t(position_b3));
        }
        assert(lim_l <= apex);
        assert(lim_r >= apex);
        assert(lim_l < coeffs->x0);
        assert(lim_r > coeffs->x0);

        *range = {lim_l, lim_r};

        return ok;
    }

    void updateRegRange(const RegCoeffs *coeffs, const double valleyPos, Range_i *regSpan);

    double apexToEdgeRatio(const RegressionGauss *mutateReg, const float *intensities);

    invalid makeValidRegression( // returns the number of the failed test
        const float *intensities,
        const std::vector<float> *intensities_log,
        const std::vector<float> *predict,
        const size_t df_sum,
        const size_t length,
        RegressionGauss *mutateReg)

    {
        assert(!mutateReg->isValid);
        const size_t scale = mutateReg->coeffs.scale;
        const RegCoeffs coeffs = mutateReg->coeffs;

        assert(scale > 1);
        assert(coeffs.x0 + scale < length);

        // for a regression to be valid, at least one coefficient must be < 0
        if (coeffs.b2 >= 0 && coeffs.b3 >= 0)
        {
            return no_apex; // b0 independent
        }

        /*
          Apex and Valley Position Filter:
          This block of code implements the apex and valley position filter.
          It calculates the apex and valley positions based on the coefficients
          matrix B. If the test fails, one of the peak halves has less than two points.
        */
        double valley_position = 0;
        int failure = calcApexAndValleyPos(mutateReg, &valley_position);
        if (failure != 0) // something went wrong
        {
            return invalid_apex; // invalid apex and valley positions, b0 independent
        }

        Range_i comp = mutateReg->regSpan;
        updateRegRange(&coeffs, valley_position, &mutateReg->regSpan);
        assert(mutateReg->regSpan.endIdx < length);
        if (coeffs.x0 - mutateReg->regSpan.startIdx < GLOBAL_MINSCALE ||
            (mutateReg->regSpan.endIdx - coeffs.x0 < GLOBAL_MINSCALE))
        {
            // only one half of the regression applies to the data, since the
            // degrees of freedom for the "squished" half results in an invalid regression
            return no_df; // b0 independent
        }
        assert(comp.startIdx == mutateReg->regSpan.startIdx);
        assert(comp.endIdx == mutateReg->regSpan.endIdx);

        // this is the error term for the corrected regression. Of the original 4 x 4 matrix,
        // only the first row is needed
        double errorMat[4] = {0};

        // @todo new error correction here. Previous covariance matrix was mse (log) * (XtX)^-1,
        // multiply with matrix U, where first four terms are partial derivative of equation in
        // correctB0 by original coefficients

        /*
          Apex to Edge Filter:
          This block of code implements the apex to edge filter. It calculates
          the ratio of the apex signal to the edge signal and ensures that the
          ratio is greater than 2. This is a pre-filter for later
          signal-to-noise ratio checkups.
          @todo this is not a relevant test
        */
        float apexToEdge = apexToEdgeRatio(mutateReg, intensities);
        if (apexToEdge < 2)
        {
            // printf("apexToEdge triggered!\n");
            return invalid_apexToEdge; // invalid apex to edge ratio // b0 independent
        }

        // everything involving the RSS is dependent on b0!
        double RSS_log = calcRSS_log(mutateReg, intensities_log); // @todo we should use the exponential for this
        assert(RSS_log > 0);
        double RSS_exp = calcRSS(predict->data(), intensities, &mutateReg->regSpan);
        assert(RSS_exp > 0);
        double mse_exp = RSS_exp / double(df_sum);
        double mse_log = RSS_log / double(df_sum);

        mutateReg->mse = mse_log; // @todo the mse should be calculated in the function that uses it

        /*
        competing regressions filter:
        If the real distribution of points could also be described as a continuum (i.e. only b0 is relevant),
        the regression does not describe a peak. This is done through a nested F-test against a constant that
        is the mean of all predicted values. @todo test this function
        */
        // bool f_ok_log = f_testRegression(intensities_log, RSS_log, &mutateReg->regSpan);
        bool f_ok = f_testRegression(intensities, RSS_exp, &mutateReg->regSpan);
        // assert(f_ok == f_ok_exp);
        if (!f_ok)
        {
            return f_test_fail; // H0 holds, the two distributions are not noticeably different
        }

        if (!isValidQuadraticTerm(&mutateReg->coeffs, mse_log, df_sum)) // @todo replace with the exponential mse?
        {
            // this should be caught by the f test (?) @todo it is not caught by the F test
            // printf("invalid quadratic triggered\n");
            return invalid_quadratic; // statistical insignificance of the quadratic term
        }
        if (!isValidPeakArea(&mutateReg->coeffs, mse_log, df_sum)) // mse used in the "calcUncertainty" function call
        {
            return invalid_area; // statistical insignificance of the area
        }

        /*
          Height Filter:
          This block of code implements the height filter. It calculates the height
          of the peak based on the coefficients matrix B. Then it calculates the
          uncertainty of the height based on the Jacobian matrix and the variance-covariance
          matrix of the coefficients. @todo this is not what the function actually does!
        */
        calcPeakHeightUncert(mutateReg);                      // mse use, again
        if (1 / mutateReg->height_uncert <= T_VALUES[df_sum]) // statistical significance of the peak height
        {
            return invalid_height;
        }

        /*
          Area Filter:
          This block of code implements the area filter. It calculates the Jacobian
          matrix for the peak area based on the coefficients matrix B. Then it calculates
          the uncertainty of the peak area based on the Jacobian matrix.
          NOTE: this function does not consider b0: i.e. to get the real uncertainty and
          area multiply both with Exp(b0) later. This is done to avoid exp function at this point
        */
        // it might be preferential to combine both functions again or store the common matrix somewhere
        calcPeakAreaUncert(mutateReg);
        mutateReg->area = peakArea(coeffs.b0, coeffs.b1, coeffs.b2, coeffs.b3, 1); // @todo scale area later
        if (mutateReg->area / mutateReg->area_uncert <= T_VALUES[df_sum])
        {
            return invalid_area; // statistical insignificance of the area
        }

        /*
          Chi-Square Filter:
          This block of code implements the chi-square filter. It calculates the chi-square
          value based on the weighted chi squared sum of expected and measured y values in
          the exponential domain. If the chi-square value is less than the corresponding
          value in the CHI_SQUARES, the regression is invalid. @todo why?
        */
        double chiSquare = calcSSE_chisqared(mutateReg, intensities, predict);
        if (chiSquare > CHI_SQUARES[df_sum])
        {
            return invalid_chisq; // statistical insignificance of the chi-square value
        }

        calcUncertaintyPos(mutateReg);
        mutateReg->df = df_sum;
        mutateReg->apex_position += mutateReg->coeffs.x0;
        assert(mutateReg->apex_position > 1); // @todo this should be superfluous
        assert(mutateReg->apex_position < length - 1);
        assert(predict->size() == length);
        mutateReg->jaccard = calcJaccardIdx(intensities, predict->data(), length);

        mutateReg->isValid = true;
        return ok;
    }

    void updateRegRange(const RegCoeffs *coeffs, const double valleyPos, Range_i *regSpan)
    {
        const size_t scale = coeffs->scale;
        const size_t idxCenter = coeffs->x0;
        size_t rangeStart = idxCenter - scale;
        size_t rangeEnd = idxCenter + scale;

        *regSpan = {rangeStart, rangeEnd};
        if (valleyPos == 0) // no valley point exists
            return;

        // set start or end to the valley point if it is within the regression span
        size_t absValley = size_t(abs(valleyPos));
        bool valleyLeft = valleyPos < 0;
        bool updateVal = absValley < scale;
        if (!updateVal)
            return;

        if (valleyLeft)
        {
            regSpan->startIdx = idxCenter - absValley;
        }
        else
        {
            regSpan->endIdx = idxCenter + absValley;
        }
        assert(regSpan->endIdx > regSpan->startIdx);
    }

    double signalToNoise(const std::vector<float> *predict,
                         const Range_i *range,
                         const double mse)
    {
        /*
          basic idea of this test:
          under the assumption that the MSE is the total variance caused by noise, we can set a critical
          signal to noise ratio for every regression. This is done via the f-test. We compare the total
          variance of the signal to the total variance of the noise (mse). Since this function only returns
          the f-value, the test itself has to be performed in addition.
        */
        double varSignal = 0; // variance of predicted intensities

        size_t start = range->startIdx;
        size_t end = range->endIdx + 1;
        size_t len = end - start;
        double meanSignal = 0;
        for (size_t i = start; i < end; i++)
        {
            meanSignal += predict->at(i);
        }
        meanSignal /= len;

        for (size_t i = start; i < end; i++)
        {
            double diff = predict->at(i) - meanSignal;
            varSignal += diff * diff;
        }
        varSignal /= len;

        return varSignal / mse;
    }
#pragma endregion "validate Regression"

#pragma endregion "validate regression"

#pragma region "create peaks"

    MeanVar weightedMeanAndVariance_EIC(const std::vector<float> *weight,
                                        const std::vector<float> *values,
                                        const Range_i regSpan)
    {
        // weighted mean using intensity as weighting factor and left_limit right_limit as range
        size_t realPoints = 0;
        double mean_weights = 0;
        double sum_weighted_x = 0; // sum of values * weight
        double sum_weight = 0;
        for (size_t j = regSpan.startIdx; j <= regSpan.endIdx; j++)
        {
            int interpolated = (*values)[j] == 0 ? 0 : 1; // this is used instead of a continue so this can be vectorised. Skips a loop if value was interpolated
            mean_weights += (*weight)[j] * interpolated;
            sum_weighted_x += (*values)[j] * (*weight)[j] * interpolated;
            sum_weight += (*weight)[j] * interpolated;
            realPoints += 1 * interpolated; // interpolated points do not count!
        }
        mean_weights /= realPoints;
        sum_weighted_x /= mean_weights;
        sum_weight /= mean_weights;

        double weighted_mean = sum_weighted_x / sum_weight;
        double sum_Qxxw = 0.0; // sum of (values - mean)^2 * weight
        for (size_t j = regSpan.startIdx; j <= regSpan.endIdx; j++)
        {
            double difference = (*values)[j] - weighted_mean;
            double interpolated = (*values)[j] == 0 ? 0 : 1; // see above, add 0 if value is not real
            sum_Qxxw += interpolated * difference * difference * (*weight)[j];
        }
        float uncertaintiy = std::sqrt(sum_Qxxw / sum_weight / realPoints);
        return {float(weighted_mean), float(uncertaintiy)};
    };

    void createCentroidPeaks(
        std::vector<CentroidPeak> *peaks,
        const std::vector<RegressionGauss> *validRegressionsVec,
        const ProfileBlock *block, // @todo replace with mz and trace range
        const size_t scanNumber,
        const size_t ID_spectrum)
    {
        assert(!validRegressionsVec->empty());
        // iterate over the validRegressions vector
        for (size_t i = 0; i < validRegressionsVec->size(); i++)
        {
            const RegressionGauss *regression = validRegressionsVec->data() + i;
            assert(regression->isValid);
            CentroidPeak peak;
            peak.number_MS1 = scanNumber;
            // add height
            RegCoeffs coeff = regression->coeffs;
            peak.height = exp_approx_d(coeff.b0 + (regression->apex_position - coeff.x0) * coeff.b1 * 0.5);
            // peak height (exp(b0 - b1^2/4/b2)) with position being -b1/2/b2
            peak.heightUncertainty = regression->height_uncert * peak.height;

            size_t offset = (size_t)std::floor(regression->apex_position);
            double mz0 = block->mz[offset];
            double delta_mz = block->mz[offset + 1] - mz0;

            // add scaled area
            double scalar = exp(coeff.b0) * delta_mz;
            peak.area = regression->area * scalar;
            peak.areaUncertainty = regression->area_uncert * scalar;

            // add scaled apex position (mz)
            peak.mz = mz0 + delta_mz * (regression->apex_position - std::floor(regression->apex_position));
            peak.mzUncertainty = regression->position_uncert * delta_mz * T_VALUES[regression->df] * sqrt(1 + 1 / (regression->df + 4));

            // quality params
            peak.DQSC = 1 - erf_approx_f(regression->area_uncert / regression->area);
            // peak.df = regression->df;

            peak.numCompetitors = regression->numCompetitors;
            peak.scale = regression->coeffs.scale;

            // traceability information
            peak.trace.ID_spectrum = ID_spectrum;
            peak.trace.start = block->startPos;
            peak.trace.end = block->startPos + block->length - 1;

            peak.ID = peaks->size();

            peaks->push_back(peak);
        }
    }

    void createFeaturePeaks(
        std::vector<FeaturePeak> *peaks,
        const std::vector<RegressionGauss> *validRegressionsVec,
        const RT_Converter *convertRT,
        const std::vector<float> *RTs) // @todo this should be handled correctly through the converter
    {
        assert(!validRegressionsVec->empty());
        assert(peaks->empty());

        for (size_t i = 0; i < validRegressionsVec->size(); i++)
        {
            const RegressionGauss *regression = validRegressionsVec->data() + i;
            assert(regression->isValid);
            size_t maxIdx = convertRT->groups.back().interpolatedIndex;
            assert(regression->regSpan.endIdx <= maxIdx);

            FeaturePeak peak;
            // add height
            RegCoeffs coeff = regression->coeffs;
            // peak height (exp(b0 - b1^2/4/b2)) with position being -b1/2/b2
            peak.height = exp_approx_d(coeff.b0 + (regression->apex_position - coeff.x0) * coeff.b1 * 0.5);
            peak.heightUncertainty = regression->height_uncert * peak.height;

            // calculate the apex position in RT
            size_t idx_leftOfApex = (size_t)regression->apex_position;
            size_t idx_leftOfApex_absolute = convertRT->groups[idx_leftOfApex].interpolatedIndex;

            assert(idx_leftOfApex != 0); // at least two points to the left if apex is > 1
            size_t idx_rightOfApex_absolute = idx_leftOfApex_absolute + 1;

            if (idx_rightOfApex_absolute > convertRT->groups.size() - 1)
            {
                continue;
            }

            assert(idx_rightOfApex_absolute < convertRT->groups.size() - 1); // at least two points to the right
            // float rt_leftOfApex = convertRT->groups[idx_leftOfApex_absolute].trueRT;
            float rt_leftOfApex_true = RTs->at(idx_leftOfApex);
            // float rt_rightOfApex = convertRT->groups[idx_rightOfApex_absolute].trueRT;
            float rt_rightOfApex_true = RTs->at(idx_leftOfApex + 1);
            assert(rt_leftOfApex_true < rt_rightOfApex_true);
            float delta_rt = rt_rightOfApex_true - rt_leftOfApex_true;
            float rt_fraction = (regression->apex_position - floor(regression->apex_position));
            assert(rt_fraction >= 0);
            assert(rt_fraction < 1);
            float rt_apex = rt_leftOfApex_true + delta_rt * rt_fraction;
            peak.retentionTime = rt_apex;
            peak.RT_Uncertainty = regression->position_uncert * delta_rt;

            // add area
            float exp_b0 = exp_approx_d(coeff.b0); // exp(b0)
            peak.area = regression->area * exp_b0 * delta_rt;
            peak.areaUncertainty = regression->area_uncert * exp_b0 * delta_rt;

            // mz, mzUncertainty, mean DQSC and meanDQSF are all calculated in after this function is called in measurement_data
            peak.DQSF = 1 - erf_approx_f(regression->area_uncert / regression->area);

            assert(regression->regSpan.endIdx - coeff.x0 > 1);
            peak.idxPeakStart = regression->regSpan.startIdx;
            peak.idxPeakEnd = regression->regSpan.endIdx;
            peak.idxCenter_offset = coeff.x0 - regression->regSpan.startIdx;
            assert(peak.idxPeakEnd > peak.idxPeakStart);
            assert(peak.idxPeakEnd > peak.idxCenter_offset);
            assert(peak.idxPeakEnd - peak.idxPeakStart >= 4); // at least five points

            peak.coefficients = coeff;
            peak.mse_base = regression->mse;

            peak.scale = regression->coeffs.scale;
            peak.interpolationCount = rangeLen(&regression->regSpan) - regression->df - 4; // -4 since the degrees of freedom are reduced by 1 per coefficient
            peak.competitorCount = regression->numCompetitors;

            peaks->push_back(peak);
        }
    }
#pragma endregion "create peaks"

#pragma region calcSSE

    double calcRSS_log(const RegressionGauss *mutateReg, const std::vector<float> *observed)
    {
        const double idxCenter = mutateReg->coeffs.x0;
        double x = double(mutateReg->regSpan.startIdx) - idxCenter;
        double RSS = 0.0;
        for (size_t i = mutateReg->regSpan.startIdx; i < mutateReg->regSpan.endIdx + 1; i++)
        {
            // double new_x = double(i) - idxCenter;
            double pred = regAt(&mutateReg->coeffs, x);
            double difference = observed->at(i) - pred;
            x += 1.0;
            RSS += difference * difference;
        }
        return RSS;
    }

    double calcRSS_H0_cf1(const float *observed, const Range_i *range)
    {
        // this function calculates the RSS for H0: y = b0 (a constant value)

        double mean = 0;
        for (size_t i = range->startIdx; i <= range->endIdx; i++)
        {
            mean += observed[i];
        }
        mean /= rangeLen(range);

        double RSS = 0;
        for (size_t i = range->startIdx; i <= range->endIdx; i++)
        {
            double difference = observed[i] - mean;
            RSS += difference * difference;
        }

        return RSS;
    }

    double calcRSS_H0_cf2(const float *observed, const Range_i *range)
    {
        // this function calculates the RSS for H0: y = b0 + x * b1 (no weights)

        double slope, intercept;
        size_t length = rangeLen(range);
        const float *start = observed + range->startIdx;
        linReg_intx(start, length, &slope, &intercept);

        double RSS = 0;
        size_t x = 0;
        for (size_t i = range->startIdx; i <= range->endIdx; i++)
        {
            double difference = observed[i] - (intercept - slope * x);
            RSS += difference * difference;
            x += 1;
        }
        assert(x == length);

        return RSS;
    }

    bool f_testRegression(const float *observed, double RSS_reg, const Range_i *range)
    {
        // during the tests, the RSS for the regression has already been calculated in calcRSS_log
        assert(RSS_reg > 0);
        const size_t length = rangeLen(range);
        double RSS_H0 = INFINITY;
        bool f_ok = false;

        RSS_H0 = calcRSS_H0_cf1(observed, range); // y = b
        f_ok = F_test_regs(RSS_reg, RSS_H0, 4, 1, length, 0.05);
        if (!f_ok)
            return false;

        RSS_H0 = calcRSS_H0_cf2(observed, range); // y = b0 + b1 * x
        f_ok = F_test_regs(RSS_reg, RSS_H0, 4, 1, length, 0.05);
        if (!f_ok)
            return false; // reject H0, significant difference from

        // no alternatives were accepted
        return true;
    }

    double calcSSE_chisqared(const RegressionGauss *mutateReg,
                             const float *observed,
                             const std::vector<float> *predict)
    {
        double result = 0.0;
        for (size_t i = mutateReg->regSpan.startIdx; i < mutateReg->regSpan.endIdx + 1; i++)
        {
            double pred = predict->at(i);
            double obs = observed[i];
            double newdiff = (obs - pred) * (obs - pred);
            result += newdiff / pred; // this part is different from the above function, do not try to merge them!
        }
        return result;
    }

#pragma endregion calcSSE

    int calcApexAndValleyPos(
        RegressionGauss *mutateReg,
        double *valley_position)
    {
        const bool valley_left = mutateReg->coeffs.b2 >= 0;
        const bool valley_right = mutateReg->coeffs.b3 >= 0;
        const bool apexLeft = mutateReg->coeffs.b1 < 0;

        // the apex cannot be on the same side as the valley point
        if (valley_left && apexLeft)
            return 1;
        if (valley_right && (!apexLeft))
            return 1;

        // position maximum / minimum of b2 or b3. This is just the frst derivative of the peak half equation
        // (d of y = b0 + b1 x + b23 x^2 => y b1 + 2 * b23 x) solved for y = 0
        double position_b2 = -mutateReg->coeffs.b1 / (2 * mutateReg->coeffs.b2);
        double position_b3 = -mutateReg->coeffs.b1 / (2 * mutateReg->coeffs.b3);

        // the apex must be at least minscale different from the valley, if it exists
        if (valley_left || valley_right)
        {
            if ((position_b3 - position_b2) <= GLOBAL_MINSCALE)
                return 2;

            *valley_position = apexLeft ? position_b3 : position_b2;
        }
        // at this stage, the apex position is relating to x0. This is adjusted at a later point
        mutateReg->apex_position = apexLeft ? position_b2 : position_b3;

        // if this difference from 0 (rounded down) is exceeded by the apex position,
        // there are not enough points to validate the peak half left
        // example: at a scale of 2, the apex position must always be smaller than 1
        double maxApexDist = double(mutateReg->coeffs.scale - GLOBAL_MINSCALE + 1);
        if (abs(mutateReg->apex_position) > maxApexDist)
            return 3;

        return 0;
    }

    double apexToEdgeRatio(const RegressionGauss *mutateReg, const float *intensities)
    {
        // is the apex at least twice as large as the outermost point?
        // assumption: outermost point is already near base level
        const size_t idxApex = size_t(mutateReg->apex_position) + mutateReg->coeffs.x0;

        double left = intensities[mutateReg->regSpan.startIdx];
        double right = intensities[mutateReg->regSpan.endIdx];
        double minIntensity = min(left, right);
        double apex = intensities[idxApex]; // @todo since this is not the actual apex height, it might be a bad idea to use it
        return apex / minIntensity;
    }

    bool isValidQuadraticTerm(const RegCoeffs *coeffs, const double mse, const size_t df_sum)
    {
        // @todo explain
        // checks if both quadratic terms are significant - should we only check the apex coeff?
        assert(mse > 0);
        const double inv_E = qalgo_matInverse[coeffs->scale].E;
        double divisor = sqrt(inv_E * mse);
        double abs2 = abs(coeffs->b2);
        double abs3 = abs(coeffs->b3);
        double tValue = max(abs2, abs3);
        return tValue / divisor > T_VALUES[df_sum] * divisor; // statistical significance of the quadratic term
    }

    void calcPeakHeightUncert(RegressionGauss *mutateReg)
    {
        double Jacobian_height[4]{1, 0, 0, 0};         // Jacobian matrix for the height
        Jacobian_height[1] = mutateReg->apex_position; // apex_position * height;
        if (mutateReg->apex_position < 0)
        {
            Jacobian_height[2] = mutateReg->apex_position * mutateReg->apex_position; // apex_position * Jacobian_height[1];
            // Jacobian_height[3] = 0;
        }
        else
        {
            // Jacobian_height[2] = 0;
            Jacobian_height[3] = mutateReg->apex_position * mutateReg->apex_position; // apex_position * Jacobian_height[1];
        }
        // at this point without height, i.e., to get the real uncertainty
        // multiply with height later. This is done to avoid exp function at this point
        // mutateReg->height_uncert = calcUncertainty(Jacobian_height, mutateReg->coeffs.scale, mutateReg->mse);
        double uncertainty = sqrt(matProductReg(Jacobian_height, mutateReg->coeffs.scale) * mutateReg->mse);
        mutateReg->height_uncert = uncertainty;
        return;
    }

    bool isValidPeakHeight(
        const RegressionGauss *mutateReg,
        const double valley_position,
        const size_t df_sum,
        const double apexToEdge)
    {
        // check if the peak height is significantly greater than edge signal - deprecated!
        double apex = mutateReg->apex_position;
        assert(apex != 0);
        const bool apexLeft = apex < 0;
        // determine left or right limit, based on apex position
        const double scale_d = double(mutateReg->coeffs.scale);
        const double altValley = scale_d * (apexLeft ? -1 : 1);
        const double valley = valley_position == 0 ? altValley : valley_position;
        // construct matrix
        const double j1 = apex - valley;
        const double j23 = apex * apex - valley * valley;
        const double j2 = apexLeft ? j23 : 0;
        const double j3 = apexLeft ? 0 : j23;
        const double jacobianHeight[4]{0, j1, j2, j3};

        // float uncertainty_apexToEdge = calcUncertainty(jacobianHeight, mutateReg->coeffs.scale, mutateReg->mse);
        double uncertainty_apexToEdge = sqrt(matProductReg(jacobianHeight, mutateReg->coeffs.scale) * mutateReg->mse);

        // @todo why -2? Just to account for position?
        bool peakHeightSignificant = (apexToEdge - 2) > T_VALUES[df_sum] * (apexToEdge * uncertainty_apexToEdge);
        return peakHeightSignificant;
    }

    void calcPeakAreaUncert(RegressionGauss *mutateReg) // @todo only take coeffs as input, do not add shadow dependency on multiplication with e^b0
    {
        double b1 = mutateReg->coeffs.b1;
        double b2 = mutateReg->coeffs.b2;
        double b3 = mutateReg->coeffs.b3;
        double _SQRTB2 = 1 / sqrt(abs(b2));
        double _SQRTB3 = 1 / sqrt(abs(b3));
        double B1_2_SQRTB2 = b1 / 2 * _SQRTB2;
        double B1_2_SQRTB3 = b1 / 2 * _SQRTB3;
        double B1_2_B2 = b1 / 2 / b2;
        double B1_2_B3 = b1 / 2 / b3;

        // error calculation uses imaginary part if coefficient is > 0
        double err_L = (b2 < 0)
                           ? experfc(B1_2_SQRTB2, -1.0) // 1 - erf(b1 / 2 / SQRTB2) // ordinary peak
                           : dawson5(B1_2_SQRTB2);      // erfi(b1 / 2 / SQRTB2);        // peak with valley point;

        double err_R = (b3 < 0)
                           ? experfc(B1_2_SQRTB3, 1.0) // 1 + erf(b1 / 2 / SQRTB3) // ordinary peak
                           : dawson5(-B1_2_SQRTB3);    // -erfi(b1 / 2 / SQRTB3);       // peak with valley point ;

        // calculate the Jacobian matrix terms
        double J_1_common_L = _SQRTB2; // SQRTPI_2 * EXP_B12 / SQRTB2;
        double J_1_common_R = _SQRTB3; // SQRTPI_2 * EXP_B13 / SQRTB3;
        double J_2_common_L = B1_2_B2 / b1;
        double J_2_common_R = B1_2_B3 / b1;
        double J_1_L = J_1_common_L * err_L;
        double J_1_R = J_1_common_R * err_R;
        double J_2_L = J_2_common_L - J_1_L * B1_2_B2;
        double J_2_R = -J_2_common_R - J_1_R * B1_2_B3;

        double J[4]; // Jacobian matrix
        J[0] = J_1_R + J_1_L;
        J[1] = J_2_R + J_2_L;
        J[2] = -B1_2_B2 * (J_2_L + J_1_L / b1);
        J[3] = -B1_2_B3 * (J_2_R + J_1_R / b1);

        // at this point the area is without exp(b0), i.e., to get the real area multiply with exp(b0) later. This is done to avoid exp function at this point
        mutateReg->area = J[0];
        // mutateReg->area_uncert = calcUncertainty(J, mutateReg->coeffs.scale, mutateReg->mse);
        double uncertainty = sqrt(matProductReg(J, mutateReg->coeffs.scale) * mutateReg->mse);
        mutateReg->area_uncert = uncertainty;

        return;
    }

    bool isValidPeakArea(const RegCoeffs *coeffs, const double mse, const size_t df_sum)
    // @todo this function re-checks regression limits, we should always use the regression range for that
    {
        // function checks for coverage of the area by the regression
        double doubleScale = double(coeffs->scale);
        double b1 = coeffs->b1;
        double b2 = coeffs->b2;
        double b3 = coeffs->b3;

        double _SQRTB2 = 1 / sqrt(abs(b2));
        double _SQRTB3 = 1 / sqrt(abs(b3));
        double B1_2_B2 = b1 / 2 / b2;
        double B1_2_B3 = b1 / 2 / b3;

        double err_L_covered = 0;
        double x_left = -doubleScale;
        {
            double B1_2_SQRTB2 = b1 / 2 * _SQRTB2;
            double EXP_B12 = exp(-b1 * B1_2_B2 / 2);

            bool valley = b2 > 0;
            if (valley)
            {
                // valley point on the left side of the apex
                // error = erfi(b1 / 2 / SQRTB2)
                double err_L = dawson5(B1_2_SQRTB2);
                // check if the valley point is inside the window for the regression and consider it if necessary
                if (-B1_2_B2 < -doubleScale)
                {
                    // valley point is outside the window, use scale as limit
                    err_L_covered = err_L - erfi_qalgo((b1 - 2 * b2 * doubleScale) / 2 * _SQRTB2) * EXP_B12;
                }
                else
                {
                    // valley point is inside the window, use valley point as limit
                    err_L_covered = err_L;
                    x_left = -B1_2_B2;
                }
            }
            else
            {
                // no valley point
                // error = 1 - erf(b1 / 2 / SQRTB2)
                double err_L = experfc(B1_2_SQRTB2, -1.0);
                // ordinary peak half, take always scale as integration limit; we use erf instead of erfi due to the sqrt of absoulte value
                // erf((b1 - 2 * b2 * scale) / 2 / SQRTB2) + err_L - 1
                err_L_covered = erf_approx_f((b1 - 2 * b2 * doubleScale) / 2 * _SQRTB2) * EXP_B12 * SQRTPI_2 + err_L - SQRTPI_2 * EXP_B12;
            }
        }

        double err_R_covered = 0;
        double x_right = doubleScale; // right limit due to the window
        {
            double B1_2_SQRTB3 = b1 / 2 * _SQRTB3;
            double EXP_B13 = exp(-b1 * B1_2_B3 / 2);

            bool valley = b3 > 0;
            if (valley)
            {
                // valley point is on the right side of the apex
                // error = - erfi(b1 / 2 / SQRTB3)
                double err_R = dawson5(-B1_2_SQRTB3);
                if (-B1_2_B3 > doubleScale)
                {
                    // valley point is outside the window, use scale as limit
                    err_R_covered = erfi_qalgo((b1 + 2 * b3 * doubleScale) / 2 * _SQRTB3) * EXP_B13 + err_R;
                }
                else
                {
                    // valley point is inside the window, use valley point as limit
                    err_R_covered = err_R;
                    x_right = -B1_2_B3;
                }
            }
            else
            {
                // no valley point
                // error = 1 + erf(b1 / 2 / SQRTB3)
                double err_R = experfc(B1_2_SQRTB3, 1.0);
                // ordinary peak half, take always scale as integration limit; we use erf instead of erfi due to the sqrt of absoulte value
                // err_R - 1 - erf((b1 + 2 * b3 * scale) / 2 / SQRTB3)
                err_R_covered = err_R - SQRTPI_2 * EXP_B13 - erf_approx_f((b1 + 2 * b3 * doubleScale) / 2 * _SQRTB3) * SQRTPI_2 * EXP_B13;
            }
        }
        assert(x_left < x_right);

        double J_covered[4]; // Jacobian matrix for the covered peak area
        {
            // b0 not used on purpose, error is height independent
            const double y_left = exp_approx_d(b1 * x_left + b2 * x_left * x_left);
            const double y_right = exp_approx_d(b1 * x_right + b3 * x_right * x_right);
            const double dX = x_right - x_left;

            // calculate the Jacobian matrix terms
            const double J_1_common_L = _SQRTB2; // SQRTPI_2 * EXP_B12 / SQRTB2;
            const double J_1_common_R = _SQRTB3; // SQRTPI_2 * EXP_B13 / SQRTB3;
            const double J_2_common_L = B1_2_B2 / b1;
            const double J_2_common_R = B1_2_B3 / b1;

            // calculate the trapzoid correction terms for the jacobian matrix
            double trpzd_b0 = (y_right + y_left) * dX / 2;
            double trpzd_b1 = (x_right * y_right + x_left * y_left) * dX / 2;
            double trpzd_b2 = (x_left * x_left * y_left) * dX / 2;
            double trpzd_b3 = (x_right * x_right * y_right) * dX / 2;

            const double J_1_L_covered = J_1_common_L * err_L_covered;
            const double J_1_R_covered = J_1_common_R * err_R_covered;
            const double J_2_L_covered = J_2_common_L - J_1_L_covered * B1_2_B2;
            const double J_2_R_covered = -J_2_common_R - J_1_R_covered * B1_2_B3;

            J_covered[0] = J_1_R_covered + J_1_L_covered - trpzd_b0;
            J_covered[1] = J_2_R_covered + J_2_L_covered - trpzd_b1;
            J_covered[2] = -B1_2_B2 * (J_2_L_covered + J_1_L_covered / b1) - trpzd_b2;
            J_covered[3] = -B1_2_B3 * (J_2_R_covered + J_1_R_covered / b1) - trpzd_b3;
        }
        if (J_covered[0] < 0)
            return false;

        // float area_uncertainty_covered = calcUncertainty(J_covered, mutateReg->coeffs.scale, mutateReg->mse);
        double area_uncertainty_covered = sqrt(matProductReg(J_covered, coeffs->scale) * mse);

        // J[0] / uncertainty > Tval
        bool J_is_significant = J_covered[0] > T_VALUES[df_sum] * area_uncertainty_covered;

        return J_is_significant;
    }
#pragma endregion isValidPeakArea

    void calcUncertaintyPos(RegressionGauss *mutateReg) // @todo make this mse independent first
    {
        const double b1 = mutateReg->coeffs.b1;
        const double b2 = mutateReg->coeffs.b2;
        const double b3 = mutateReg->coeffs.b3;
        const double apex = mutateReg->apex_position;
        double J[4] = {0}; // Jacobian matrix

        J[1] = apex / b1;
        if (mutateReg->apex_position < 0)
        {
            J[2] = -apex / b2;
            // J[3] = 0;
        }
        else
        {
            // J[2] = 0;
            J[3] = -apex / b3;
        }

        // mutateReg->position_uncert = calcUncertainty(J, mutateReg->coeffs.scale, mutateReg->mse);
        double uncertainty = matProductReg(J, mutateReg->coeffs.scale);
        mutateReg->position_uncert = sqrt(uncertainty * mutateReg->mse);
        return;
    }

#pragma region "convolve regression"

    double matProductReg(const double J[4], const size_t scale)
    {
        // Calculate the Matrix Product of J * Xinv * J^T for uncertainty calculation
        const MatInverse inv = qalgo_matInverse[scale];
        double vecMatrxTranspose = J[0] * J[0] * inv.A +
                                   J[1] * J[1] * inv.C +
                                   (J[2] * J[2] + J[3] * J[3]) * inv.E +
                                   2 * (J[2] * J[3] * inv.F +
                                        J[0] * (J[1] + J[3]) * inv.B +
                                        J[1] * (J[2] - J[3]) * inv.D);
        return vecMatrxTranspose;
    }

#pragma endregion "convolve regression"

#pragma region "Feature Detection"

    double medianVec(const std::vector<float> *vec)
    {
        std::vector<float> tmpDiffs = *vec;
        std::sort(tmpDiffs.begin(), tmpDiffs.end());

        size_t size = tmpDiffs.size();
        assert(size > 0);
        if (size % 2 == 0)
        {
            return (tmpDiffs[size / 2 - 1] + tmpDiffs[size / 2]) / 2;
        }
        else
        {
            return tmpDiffs[size / 2];
        }
    }

    // RT_Converter interpolateScanNumbers(const std::vector<float> *retentionTimes)
    // {
    //     // This function interpolates the existing RTs of MS1 spectra and produces a vector that contains,
    //     // for the index of every real MS1 spectrum, the scan number with interpolations for that spectrum.
    //     // Since we work with an integer scale, interpolations are rounded at 0.5.
    //     // It is assumed that the RT vector is sorted.

    //     assert(!retentionTimes->empty());
    //     std::vector<float> diffs(retentionTimes->size() - 1, 0);
    //     for (size_t i = 0; i < diffs.size(); i++)
    //     {
    //         diffs[i] = retentionTimes->at(i + 1) - retentionTimes->at(i);
    //     }
    //     assert(!diffs.empty());

    //     // @todo this should be the mode, not the median
    //     float expectedDiff = medianVec(&diffs);
    //     // if this is not given, there are severe distortions at some point in the data. We accept one non-compliance,
    //     // at the first scan in the spectrum only, which is generally one of the first recorded scans
    //     // incorrectness here could be due to a vendor-specific instrument checkup procedure, which has been observed once at least
    //     // assert(tmpDiffs[1] > expectedDiff / 2);

    //     size_t scanCount = retentionTimes->size();
    //     std::vector<RT_Grouping> totalRTs;
    //     totalRTs.reserve(scanCount * 2);
    //     std::vector<size_t> idxToGrouping(scanCount, -1);

    //     float critDiff = expectedDiff * 1.5; // if the difference is 1.5 times greater than the critDiff, there should be at least one interpolation

    //     for (size_t i = 0; i < diffs.size(); i++)
    //     {
    //         size_t interpScan = totalRTs.size();
    //         totalRTs.push_back({i, interpScan, (*retentionTimes)[i], false});
    //         idxToGrouping[i] = interpScan;

    //         if (diffs[i] > critDiff + FLT_EPSILON) // this is necessary since for truly equidistant data, critDiff can be greater than diff even if they should be identical
    //         {
    //             // interpolate at least one point
    //             size_t numInterpolations = size_t(diffs[i] / expectedDiff + 0.5 - FLT_EPSILON); // + 0.5 since value is truncated (round up), see above for epsilon
    //             assert(numInterpolations != 0);

    //             float RTstep = diffs[i] / (numInterpolations);
    //             for (size_t j = 1; j < numInterpolations; j++) // +1 since the span is between two points
    //             {
    //                 interpScan = totalRTs.size();
    //                 float newRT = (*retentionTimes)[i] + RTstep * j;
    //                 totalRTs.push_back({size_t(-1), interpScan, newRT, true});
    //             }
    //         }
    //     }
    //     // add in the last remaining point
    //     size_t lastScan = totalRTs.size();
    //     idxToGrouping.back() = lastScan;
    //     totalRTs.push_back({diffs.size(), lastScan, retentionTimes->back(), false});

    //     return {totalRTs, idxToGrouping};
    // }

    void fillPeakVals(EIC *eic, FeaturePeak *currentPeak)
    {
        currentPeak->scanPeakStart = eic->scanNumbers.front();
        currentPeak->scanPeakEnd = eic->scanNumbers.back();

        // the correct limits in the non-interpolated EIC need to be determined. They are already included
        // in the cumulative degrees of freedom, but since there, df 0 is outside the EIC, we need to
        // use the index df[limit] - 1 into the original, non-interpolated vector

        size_t limit_L = eic->df[currentPeak->idxPeakStart];
        limit_L = min(limit_L, limit_L - 1); // uint underflows, so no issues.
        size_t limit_R = eic->df[currentPeak->idxPeakEnd] - 1;
        assert(limit_L < limit_R);

        currentPeak->idxBinStart = limit_L;
        currentPeak->idxBinEnd = limit_R;

        Range_i regSpan = {limit_L, limit_R};

        auto tmp = weightedMeanAndVariance_EIC(&eic->ints_area, &eic->mz, regSpan);
        currentPeak->mz = tmp.mean;
        currentPeak->mzUncertainty = tmp.var;
        currentPeak->DQSC = weightedMeanAndVariance_EIC(&eic->ints_area, &eic->DQSC, regSpan)
                                .mean;
        currentPeak->DQSB = weightedMeanAndVariance_EIC(&eic->ints_area, &eic->DQSB, regSpan)
                                .mean;
    }

    std::vector<FeaturePeak> findFeatures(std::vector<EIC> &EICs,
                                          const RT_Converter *convertRT)
    {
        // reset the logger, not the upcount
        resetLogger(&globalLogStruct);
        globalLogStruct.features = true;
        globalLogStruct.checkedRegionCount = EICs.size();

        // @todo this is not a universal limit and only chosen for computational speed at the moment
        // with an estimated scan difference of 0.6 s this means the maximum peak width is 61 * 0.6 = 36.6 s
        static const size_t GLOBAL_MAXSCALE_FEATURES = 30;
        assert(GLOBAL_MAXSCALE_FEATURES <= MAXSCALE);

        std::vector<FeaturePeak> peaks;    // return vector for feature list
        peaks.reserve(EICs.size() / 4);    // should be enough to fit all features without reallocation
        std::vector<FeaturePeak> tmpPeaks; // add features to this before pasting into FL

        std::vector<RegressionGauss> validRegressions;
        validRegressions.reserve(512);

        std::vector<float> logIntensity;

        for (size_t i = 0; i < EICs.size(); ++i)
        {
            auto currentEIC = EICs[i];
            if (currentEIC.scanNumbers.size() < 5)
            {
                continue; // skip due to lack of data, i.e., degrees of freedom will be zero
            }

            validRegressions.clear();
            size_t length = currentEIC.df.size();
            assert(length > 4); // data must contain at least five points

            logIntensity.resize(length);
            logIntensity.clear();

            size_t maxScale = min(GLOBAL_MAXSCALE_FEATURES, (length - 1) / 2);

            runningRegression(
                currentEIC.ints_area.data(),
                &logIntensity,
                currentEIC.df.data(),
                currentEIC.ints_area.size(),
                maxScale,
                &validRegressions);

            if (!validRegressions.empty())
            {
                createFeaturePeaks(&tmpPeaks, &validRegressions, convertRT, &currentEIC.RT);

                for (auto peak : tmpPeaks)
                {
                    assert(peak.retentionTime > currentEIC.RT.front());
                    assert(peak.retentionTime < currentEIC.RT.back());
                }
            }
            // @todo extract the peak construction here and possibly extract findFeatures into a generic function

            if (tmpPeaks.empty())
            {
                continue;
            }
            for (size_t j = 0; j < tmpPeaks.size(); j++)
            {
                FeaturePeak currentPeak = tmpPeaks[j];

                // the correct limits in the non-interpolated EIC need to be determined. They are already included
                // in the cumulative degrees of freedom, but since there, df 0 is outside the EIC, we need to
                // use the index df[limit] - 1 into the original, non-interpolated vector

                currentPeak.idxBin = i;

                fillPeakVals(&currentEIC, &currentPeak);
                assert(currentPeak.scanPeakEnd < convertRT->groups.size());

                peaks.push_back(currentPeak);
            }

            tmpPeaks.clear();
        }

        // advance logger by one dataset
        globalLogStruct.totalPeakCount = peaks.size();
        loggerPrint(&globalLogStruct);
        globalLogStruct.dataID += 1;

        // peaks are sorted here so they can be treated as const throughout the rest of the program
        std::sort(peaks.begin(), peaks.end(), [](FeaturePeak lhs, FeaturePeak rhs)
                  { return lhs.retentionTime < rhs.retentionTime; });
        return peaks;
    }

#pragma endregion "Feature Detection"

#pragma region "find centroids"

    size_t getProfileRegions(
        std::vector<ProfileBlock> *groupedData,
        const std::vector<float> *spectrum_mz,
        const std::vector<float> *spectrum_int);

    // int findCentroids(XML_File &data, // @todo get rid of the direct coupling to pugixml
    //                   const std::vector<unsigned int> *selectedIndices,
    //                   std::vector<CentroidPeak> *centroids)
    // {
    //     // resetting the logger does not change the number of runs or the features field
    //     resetLogger(&globalLogStruct);
    //     globalLogStruct.features = false;

    //     const size_t countSelected = selectedIndices->size();

    //     std::vector<float> spectrum_mz(1000);
    //     std::vector<float> spectrum_int(1000);

    //     std::vector<ProfileBlock> subprofiles(500);
    //     for (size_t scanNum = 0; scanNum < countSelected; ++scanNum)
    //     {
    //         // avoid needless allocation / deallocation of otherwise scope-local vectors
    //         spectrum_mz.clear();
    //         spectrum_int.clear();
    //         subprofiles.clear();

    //         size_t ID_spectrum = selectedIndices->at(scanNum);
    //         data.get_spectrum(&spectrum_mz, &spectrum_int, ID_spectrum);

    //         size_t maxWindowSize = getProfileRegions(&subprofiles, &spectrum_mz, &spectrum_int);
    //         if (maxWindowSize == 0) // this is also relevant to filtering, add a warning if no filter?
    //             continue;

    //         findCentroidPeaks(centroids, &subprofiles, scanNum, ID_spectrum, maxWindowSize);

    //         // add sum statistics to log struct
    //         globalLogStruct.checkedRegionCount += subprofiles.size();
    //     }

    //     globalLogStruct.totalPeakCount += centroids->size();

    //     // log printing has to be handled by the logged function itself
    //     // data ID is only incremented after feature processing
    //     loggerPrint(&globalLogStruct);

    //     for (unsigned int i = 0; i < centroids->size(); i++)
    //     {
    //         centroids->at(i).ID = i;
    //     }

    //     return centroids->size();
    // }

    int getNextProfileRegion(
        const std::vector<float> *spectrum_mz,
        const std::vector<float> *spectrum_int,
        ProfileBlock *block);

    CentroidPeak peakToCen(const PeakFit *peak, size_t id, size_t specNum)
    {
        CentroidPeak cen = {0};

        cen.area = peak->area;
        cen.areaUncertainty = peak->area_uncert;
        cen.DQSC = peak->dqs;
        cen.height = peak->height;
        cen.heightUncertainty = peak->height_uncert;
        cen.ID = id;
        cen.mz = peak->position;
        cen.mzUncertainty = peak->position_uncert;
        cen.number_MS1 = specNum;
        cen.scale = peak->coeffs.scale;
        cen.width = peak->fwhm;

        return cen;
    }

    // int findCentroids_new(XML_File &data, // @todo move this into the xml file specific part of qAlgorithms
    //                       const std::vector<unsigned int> *selectedIndices,
    //                       std::vector<CentroidPeak> *centroids)
    // {
    //     const size_t countSelected = selectedIndices->size();

    //     std::vector<float> spectrum_mz;
    //     spectrum_mz.reserve(1000);
    //     std::vector<float> spectrum_int;
    //     spectrum_int.reserve(1000);

    //     size_t scanNum = 0;
    //     std::vector<PeakFit> ret; // @todo add tracking stuff for retaining spectrum / startpoint information
    //     ProfileBlock block;
    //     size_t ID_spectrum = selectedIndices->front();
    //     while (scanNum < countSelected)
    //     {
    //         if (getNextProfileRegion(&spectrum_mz, &spectrum_int, &block))
    //         {
    //             // avoid needless allocation / deallocation of otherwise scope-local vectors
    //             spectrum_mz.clear();
    //             spectrum_int.clear();

    //             ID_spectrum = selectedIndices->at(scanNum);
    //             data.get_spectrum(&spectrum_mz, &spectrum_int, ID_spectrum);

    //             scanNum += 1;
    //         }
    //         // at this point, the block contains one continuus region of points in the source spectrum

    //         qpeaks_find(block.intensity, block.mz, nullptr, block.length, GLOBAL_MAXSCALE_CENTROID, &ret);
    //     }
    //     // move peak data structure into centroid struct @todo remember spectrum id

    //     for (unsigned int i = 0; i < centroids->size(); i++)
    //     {
    //         centroids->at(i).ID = i;
    //     }

    //     return centroids->size(); // @todo rework centroids to be multiple separate arrays
    // }

    // size_t removedPoints = 0;

    // size_t getProfileRegions(
    //     std::vector<ProfileBlock> *groupedData,
    //     const std::vector<float> *spectrum_mz,
    //     const std::vector<float> *spectrum_int)
    // {
    //     assert(groupedData->empty());
    //     assert(spectrum_int->size() == spectrum_mz->size());
    //     // this function walks through both spectra and fills a vector of blocks
    //     // every block is a region which is searched for centroids. Further, the
    //     // edges of the block are interpolated such that the third real point (which
    //     // is the start of valid regressions) can be occupied with a symmetric regression
    //     // of the maximum scale (currently 8)
    //     // the function returns the largest created block.

    //     const float minIntensity = 10; // @todo bad solution?

    //     size_t maxLength = 0;
    //     size_t currLen = 0; // current length
    //     size_t lastZero = 0;
    //     for (size_t pos = 0; pos < spectrum_int->size(); pos++)
    //     {
    //         if (spectrum_int->at(pos) < minIntensity)
    //         {
    //             if (currLen >= 5)
    //             {
    //                 size_t idxL = lastZero + 1;
    //                 size_t idxR = pos - 1;

    //                 ProfileBlock b = {
    //                     spectrum_int->data() + idxL,
    //                     spectrum_mz->data() + idxL,
    //                     idxL,
    //                     idxR - idxL + 1};
    //                 groupedData->push_back(b);

    //                 maxLength = maxLength > currLen ? maxLength : currLen;
    //             }
    //             currLen = 0;
    //             lastZero = pos;
    //         }
    //         else
    //         {
    //             currLen += 1;
    //         }
    //     }
    //     if (currLen >= 5)
    //     {
    //         size_t idxL = lastZero + 1;
    //         size_t idxR = spectrum_int->size() - 1;
    //         ProfileBlock b = {
    //             spectrum_int->data() + idxL,
    //             spectrum_mz->data() + idxL,
    //             idxL,
    //             idxR - idxL + 1};
    //         groupedData->push_back(b);

    //         maxLength = maxLength > currLen ? maxLength : currLen;
    //     }

    //     assert(maxLength > 4);
    //     return maxLength;
    // }

    // int getNextProfileRegion(
    //     const std::vector<float> *spectrum_mz,
    //     const std::vector<float> *spectrum_int,
    //     ProfileBlock *block)
    // {
    //     size_t idx = block->startPos + block->length; // first element past the previous block
    //     size_t len = spectrum_int->size();

    //     if (idx >= len || len - idx < 5)
    //     {
    //         block->intensity = nullptr;
    //         block->mz = nullptr;
    //         block->length = 0;
    //         block->startPos = 0;
    //         return 1;
    //     }

    //     const float minIntensity = 10; // @todo bad solution?

    //     for (; idx < len; idx++)
    //     {
    //         if (spectrum_int->at(idx) > minIntensity)
    //             break;
    //     }
    //     const float *intensity = spectrum_int->data() + idx;
    //     const float *mz = spectrum_mz->data() + idx;

    //     size_t length = 0;
    //     for (size_t i = idx; i < len; i++)
    //     {
    //         length += 1;
    //         if (spectrum_int->at(i) < minIntensity)
    //             break;
    //     }

    //     block->intensity = intensity;
    //     block->mz = mz;
    //     block->startPos = idx;
    //     block->length = length;
    //     return 0;
    // }

    inline double regAt(const RegCoeffs *coeff, const double x)
    {
        const double b23 = x < 0 ? coeff->b2 : coeff->b3;
        return coeff->b0 + (coeff->b1 + x * b23) * x;
    }

    double regExp_fac(const RegCoeffs *coeff, const double x)
    {
        double b23 = x < 0 ? coeff->b2 : coeff->b3;
        return exp((coeff->b1 + x * b23) * x);
    }

    double fullWidthHalfMax(const RegCoeffs *coeff, const double height, const double delta_x)
    {
        // solve height / 2 = exp(b0 + x b1 + x^2 b2)
        double y = log(height / 2);
        // use quadratic formula for a x^2 + b x + c = 0
        double a_l = coeff->b2;
        double a_r = coeff->b3;
        double b = coeff->b1;
        double c = coeff->b0 - y;

        double x_l, x_r, dummy;
        solveQuadratic(a_l, b, c, &x_l, &dummy);
        solveQuadratic(a_r, b, c, &dummy, &x_r);

        assert(x_l < x_r);

        return (x_r - x_l) * delta_x;
    }

    double FWHM_to_sdev(const double fwhm)
    {
        // equivalent standard deviation for a symmetrical gaussian peak
        // for the gaussian peak, fwhm only depends on the standard deviation.
        // this function returns the standard deviation based on that model.

        // FWHM = 2 * sqrt(-log(0.5) * 2 * sdev^2)
        // (FWHM / 2)^2 = -log(0.5) * 2 * sdev^2
        // sdev^2 = (FWHM / 2)^2 / (-log(0.5) * 2)
        // sdev = sqrt((FWHM / 2)^2 / (-log(0.5) * 2))
        // sdev = (FWHM / 2) / sqrt(-log(0.5) * 2))
        // sdev = FWHM / (sqrt(-log(0.5) * 2)) * 2)

        // constexpr double divisor = sqrt(-log(0.5) * 2.0) * 2.0;
        const double divisor = 2.354820045030949;
        return fwhm / divisor;
    }

    // base function: integral e^(b0 + b1 x + b2 x^2) dx =
    // [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * erfi( (b1 + 2 b2 x) / (2 sqrt(b2)) ) ] / (2 sqrt(b2))   // source: wolfram alpha
    // erfi(x) = i * erf(i * x)
    // under the condition b2 < 0: sqrt(b2) = i * sqrt(-b2), where sqrt(-b2) is a real number
    // [(real part) * i * erf( i * (real) / (i * 2 * sqrt(-b2)) ) ] / (i * 2 * sqrt(-b2))
    // i within and outside of the error function cancel each other out: i / i = 1
    // erf(0) = 0 ; erf(-inf) = -1 ; erf(inf) = 1
    // assuming b2 < 0:
    // F(-inf) = [ (...) * erf( (b1 + 2 b2 * -inf) / (> 0) ) ] / (...)
    // b2 * -inf = +inf ; b2 * +inf = -inf
    // F(-inf) = [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * +1 ] / (2 sqrt(-b2))
    // F(+inf) = [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * -1 ] / (2 sqrt(-b2))
    // F(0) =    [ sqrt(pi) * e^( b0 - b1^2/(4 b2) ) * erf( b1 / (2 sqrt(-b2)) ) ] / (2 sqrt(-b2))
    double peakArea(const double b0, const double b1, const double b2, const double b3, const double delta_x)
    {
        // reasoning: the regression operates on the transformed x-axis x' = (x - x0) / delta_x.
        // In order for the area to be correct, the position of the regression is irrelevant, so
        // we can assume x0 = 0. This means that for the retransformation, we must use:
        // y = exp(b0 + b1 * (x / dx) + b2 * (x / dx)^2)
        // This can be factored out to modified coefficients (b0 is not modified).
        // double b1 = b1_in * delta_x;
        // double b2 = b2_in * delta_x * delta_x;
        // double b3 = b3_in * delta_x * delta_x;
        // @todo this is not working in practice. Observation shows that the transformation from
        // the predicted area at delta_x = 1 (A1) to the area in the untransformed x axis is A = A1 * delta_x
        // Since that is a sensible result, we will stick to the single multiplication and unmodified coeffs

        bool b2_pos = b2 > 0;
        bool b3_pos = b3 > 0;
        assert(!(b2_pos && b3_pos));

        // only relevant if one of the coefficients is > 0: Calculate the valley position for a given
        // peak. This is required during the calculation of the erfi term. Always choosing the valley
        // will lead to distortion, but it should overall be lower than if the peak is always cut off
        // at the end of the peak window.
        double x = -b1 / (2 * (b2_pos ? b2 : b3));

        // note: this is not factored out into a separate function since + and - inf use different terms
        double dsqrt_b2 = 2 * sqrt(-b2);
        double dsqrt_b3 = 2 * sqrt(-b3);
        double eterm_b2 = exp(b0 - (b1 * b1) / (4 * b2));
        double eterm_b3 = exp(b0 - (b1 * b1) / (4 * b3));
        const double sqrt_pi = 1.7724538509055158819; // sqrt(M_PI);

        double F_b2_ninf = (sqrt_pi * eterm_b2 * 1) / dsqrt_b2;
        // double F_b2_inf = (sqrt_pi * eterm_b2 * -1) / dsqrt_b2;
        double error_b2 = b2_pos ? liberfc::erfi((b1 + 2 * b2 * x) / dsqrt_b2) : erf(b1 / dsqrt_b2);
        double F_b2_zero = (sqrt_pi * eterm_b2 * error_b2) / dsqrt_b2;
        double area_L = F_b2_zero - F_b2_ninf;

        // double F_b3_ninf = (sqrt_pi * eterm_b3 * 1) / dsqrt_b3;
        double F_b3_inf = (sqrt_pi * eterm_b3 * -1) / dsqrt_b3;
        double error_b3 = b3_pos ? liberfc::erfi((b1 + 2 * b3 * x) / dsqrt_b3) : erf(b1 / dsqrt_b3);
        double F_b3_zero = (sqrt_pi * eterm_b3 * error_b3) / dsqrt_b3;
        double area_R = F_b3_inf - F_b3_zero;

        assert((area_L > 0) == (area_R > 0));
        // there are cases where both areas are negative, but result in the correct value when summed anyway
        // assert(area_L > 0);
        // assert(area_R > 0);
        area_L = abs(area_L);
        area_R = abs(area_R);
        double area_F = (area_L + area_R) * delta_x;
        return area_F;
    }

    double peakAreaUncert(const RegCoeffs *coeffs, const double delta_x, const double area)
    {
        // uses the jacobi matrix for uncertainty calculation of the area
        // refer to the supplement to qPeaks
        double jacobi[4] = {area / delta_x, 0, 0, 0};
        jacobi[1] = (exp(coeffs->b0) / 2 - jacobi[0] * coeffs->b1 / 2) * (1 / coeffs->b2 + 1 / coeffs->b3);
        jacobi[2] = -1 / (2 * coeffs->b2) * (jacobi[0] + jacobi[1] * coeffs->b1);
        jacobi[3] = -1 / (2 * coeffs->b3) * (jacobi[0] + jacobi[1] * coeffs->b1);

        // @todo
    }
}