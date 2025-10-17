/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#ifndef PATROON_MSTK_H
#define PATROON_MSTK_H

#ifdef WITH_MSTK

#include <Rcpp.h>

#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

#include "msdata.h"

class MSReadBackendMSTK: public MSReadBackend
{
    static int backends; // UNDONE: for debugging

    void doOpen(const std::string &) override { }
    void doClose(void) override { }
    ThreadDataType doGetThreadData(void) const override;
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                               const SpectrumRawSelection &scanSel,
                               const SpectrumRawTypes::MobilityRange &mobRange,
                               SpectrumRawTypes::Intensity minIntensityIMS) const override;
    
public:
    MSReadBackendMSTK(void) { ++backends; Rcpp::Rcout << "constr: backends:" << backends << "\n"; }
    ~MSReadBackendMSTK(void) { --backends; Rcpp::Rcout << "destr: backends:" << backends << "\n"; }

    void generateSpecMetadata(void);
    int getBackends(void) const { return backends; }
};

//void writeMS1SpectraMSTK(std::string path, const std::vector<SpectrumRaw> &spectra, const SpectrumRawMetadataMS &meta);

#else

#include "msdata.h"

class MSReadBackendMSTK: public MSReadBackendUnavailable
{
public:
    MSReadBackendMSTK(void) : MSReadBackendUnavailable("mstoolkit") { }
    ~MSReadBackendMSTK(void) { }
    void generateSpecMetadata(void) { };
    int getBackends(void) const { return 0; }
};

#endif // WITH_MSTK

RCPP_EXPOSED_CLASS(MSReadBackendMSTK)


#endif
