/*
 * SPDX-FileCopyrightText: 2016-2025 Rick Helmus <r.helmus@uva.nl>
 *
 * SPDX-License-Identifier: GPL-3.0-only
 */

#ifndef PATROON_SC_H
#define PATROON_SC_H

#include <Rcpp.h>

#include "msdata.h"

#include "utils-xml.h" // for pugixml
#include "StreamCraft/StreamCraft_lib.h"

class MSReadBackendSC: public MSReadBackend
{
    std::unique_ptr<sc::MS_FILE> handle;
    
    void doOpen(const std::string &) override;
    void doClose(void) override;
    ThreadDataType doGetThreadData(void) const override { return nullptr; }
    SpectrumRaw doReadSpectrum(const ThreadDataType &tdata, SpectrumRawTypes::MSLevel MSLevel,
                               const SpectrumRawSelection &scanSel,
                               const SpectrumRawTypes::MobilityRange &mobRange,
                               SpectrumRawTypes::Intensity minIntensityIMS) const override;
    
public:
    MSReadBackendSC(void) { }
    ~MSReadBackendSC(void) { }
    
    void generateSpecMetadata(void);
};

RCPP_EXPOSED_CLASS(MSReadBackendSC)

#endif
