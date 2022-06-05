/*
 *   OpenTIMS: a fully open-source library for opening Bruker's TimsTOF data files.
 *   Copyright (C) 2020-2021 Michał Startek and Mateusz Łącki
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License, version 3 only,
 *   as published by the Free Software Foundation.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <memory>
#include "bruker_api.h"
#include "platform.h"

#ifdef OPENTIMS_BUILDING_R
#include "opentimsr_types.h"
#else
#include "opentims.h"
#endif

#include "so_manager.h"

class Scan2InvIonMobilityConverter
{
 public:
    virtual void convert(uint32_t frame_id,
                         double* inv_ion_mobilities,
                         const double* scans,
                         uint32_t size) = 0;

    virtual void convert(uint32_t frame_id,
                         double* inv_ion_mobilities,
                         const uint32_t* scans,
                         uint32_t size) = 0;

    virtual void inverse_convert(uint32_t frame_id,
                                 uint32_t* scans,
                                 const double* inv_ion_mobilities,
                                 uint32_t size) = 0;

    virtual ~Scan2InvIonMobilityConverter();
    virtual std::string description() const;
};

class ErrorScan2InvIonMobilityConverter : public Scan2InvIonMobilityConverter
{
 public:
    ErrorScan2InvIonMobilityConverter(TimsDataHandle&);
    void convert(uint32_t, double*, const double*, uint32_t) override final;
    void convert(uint32_t, double*, const uint32_t*, uint32_t) override final;
    void inverse_convert(uint32_t frame_id, uint32_t* scans, const double* inv_ion_mobilities, uint32_t size) override final;
    std::string description() const override final;
};

class BrukerScan2InvIonMobilityConverter final : public Scan2InvIonMobilityConverter
{
    const LoadedLibraryHandle lib_handle;
    uint64_t bruker_file_handle;

    tims_open_fun_t* tims_open;
    tims_get_last_error_string_fun_t* tims_get_last_error_string;
    tims_close_fun_t* tims_close;
    tims_convert_fun_t* tims_scannum_to_inv_ion_mobility;
    tims_convert_fun_t* tims_inv_ion_mobility_to_scannum;

    std::string get_tims_error();

 public:
    BrukerScan2InvIonMobilityConverter(TimsDataHandle& TDH, const std::string& lib_path);
    ~BrukerScan2InvIonMobilityConverter();

    void convert(uint32_t frame_id, double* inv_ion_mobilities, const double* scans, uint32_t size) override final;
    void convert(uint32_t frame_id, double* inv_ion_mobilities, const uint32_t* scans, uint32_t size) override final;
    void inverse_convert(uint32_t frame_id, uint32_t* scans, const double* inv_ion_mobilities, uint32_t size) override final;

    std::string description() const override final;
};


/* ============================================================================================= */

class Scan2InvIonMobilityConverterFactory
{
 public:
    virtual std::unique_ptr<Scan2InvIonMobilityConverter> produce(TimsDataHandle& TDH) = 0;
    virtual ~Scan2InvIonMobilityConverterFactory();
};

class ErrorScan2InvIonMobilityConverterFactory final : public Scan2InvIonMobilityConverterFactory
{
 public:
    std::unique_ptr<Scan2InvIonMobilityConverter> produce(TimsDataHandle& TDH) override final;
};

class BrukerScan2InvIonMobilityConverterFactory final : public Scan2InvIonMobilityConverterFactory
{
    const std::string dll_path;
    const LoadedLibraryHandle lib_hndl;
 public:
    BrukerScan2InvIonMobilityConverterFactory(const char* _dll_path);
    BrukerScan2InvIonMobilityConverterFactory(const std::string& _dll_path);
    std::unique_ptr<Scan2InvIonMobilityConverter> produce(TimsDataHandle& TDH) override final;
};

class DefaultScan2InvIonMobilityConverterFactory final
{
    static std::unique_ptr<Scan2InvIonMobilityConverterFactory> fac_instance;
 public:
    static std::unique_ptr<Scan2InvIonMobilityConverter> produceDefaultConverterInstance(TimsDataHandle& TDH);

    template<class FactoryType, class... Args> static void setAsDefault(Args&& ... args)
    {
        static_assert(std::is_base_of<Scan2InvIonMobilityConverterFactory, FactoryType>::value, "FactoryType must be a subclass of Scan2InvIonMobilityConverterFactory");
        fac_instance = std::make_unique<FactoryType>(std::forward<Args...>(args...));
    }
};
