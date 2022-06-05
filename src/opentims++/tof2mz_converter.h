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

class Tof2MzConverter
{
 public:
    virtual void convert(uint32_t frame_id, double* mzs, const double* tofs, uint32_t size) = 0;
    virtual void convert(uint32_t frame_id, double* mzs, const uint32_t* tofs, uint32_t size) = 0;
    virtual void inverse_convert(uint32_t frame_id, uint32_t* tofs, const double* mzs, uint32_t size) = 0;
    virtual ~Tof2MzConverter();
    virtual std::string description();
};

class ErrorTof2MzConverter : public Tof2MzConverter
{
 public:
    ErrorTof2MzConverter(TimsDataHandle&);
    void convert(uint32_t, double*, const double*, uint32_t) override final;
    void convert(uint32_t, double*, const uint32_t*, uint32_t) override final;
    void inverse_convert(uint32_t frame_id, uint32_t* tofs, const double* mzs, uint32_t size) override final;
    std::string description() override final;
};

class BrukerTof2MzConverter final : public Tof2MzConverter
{
    const LoadedLibraryHandle lib_handle;
    uint64_t bruker_file_handle;

    tims_open_fun_t* tims_open;
    tims_get_last_error_string_fun_t* tims_get_last_error_string;
    tims_close_fun_t* tims_close;
    tims_convert_fun_t* tims_index_to_mz;
    tims_convert_fun_t* tims_mz_to_index;

    std::string get_tims_error();

 public:
    BrukerTof2MzConverter(TimsDataHandle& TDH, const std::string& lib_path);
    ~BrukerTof2MzConverter();

    void convert(uint32_t frame_id, double* mzs, const double* tofs, uint32_t size) override final;
    void convert(uint32_t frame_id, double* mzs, const uint32_t* tofs, uint32_t size) override final;
    void inverse_convert(uint32_t frame_id, uint32_t* tofs, const double* mzs, uint32_t size) override final;

    std::string description() override final;
};

class Tof2MzConverterFactory
{
 public:
    virtual std::unique_ptr<Tof2MzConverter> produce(TimsDataHandle& TDH) = 0;
    virtual ~Tof2MzConverterFactory();
};

class ErrorTof2MzConverterFactory final : public Tof2MzConverterFactory
{
 public:
    std::unique_ptr<Tof2MzConverter> produce(TimsDataHandle& TDH) override final;
};

class BrukerTof2MzConverterFactory final : public Tof2MzConverterFactory
{
    const std::string dll_path;
    const LoadedLibraryHandle lib_hndl;
 public:
    BrukerTof2MzConverterFactory(const char* _dll_path);
    BrukerTof2MzConverterFactory(const std::string& _dll_path);
    std::unique_ptr<Tof2MzConverter> produce(TimsDataHandle& TDH) override final;
};

class DefaultTof2MzConverterFactory final
{
    static std::unique_ptr<Tof2MzConverterFactory> fac_instance;
 public:
    static std::unique_ptr<Tof2MzConverter> produceDefaultConverterInstance(TimsDataHandle& TDH);

    template<class FactoryType, class... Args> static void setAsDefault(Args&& ... args)
    {
        static_assert(std::is_base_of<Tof2MzConverterFactory, FactoryType>::value, "FactoryType must be a subclass of Tof2MzConverterFactory");
        fac_instance = std::make_unique<FactoryType>(std::forward<Args...>(args...));
    }
};
