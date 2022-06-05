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

#include "tof2mz_converter.h"

std::unique_ptr<Tof2MzConverterFactory> DefaultTof2MzConverterFactory::fac_instance;

Tof2MzConverter::~Tof2MzConverter() {}

std::string Tof2MzConverter::description() { return "Tof2MzConverter default"; }

/*
 * ErrorTof2MzConverter implementation
 */
ErrorTof2MzConverter::ErrorTof2MzConverter(TimsDataHandle&) {}

void ErrorTof2MzConverter::convert(uint32_t, double*, const double*, uint32_t)
{
    throw std::logic_error("Default conversion method must be selected BEFORE opening any TimsDataHandles - or it must be passed explicitly to the constructor");
}

void ErrorTof2MzConverter::convert(uint32_t, double*, const uint32_t*, uint32_t)
{
    throw std::logic_error("Default conversion method must be selected BEFORE opening any TimsDataHandles - or it must be passed explicitly to the constructor");
}

void ErrorTof2MzConverter::inverse_convert(uint32_t, uint32_t*, const double*, uint32_t)
{
    throw std::logic_error("Default conversion method must be selected BEFORE opening any TimsDataHandles - or it must be passed explicitly to the constructor");
}

std::string ErrorTof2MzConverter::description() { return "ErrorTof2MzConverter default"; }

/*
 * BrukerTof2MzConverter implementation
 */

std::string BrukerTof2MzConverter::get_tims_error()
{
    const size_t buf_size = 10000;
    std::unique_ptr<char[]> buf = std::make_unique<char[]>(buf_size);
    tims_get_last_error_string(buf.get(), buf_size-1);
    buf[buf_size-1] = '\0';
    return std::string(buf.get());
}

BrukerTof2MzConverter::BrukerTof2MzConverter(TimsDataHandle& TDH, const std::string& lib_path) : lib_handle(lib_path), bruker_file_handle(0)
{
    tims_open = lib_handle.symbol_lookup<tims_open_fun_t>("tims_open");
    tims_get_last_error_string = lib_handle.symbol_lookup<tims_get_last_error_string_fun_t>("tims_get_last_error_string");
    tims_close = lib_handle.symbol_lookup<tims_close_fun_t>("tims_close");
    tims_index_to_mz = lib_handle.symbol_lookup<tims_convert_fun_t>("tims_index_to_mz");
    tims_mz_to_index = lib_handle.symbol_lookup<tims_convert_fun_t>("tims_mz_to_index");

    bruker_file_handle = (*tims_open)(TDH.tims_dir_path.c_str(), 0); // Recalibrated states not supported

    if(bruker_file_handle == 0)
        throw std::runtime_error("tims_open(" + TDH.tims_dir_path + ") failed. Reason: " + get_tims_error());
}

BrukerTof2MzConverter::~BrukerTof2MzConverter()
{
    if(bruker_file_handle != 0) tims_close(bruker_file_handle);
}

void BrukerTof2MzConverter::convert(uint32_t frame_id, double* mzs, const double* tofs, uint32_t size)
{
    tims_index_to_mz(bruker_file_handle, frame_id, tofs, mzs, size);
}

void BrukerTof2MzConverter::convert(uint32_t frame_id, double* mzs, const uint32_t* tofs, uint32_t size)
{
    std::unique_ptr<double[]> dbl_tofs = std::make_unique<double[]>(size);
    for(uint32_t idx = 0; idx < size; idx++)
        dbl_tofs[idx] = static_cast<double>(tofs[idx]);
    tims_index_to_mz(bruker_file_handle, frame_id, dbl_tofs.get(), mzs, size);
}

void BrukerTof2MzConverter::inverse_convert(uint32_t frame_id, uint32_t* tofs, const double* mzs, uint32_t size)
{
    std::unique_ptr<double[]> dbl_tofs = std::make_unique<double[]>(size);
    tims_mz_to_index(bruker_file_handle, frame_id, mzs, dbl_tofs.get(), size);
    for(uint32_t idx = 0; idx < size; idx++)
        tofs[idx] = static_cast<uint32_t>(dbl_tofs[idx]);
}

std::string BrukerTof2MzConverter::description()
{
    return "BrukerTof2MzConverter";
}

/*
 * Tof2MzConverterFactory implementation
 */

Tof2MzConverterFactory::~Tof2MzConverterFactory() {}

/*
 * ErrorTof2MzConverterFactory
 */
std::unique_ptr<Tof2MzConverter> ErrorTof2MzConverterFactory::produce(TimsDataHandle& TDH)
{
    return std::make_unique<ErrorTof2MzConverter>(TDH);
}

/*
 * BrukerTof2MzConverterFactory implementation
 */

BrukerTof2MzConverterFactory::BrukerTof2MzConverterFactory(const char* _dll_path) :
    dll_path(_dll_path),
    lib_hndl(_dll_path)
    {}

BrukerTof2MzConverterFactory::BrukerTof2MzConverterFactory(const std::string& _dll_path) :
    dll_path(_dll_path),
    lib_hndl(_dll_path)
    {}

std::unique_ptr<Tof2MzConverter> BrukerTof2MzConverterFactory::produce(TimsDataHandle& TDH)
{
    return std::make_unique<BrukerTof2MzConverter>(TDH, dll_path.c_str());
}

/*
 * DefaultTof2MzConverterFactory implementation
 */

std::unique_ptr<Tof2MzConverter> DefaultTof2MzConverterFactory::produceDefaultConverterInstance(TimsDataHandle& TDH)
{
    if(!fac_instance)
        fac_instance = std::make_unique<ErrorTof2MzConverterFactory>();

    return fac_instance->produce(TDH);
}
