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

#include "scan2inv_ion_mobility_converter.h"

std::unique_ptr<Scan2InvIonMobilityConverterFactory> DefaultScan2InvIonMobilityConverterFactory::fac_instance;

Scan2InvIonMobilityConverter::~Scan2InvIonMobilityConverter() {}

std::string Scan2InvIonMobilityConverter::description() const
{
    return "Scan2InvIonMobilityConverter default";
}

/*
 * ErrorScan2InvIonMobilityConverter implementation
 */

ErrorScan2InvIonMobilityConverter::ErrorScan2InvIonMobilityConverter(TimsDataHandle&) {}

void ErrorScan2InvIonMobilityConverter::convert(uint32_t, double*, const double*, uint32_t)
{
    throw std::logic_error("Default conversion method must be selected BEFORE opening any TimsDataHandles - or it must be passed explicitly to the constructor");
}

void ErrorScan2InvIonMobilityConverter::convert(uint32_t, double*, const uint32_t*, uint32_t)
{
    throw std::logic_error("Default conversion method must be selected BEFORE opening any TimsDataHandles - or it must be passed explicitly to the constructor");
}

void ErrorScan2InvIonMobilityConverter::inverse_convert(uint32_t, uint32_t*, const double*, uint32_t)
{
    throw std::logic_error("Default conversion method must be selected BEFORE opening any TimsDataHandles - or it must be passed explicitly to the constructor");
}

std::string ErrorScan2InvIonMobilityConverter::description() const
{
    return "ErrorScan2InvIonMobilityConverter default";
}

/*
 * BrukerScan2InvIonMobilityConverter implementation
 */

std::string BrukerScan2InvIonMobilityConverter::get_tims_error()
{
    const size_t buf_size = 10000;
    std::unique_ptr<char[]> buf = std::make_unique<char[]>(buf_size);
    tims_get_last_error_string(buf.get(), buf_size-1);
    buf[buf_size-1] = '\0';
    return std::string(buf.get());
}

BrukerScan2InvIonMobilityConverter::BrukerScan2InvIonMobilityConverter(TimsDataHandle& TDH, const std::string& lib_path) : lib_handle(lib_path), bruker_file_handle(0)
{
    tims_open = lib_handle.symbol_lookup<tims_open_fun_t>("tims_open");
    tims_get_last_error_string = lib_handle.symbol_lookup<tims_get_last_error_string_fun_t>("tims_get_last_error_string");
    tims_close = lib_handle.symbol_lookup<tims_close_fun_t>("tims_close");
    tims_scannum_to_inv_ion_mobility = lib_handle.symbol_lookup<tims_convert_fun_t>("tims_scannum_to_oneoverk0");
    tims_inv_ion_mobility_to_scannum = lib_handle.symbol_lookup<tims_convert_fun_t>("tims_oneoverk0_to_scannum");

    bruker_file_handle = (*tims_open)(TDH.tims_dir_path.c_str(), 0); // Recalibrated states not supported

    if(bruker_file_handle == 0)
        throw std::runtime_error("tims_open(" + TDH.tims_dir_path + ") failed. Reason: " + get_tims_error());
}

BrukerScan2InvIonMobilityConverter::~BrukerScan2InvIonMobilityConverter()
{
    if(bruker_file_handle != 0) tims_close(bruker_file_handle);
}

void BrukerScan2InvIonMobilityConverter::convert(uint32_t frame_id,
             double* inv_ion_mobilities,
             const double* scans,
             uint32_t size)
{
    tims_scannum_to_inv_ion_mobility(bruker_file_handle, frame_id, scans, inv_ion_mobilities, size);
}


void BrukerScan2InvIonMobilityConverter::convert(uint32_t frame_id,
             double* inv_ion_mobilities,
             const uint32_t* scans,
             uint32_t size)
{
    std::unique_ptr<double[]> dbl_scans = std::make_unique<double[]>(size);
    for(uint32_t idx = 0; idx < size; idx++)
        dbl_scans[idx] = static_cast<double>(scans[idx]);
    tims_scannum_to_inv_ion_mobility(bruker_file_handle, frame_id, dbl_scans.get(), inv_ion_mobilities, size);
}

void BrukerScan2InvIonMobilityConverter::inverse_convert(uint32_t frame_id,
             uint32_t* scans,
             const double* inv_ion_mobilities,
             uint32_t size)
{
    std::unique_ptr<double[]> dbl_scans = std::make_unique<double[]>(size);
    tims_inv_ion_mobility_to_scannum(bruker_file_handle, frame_id, inv_ion_mobilities, dbl_scans.get(), size);
    for(uint32_t idx = 0; idx < size; idx++)
        scans[idx] = static_cast<double>(dbl_scans[idx]);
}

std::string BrukerScan2InvIonMobilityConverter::description() const { return "BrukerScan2InvIonMobilityConverter"; }

/*
 * Scan2InvIonMobilityConverterFactory implementation
 */

Scan2InvIonMobilityConverterFactory::~Scan2InvIonMobilityConverterFactory() {}

/*
 * ErrorScan2InvIonMobilityConverterFactory implementation
 */

std::unique_ptr<Scan2InvIonMobilityConverter> ErrorScan2InvIonMobilityConverterFactory::produce(TimsDataHandle& TDH)
{
    return std::make_unique<ErrorScan2InvIonMobilityConverter>(TDH);
}

/*
 * BrukerScan2InvIonMobilityConverterFactory implementation
 */

BrukerScan2InvIonMobilityConverterFactory::BrukerScan2InvIonMobilityConverterFactory(const char* _dll_path) : dll_path(_dll_path), lib_hndl(_dll_path) {}

BrukerScan2InvIonMobilityConverterFactory::BrukerScan2InvIonMobilityConverterFactory(const std::string& _dll_path) : dll_path(_dll_path), lib_hndl(_dll_path) {}

std::unique_ptr<Scan2InvIonMobilityConverter> BrukerScan2InvIonMobilityConverterFactory::produce(TimsDataHandle& TDH)
{
    return std::make_unique<BrukerScan2InvIonMobilityConverter>(TDH, dll_path.c_str());
}

/*
 * DefaultScan2InvIonMobilityConverterFactory implementation
 */
std::unique_ptr<Scan2InvIonMobilityConverter> DefaultScan2InvIonMobilityConverterFactory::produceDefaultConverterInstance(TimsDataHandle& TDH)
{
    if(!fac_instance)
        fac_instance = std::make_unique<ErrorScan2InvIonMobilityConverterFactory>();

    return fac_instance->produce(TDH);
}
