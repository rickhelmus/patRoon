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

typedef uint64_t tims_open_fun_t(const char *path, uint32_t recalibration);

typedef uint32_t tims_convert_fun_t(uint64_t file_hndl, int64_t frame_id, const double *tofs, double *mzs, uint32_t arr_size);

typedef uint32_t tims_scan2inv_ion_mobility_fun_t(uint64_t file_hndl, int64_t frame_id, const double *tofs, double *mzs, uint32_t arr_size);

typedef uint32_t tims_get_last_error_string_fun_t(char *target, uint32_t str_length);

typedef void tims_close_fun_t(uint64_t file_hndl);

typedef void tims_set_num_threads_t(uint32_t n);
