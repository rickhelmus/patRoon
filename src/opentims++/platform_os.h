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

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
    #define OPENTIMS_UNIX
    #include <unistd.h>
    #ifndef _POSIX_VERSION
        #warning "Seems we're on Unix but not POSIX. Things might break."
    #endif
#elif defined(_WIN32) || defined(_WIN64)
    #define OPENTIMS_WINDOWS
#endif
