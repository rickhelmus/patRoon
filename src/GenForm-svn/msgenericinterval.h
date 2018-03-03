/*
 *  msgenericinterval.h, part of GenForm by M. Meringer Copyright (C) 2015
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
*/

#ifndef	MS_GENERIC_INTERVAL_H
#define MS_GENERIC_INTERVAL_H

#include <utility>

template<class T>
class GenericInterval : public std::pair<T,T>
{
public:
			GenericInterval()
			{ this->first = this->second = 0; }
			GenericInterval(T first, T second)
			{ this->first = first; this->second = second; }
			GenericInterval(T value)
			{ this->first = this->second = value;  }

			bool	IsValid()							const
			{ return(this->first<=this->second); }
};

#endif
