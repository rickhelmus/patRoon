/*
 *  msexception.h, part of GenForm by M. Meringer Copyright (C) 2015
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

#ifndef MS_EXCEPTION_H
#define MS_EXCEPTION_H

#include <stdexcept>
#include <string>

class MsError : public std::exception {
 private:
  std::string m_strMsg;

 public:
  MsError( const std::string& strMsg)        { m_strMsg = strMsg; }
  ~MsError() throw()                                              { }
  virtual const char* what() const throw()  { return(m_strMsg.data()); }
};

#endif // MS_EXCEPTION_H
