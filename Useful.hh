/*===================================================================================

  File: Useful.hh

  ===================================================================================

    Copyright (C) 2016 Vicente Palazón-González

    This program is free software; you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation; either version 3 of the License, or (at your option) any later
    version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this
    program (file LICENSE.txt); if not, write to the Free Software Foundation, Inc.,
    675 Mass Ave, Cambridge, MA 02139, USA.

    You can contact the author at:

    palazon@uji.es
    
  =================================================================================== */

#ifndef _USEFUL_HH
#define _USEFUL_HH 

#include <sys/stat.h>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

bool fileExists(std::string strFilename);

class BadConversion : public std::runtime_error {
public:
    BadConversion(const std::string& s)
        : std::runtime_error(s)
    { }
};
 
float convertToFloat(const std::string& s);

double convertToDouble(const std::string& s);

float myatof(char *s);

#endif

