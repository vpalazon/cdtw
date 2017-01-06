/*===================================================================================

  File: Useful.cc

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

#include "Useful.hh"

bool fileExists(std::string strFilename) {
    struct stat stFileInfo;
    bool blnReturn;
    int intStat;

    intStat = stat(strFilename.c_str(), &stFileInfo);
    if (intStat == 0) {
        blnReturn = true;
    } else {
        blnReturn = false;
    }

    return(blnReturn);
}

float convertToFloat(const std::string& s) {
    std::istringstream i(s);
    float x;
    if (!(i >> x))
        throw BadConversion("convertToFloat(\"" + s + "\")");
    return x;
} 

double convertToDouble(const std::string& s) {
    std::istringstream i(s);
    double x;
    if (!(i >> x))
        throw BadConversion("convertToDouble(\"" + s + "\")");
    return x;
}

float myatof(char *s) {
    float aux;
    aux = convertToFloat(std::string(s));
    return aux;
}


