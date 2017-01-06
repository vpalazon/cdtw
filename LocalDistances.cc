/*===================================================================================

  File: LocatDistances.cc

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

#include "LocalDistances.hh"
#include "Samples.hh"
#include <cmath>
#include "Defs.hh"

float d1vsqrt(t_vvalue a, t_vvalue b) {
    float sum = 0;
    for (unsigned k=0; k<a.size(); k++)
        sum += abs(float(b[k] - a[k]));
    return sqrt(sum);
}

float d2v(t_vvalue a, t_vvalue b) {
    float sum = 0;
    float aux;
    for (unsigned k=0; k<a.size(); k++) {
        aux = b[k] - a[k];
        sum += aux*aux;
    }
    return sqrt(sum);
}

float chiSquared(t_vvalue a, t_vvalue b) {
  float sum = 0;
  float aux;
  for (unsigned k=0; k<a.size(); k++) {
    aux = a[k] - b[k];
    sum += ((aux * aux) / (a[k] + b[k] + eps));
  }
  return 0.5*sum;
}
