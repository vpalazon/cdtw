/*===================================================================================

  File: test_dtw.cc

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

#include <iostream>
#include <string>
#include <cstdlib>
#include <set>

#include "Dtw.hh"
#include "Samples.hh"
#include "LocalDistances.hh"
#include "Defs.hh"
#include "Useful.hh"

extern "C" {
#include "chronometer.h"
}

using namespace std;

int main(int argc, char **argv) {

    float (*dist) (t_vvalue, t_vvalue) = d1vsqrt;
    Samples<t_vvalue> seqs;
    seqs.setVectorSize(1);
    seqs.read("mpeg7bbr_s-cchd100_f-cur.vectors");

    ClockReset();

    for (unsigned i=0; i<seqs.samples->size(); ++i) {
        
        float min_value = INF;
        unsigned min_j = 0;
        for (unsigned j=0; j<seqs.samples->size(); ++j) {
                
            if (i == j) continue;
            
            float cost = dtw((*seqs.samples)[i],
                             (*seqs.samples)[j], dist);
            
            if (cost < min_value) {
                min_value = cost;
                min_j = j;
            }
        }

        cout << i << " " << min_j << endl;

    }
    
    float t1 = ClockTotal();

    cout << "time: " << t1 << endl;
    
    return 0;
}

