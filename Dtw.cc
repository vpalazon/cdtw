/*===================================================================================

  File: Dtw.cc

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

#include "Dtw.hh"
#include <cmath>

#include "Samples.hh"
#include "matrix.hh"
#include "LocalDistances.hh"

float xdtw(unsigned lenx, unsigned leny, unsigned left, unsigned right, 
	   matrix<float> &D, const matrix<float> &d, matrix<unsigned> &mins, matrix<unsigned> &maxs) {
    unsigned s = (left+right)/2;
    vector<unsigned> &lmin = mins[left], &lmax = maxs[left], &rmin = mins[right], &rmax = maxs[right];
    D[s][0] = d[s][0];
    unsigned fin = min(s+lenx+1, rmax[0]+1);
    for (unsigned i=s+1; i<fin; ++i) D[i][0] = D[i-1][0] + d[i][0];
    for (unsigned j=1; j<leny; ++j) {
        unsigned ini = max(lmin[j], s), fin = min(s+lenx+1, rmax[j]+1);
        D[ini-1][j] = INF;
        if (ini == s or lmax[j-1] == lmin[j]) D[ini-1][j-1] = INF; 
        if (fin > rmin[j]) {
            if (rmax[j-1] < rmin[j]) D[rmin[j]][j-1] = INF;
            for (unsigned i=rmin[j]+1; i<=fin; ++i) D[i][j-1] = INF;
        }
        else if (fin == rmin[j] and rmax[j-1] < fin) D[fin][j-1] = INF;
        for (unsigned i= ini; i<fin; ++i) {
            float a=D[i-1][j-1], b=D[i-1][j], c=D[i][j-1];
            if (a<=b)
                if (a<=c) D[i][j] = a + d[i][j];
                else D[i][j] = c + d[i][j];
            else
                if (b<=c) D[i][j] = b + d[i][j];
                else D[i][j] = c + d[i][j];
        }
    }
    unsigned i = s+lenx-1, j = leny-1;
    mins[s][j] = maxs[s][j] = i;
    while (i and j) {
        float a=D[i-1][j-1], b=D[i-1][j], c=D[i][j-1];
        if (a<=b)
            if (a<=c) --i, --j;
            else --j;
        else
            if (b<=c) --i;
            else --j;
        if (i < mins[s][j]) mins[s][j] = i;
        if (i > maxs[s][j]) maxs[s][j] = i;
    }
    for (int k=j; k>=0; --k) {
        mins[s][k] = s;
        maxs[s][k] = max(maxs[s][k], i);
    }
    return min(D[s+lenx-1][leny-1], D[s+lenx][leny-1]);
}

float maesStep(unsigned lenx, unsigned leny, unsigned left, unsigned right, 
	       matrix<float> &D, const matrix<float> &d, matrix<unsigned> &mins, matrix<unsigned> &maxs) {
    float r = INF;
    unsigned shift = (left+right)/2;
    if (left < shift and shift < right) {
        r = xdtw(lenx, leny, left, right, D, d, mins, maxs);
    }
    if (left < shift-1) {
        float rr = maesStep(lenx, leny, left, shift, D, d, mins, maxs);
        if (rr < r) r = rr;
    }
    if (shift < right-1) {
        float rr = maesStep(lenx, leny, shift, right, D, d, mins, maxs);
        if (rr < r) r = rr;
    }
    return r;
}

