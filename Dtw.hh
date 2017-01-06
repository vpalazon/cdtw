/*===================================================================================

  File: Dtw.hh

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

#ifndef _DTW_HH
#define _DTW_HH

#include "Samples.hh"
#include <vector>
#include <list>
#include <utility>
#include <cmath>
#include "matrix.hh"

#include "Defs.hh"
#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>

template <typename T> float dtw(Sample<T> &sx, Sample<T> &sy, float (*d) (T, T));

template <typename T> float bruteForceCycDtw(Sample<T> &sx, Sample<T> &sy, float (*d) (T, T));

template <typename T> float maes(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T));

float xdtw(unsigned lenx, unsigned leny, unsigned left, unsigned right, 
	   matrix<float> &D, const matrix<float> &d, matrix<unsigned> &mins, matrix<unsigned> &maxs);
float maesStep(unsigned lenx, unsigned leny, unsigned left, unsigned right, 
	       matrix<float> &D, const matrix<float> &d, matrix<unsigned> &mins, matrix<unsigned> &maxs);

template <typename T>
float dtw(Sample<T> &sx, Sample<T> &sy, float (*d) (T, T)) {
    vector<T> &x = sx.values;
    vector<T> &y = sy.values;
    matrix<float> D = matrix<float>(x.size(), y.size());
    float distance;

    D[0][0] = d(x[0], y[0]);

    for (unsigned i=1; i<x.size(); ++i) 
        D[i][0] = D[i-1][0] + d(x[i], y[0]);

    for (unsigned j=1; j<y.size(); ++j) 
        D[0][j] = D[0][j-1] + d(x[0], y[j]);

    for (unsigned i=1; i<x.size(); ++i)
        for (unsigned j=1; j<y.size(); ++j) {
            float a=D[i-1][j-1], b = D[i-1][j], c = D[i][j-1];
            distance = d(x[i], y[j]);
            if (a <= b)
                if (a <= c) D[i][j] = a+distance;
                else D[i][j] = c+distance;
            else
                if (b <= c) D[i][j] = b+distance;
                else D[i][j] = c+distance;
        }

    return D[x.size()-1][y.size()-1];
}

template <typename T>
float bruteForceCycDtw(Sample<T> &sx, Sample<T> &sy, float (*d) (T, T)) {
    vector<T> &x = sx.values;
    vector<T> &y = sy.values;
    vector<float> distances;
    distances.reserve(x.size()*y.size());
    float distance;
    for (unsigned i=0; i<x.size(); i++) {
        x.insert(x.begin(), *(x.end()-1));
        x.erase(x.end()-1);
        for (unsigned int j=0; j<y.size(); j++) { 
            y.insert(y.begin(), *(y.end()-1));
            y.erase(y.end()-1);
            distance = dtw(sx, sy, d);
            distances.push_back(distance);
        }
    }
    return *min_element(distances.begin(), distances.end());
}


template <typename T>
float maes(Sample<T> &sx, Sample<T> &sy, float (*distance) (T, T)) {
    vector<T> &x = sx.values;
    vector<T> &y = sy.values;

    unsigned lenx = x.size(), leny = y.size();
    matrix<float> d(2*lenx+1, leny);
    for (unsigned i = 0; i < lenx; i++)
        for (unsigned j = 0; j < leny; j++)
            d[i+lenx][j] = d[i][j] = distance(x[i], y[j]);
  
    for (unsigned j=0; j<leny; j++)
        d[2*lenx][j] = d[0][j];
  
    matrix<float> D(2*lenx+1, leny);
    D[0][0] = d[0][0];
    for (unsigned i=1; i<lenx+1; ++i) D[i][0] = D[i-1][0] + d[i][0];
    for (unsigned j=1; j<leny; ++j) D[0][j] = D[0][j-1] + d[0][j];
    for (unsigned i=1; i<lenx+1; ++i)
        for (unsigned j=1; j<leny; ++j) {
            float a=D[i-1][j-1], b = D[i-1][j], c = D[i][j-1];
            if (a <= b)
                if (a <= c) D[i][j] = a+d[i][j];
                else D[i][j] = c+d[i][j];
            else
                if (b <= c) D[i][j] = b+d[i][j];
                else D[i][j] = c+d[i][j];
        }
    float dist = min(D[lenx-1][leny-1], D[lenx][leny-1]);
    matrix<unsigned> mins(lenx+1, leny);
    matrix<unsigned> maxs(lenx+1, leny);
    for (unsigned i=0; i<=lenx; ++i)
        for (unsigned j=0; j<leny; ++j) {
            mins[i][j] = 2*lenx;
            maxs[i][j] = 0;
        }
    unsigned i = lenx-1, j = leny-1;
    mins[0][j] = maxs[0][j] = i;
    while (i and j) {
        float a=D[i-1][j-1], b = D[i-1][j], c = D[i][j-1];
        if (a <= b)
            if (a <= c) --i,--j;
            else --j;
        else
            if (b <= c) --i;
            else --j;
        if (i < mins[0][j]) mins[0][j] = i;
        if (i > maxs[0][j]) maxs[0][j] = i;
    }
    for (int k=j; k>=0; --k) {
        mins[0][k] = 0;
        maxs[0][k] = max(maxs[0][k], i);
    }
    for (unsigned k=0; k<leny; k++) {
        mins[lenx][k] = mins[0][k] + lenx;
        maxs[lenx][k] = maxs[0][k] + lenx;
    }
    float r =  maesStep(lenx, leny, 0, lenx, D, d, mins, maxs);
    return min(dist, r);
}


#endif

