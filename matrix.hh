/*===================================================================================

  File: matrix.hh

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

#ifndef MATRIX_HH
#define MATRIX_HH

#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

template <class T>
class matrix {
public:
    matrix(unsigned size_x, unsigned size_y): array(size_x) {
        for(unsigned i=0; i<size_x; i++)
            array[i].resize(size_y);
    }
    matrix(const matrix &rhs): array(rhs.array) {}
  
    const vector<T> &operator[](int y) const
    { return array[y]; }
    vector<T> &operator[](int y)
    { return array[y]; }
  
    int size_x() const
    { return array.size(); }
  
    int size_y() const
    { return size_x() ? array[0].size() : 0; }
  
    void print() {
        for (int j=((*this).size_y()-1); j>=0; j--) {
            for (int i=0; i<(*this).size_x(); i++)
                cout << setw(2) << (*this)[i][j] << " ";
            cout << endl;
        }
    }
  
    void print2() {
        for (int i=((*this).size_x()-1); i>=0; i--) {
            for (int j=0; j<(*this).size_y(); j++)
                cout << setw(10) << scientific << setprecision(2) 
                     << (*this)[i][j] << " ";
            cout << endl;
        }
    }
  
private:
    vector< vector<T> > array;
};

#endif

