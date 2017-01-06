/*===================================================================================

  File: Samples.hh

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

#ifndef _SAMPLES_HH
#define _SAMPLES_HH

#include <vector>

#include <string>
#include <fstream>
#include <strstream>
#include <iostream>
#include <cstdlib>

#include "Useful.hh"

using namespace std;

#define MAX_LINE_SZ 200000 // Tamaño máximo de línea del fichero (en caracteres).
#define MAX_VALUE_SZ 100 	// Tamaño máximo de valor

typedef float t_value;

typedef pair<float, float> t_pvalue;

typedef vector<float> t_vvalue;

template <typename T> class Sample {
public:
    int label;
    vector<T> values;
    void init(int label, vector<T> &v);
    unsigned size(void);
    T &operator [] (unsigned i);
    vector<T> *getValues(void);
    int getLabel(void) {return label;}
    void cyclicShift(int r);
    void clear(void);
    
private:
    void initT(int label, vector<t_value> &v);
    void initT(int label, vector<t_pvalue> &v);
    void initT(int label, vector<t_vvalue> &v);
};

template <typename T> class Samples {
public:
    vector< Sample <T> > *samples;
    unsigned vector_size;
    Samples();
    ~Samples();
    void setVectorSize(unsigned size);
    void read(const char *filename);
    void push_back(Sample<T> &s);
    void show(void);
    Sample<T> *get(unsigned i);
    unsigned size(void);
private:
    void readT(vector< Sample<t_value> > *p, const char *filename);
    void readT(vector< Sample<t_pvalue> > *p, const char *filename);
    void readT(vector< Sample<t_vvalue> > *p, const char *filename);
    void showT(vector< Sample<t_value> > *p);
    void showT(vector< Sample<t_pvalue> > *p);
    void showT(vector< Sample<t_vvalue> > *p);
};

template <typename T>
void Sample<T>::init(int label, vector<T> &v) {
    initT(label, v);
}

template <typename T>
void Sample<T>::initT(int label, vector<t_value> &v) {
    this->label = label;
    this->values.resize(v.size());
    for (unsigned i = 0; i < v.size(); i++)
        this->values[i] = v[i];
}

template <typename T>
void Sample<T>::initT(int label, vector<t_pvalue> &v) {
    this->label = label;
    this->values.resize(v.size());
    for (unsigned i = 0; i < v.size(); i++) {
    this->values[i].first = v[i].first;
    this->values[i].second = v[i].second;
  }
}

template <typename T>
void Sample<T>::initT(int label, vector<t_vvalue> &v) {
    this->label = label;
    this->values.resize(v.size());
    for (unsigned i = 0; i < v.size(); i++) {
        this->values[i].resize(v[i].size());
        for (unsigned j = 0; j < v[i].size(); j++)
            this->values[i][j] = v[i][j];
    }
}

template <typename T>
unsigned Sample<T>::size() {
    return this->values.size();
}

template <typename T>
T &Sample<T>::operator [] (unsigned i) {
    return values[i];
}

template <typename T>
vector<T> *Sample<T>::getValues(void) {
    return &(this->values);
}

template <typename T>
void Sample<T>::cyclicShift(int r) {
    if (r != 0) {
        for (int i=0; i<r; i++) {
            values.insert(values.end(), *(values.begin()));
            values.erase(values.begin());
        }
    }
}

template <typename T>
void Sample<T>::clear(void) {
    values.clear();
}

template <typename T>
Samples<T>::Samples() {
    this->samples = 0;
    this->vector_size = 0;
}

template <typename T>
Samples<T>::~Samples() {
    delete this->samples;
}

template <typename T>
void Samples<T>::setVectorSize(unsigned size) {
    this->vector_size = size;
}

template <typename T>
void Samples<T>::read(const char *filename) {
    readT(this->samples, filename);
}

template <typename T>
void Samples<T>::readT(vector< Sample<t_value> > *p, const char *filename) {
    char buf[MAX_LINE_SZ];
    string aux;
    string svalue;
    fstream file_op(filename, ios::in);
    Sample<t_value> s;
    vector< Sample<t_value> > *ss = new vector< Sample<t_value> >();
    char label[MAX_VALUE_SZ];
    char value[MAX_VALUE_SZ];

    while(!file_op.eof()) {
        file_op.getline(buf, MAX_LINE_SZ);
        svalue = buf;
        if (svalue == "") continue; // la última línea puede contener sólo un '\n'
        istrstream isvalue(svalue.c_str(), svalue.size());
        isvalue.get(label, MAX_VALUE_SZ, ' ');
        s.label = atoi(label);
        isvalue.get();		// quitamos el ' '
        while(isvalue.get(value, MAX_VALUE_SZ, ' ')) {
            isvalue.get();
            //s.values.push_back(atof(value));
            s.values.push_back(myatof(value));
        }
        
        ss->push_back(s);
        s.values.erase(s.values.begin(), s.values.end());
    }
    
    file_op.close();
    this->samples = ss;
}

template <typename T>
void Samples<T>::readT(vector< Sample<t_pvalue> > *p, const char *filename) {
    char buf[MAX_LINE_SZ];
    string aux;
    string svalue;
    fstream file_op(filename, ios::in);
    Sample<t_pvalue> ps;

    vector< Sample<t_pvalue> > *pss = new vector< Sample<t_pvalue> >();
    char label[MAX_VALUE_SZ];
    char value[MAX_VALUE_SZ];
    float f1;
    float f2;

    while(!file_op.eof()) {
        file_op.getline(buf, MAX_LINE_SZ);
        svalue = buf;
        if (svalue == "") continue; // la última línea puede contener sólo un '\n'
        istrstream isvalue(svalue.c_str(), svalue.size());
        isvalue.get(label, MAX_VALUE_SZ, ' ');
        ps.label = atoi(label);
        isvalue.get();		// quitamos el ' '
        while(isvalue.get(value, MAX_VALUE_SZ, ' ')) {
            isvalue.get();
            //f1 = atof(value);
            f1 = myatof(value);
            isvalue.get(value, MAX_VALUE_SZ, ' ');
            isvalue.get();
            //f2 = atof(value);
            f2 = myatof(value);
            ps.values.push_back(pair<float, float>(f1, f2));
        }

        pss->push_back(ps);
        ps.values.erase(ps.values.begin(), ps.values.end());
    }

    file_op.close();
    this->samples = pss;
}

template <typename T>
void Samples<T>::readT(vector< Sample<t_vvalue> > *p, const char *filename) {
    if (this->vector_size == 0) {
        cerr << "Error: vector size not initialized\n" << endl;
        return;
    }

    char buf[MAX_LINE_SZ];
    string aux;
    string svalue;
    fstream file_op(filename, ios::in);
    Sample<t_vvalue> vs;
    t_vvalue vaux;
    t_vvalue v(this->vector_size);
    vector< Sample<t_vvalue> > *vss = new vector< Sample<t_vvalue> >();
    char label[MAX_VALUE_SZ];
    char value[MAX_VALUE_SZ];

    // unsigned num = 0;

    while(!file_op.eof()) {
        file_op.getline(buf, MAX_LINE_SZ);
        svalue = buf;
        if (svalue == "") continue; // la última línea puede contener sólo un '\n'
        istrstream isvalue(svalue.c_str(), svalue.size());
        isvalue.get(label, MAX_VALUE_SZ, ' ');
        vs.label = atoi(label);
        isvalue.get();		// quitamos el ' '
        while(isvalue.get(value, MAX_VALUE_SZ, ' ')) {
            isvalue.get();
            //vaux.push_back(atof(value));
            vaux.push_back(myatof(value));
        }

        for (unsigned i=0; i<vaux.size()/this->vector_size; i++) {
            for (unsigned j=0; j<this->vector_size; j++) {
                v[j] = vaux[i*this->vector_size+j];
            }
    
            vs.values.push_back(v);
        }
        vss->push_back(vs);
        vaux.erase(vaux.begin(), vaux.end());
        vs.values.erase(vs.values.begin(), vs.values.end());

        // cout << num << endl;
        // num++;
        // if (num == 100)
        //     break;
    }

    file_op.close();
    this->samples = vss;
}

template <typename T>
void Samples<T>::push_back(Sample<T> &s) {
    if (this->samples == 0)
        this->samples = new vector< Sample<T> >();
    this->samples->push_back(s);
}

template <typename T>
void Samples<T>::show(void) {
    showT(this->samples);
}

template <typename T>
void Samples<T>::showT(vector< Sample<t_value> > *p) {
    for(unsigned int i=0; i<p->size(); i++) {
        cout << (*p)[i].label << " ";
        for(unsigned int j=0; j<(*p)[i].values.size(); j++) {
            cout << (*p)[i].values[j] << " ";
        }
        cout << endl;
    }
}

template <typename T>
void Samples<T>::showT(vector< Sample<t_pvalue> > *p) {
    for(unsigned int i=0; i<p->size(); i++) {
        cout << (*p)[i].label << " ";
        for(unsigned int j=0; j<(*p)[i].values.size(); j++) {
            cout << (*p)[i].values[j].first
                 << ", " 
                 << (*p)[i].values[j].second
                 << " ";
        }
        cout << endl;
    }
}

template <typename T>
void Samples<T>::showT(vector< Sample<t_vvalue> > *p) {
    for(unsigned int i=0; i<p->size(); i++) {
        cout << (*p)[i].label << " ";
        for(unsigned int j=0; j<(*p)[i].values.size(); j++) {
            for (unsigned int k=0; k<(*p)[i].values[j].size(); k++) {
                cout << (*p)[i].values[j][k] << " ";
            }
        }
        cout << endl;
    }
}

template <typename T>
Sample<T> *Samples<T>::get(unsigned i) {
    return &((*this->samples)[i]);
}

template <typename T>
unsigned Samples<T>::size(void) {
    return this->samples->size();
}

#endif

