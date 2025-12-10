//
// Created by lucile.schulthess on 13.08.2021.
//

#include <cmath>
#pragma once

#include "../math/maths.h"
#include <string>
#include "../iofile/iofile.h"
#include <vector>

/****************************************************************************/


/****************************************************************************/

class environment {
public:
    environment(const char* ifile_in, vector<double> &data,const char* dataName, char isep = '\t', int inTline = 1, int inTcol = 1);
    environment(const char* ifile_in,vector<double>	&data, int index=1, char isep = '\t', int inTline = 1, int inTcol = 1);
    environment()=default;
    ~environment()=default;
    environment(const char* city, const char* year);

    void Moredata(const char* file_in,vector<double>	&data, int index);
    void Moredata(vector<double> &data,const char* city, const char* year, const char* type);
    void Moredata(vector<double> &data, int index);
    void Moredata(const char* file_in,vector<double>	&data,const char *dataName);
    void Moredata(vector<double>	&data,const char *dataName);


 //   environment(const char* file_in,const char *dataName);
 //   environment(const char* file_in,const int index);

    void init(){};
    void init(const char* city, const char* year);
    void init(const char* ifile_in, char isep = '\t', int inTline = 1, int inTcol = 1);
    void init(const char* ifile_in, vector<double> &data, const char *dataName, char isep = '\t', int inTline = 1, int inTcol = 1);
    void init(const char* ifile_in, vector<double>&data, int index=1, char isep = '\t', int inTline = 1, int inTcol = 1);


    void Display();

    vector<double>	T,Irr,hum,flow,Irr_diff,Tsol; // est-ce utile ?
    vector<int>	a;
    int Time=0;
    double Tmean=0;

private:
    //const char* file_in;
    string file_in;
    char sep;
    int nbTitleLines,nbTitleColumns;
};
