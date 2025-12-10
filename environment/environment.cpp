//
// Created by Lucile Schulthess on 16.08.2021.
//

#include "environment.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include "../iofile/iofile.h"

using namespace std;

/****************************************************************************/


environment::environment(const char* ifile_in,vector<double>&data,const char *dataName,const char isep,const int inTline,const int inTcol){
    file_in=ifile_in;
    sep = isep;
    // Nombre de lignes / Colones de titres
    nbTitleLines = inTline;
    nbTitleColumns = inTcol;
    Getdata(file_in.c_str(), data,Time, dataName,nbTitleLines, nbTitleColumns, sep);
    T = data;
    Tmean=sum(T) / (int)T.size();
    init();
}

environment::environment(const char* ifile_in,vector<double>&data,const int index,const char isep,const int inTline,const int inTcol){
    file_in=ifile_in;
    sep = isep;
    // Nombre de lignes / Colones de titres
    nbTitleLines = inTline;
    nbTitleColumns = inTcol;
    Getdata(file_in.c_str(), data, Time, index, nbTitleLines, nbTitleColumns, sep);
    T = data;
    Tmean=sum(T) / (int)T.size();
    init();
}

environment::environment(const char* city, const char* year){
    string s;
    s ="dataBase/Meteo_";
    s +=  city;
    s += ".txt";
    file_in= s;
    init(file_in.c_str(),T,year);
}

void environment::init(const char* city, const char* year){
    string s;
    s ="dataBase/Meteo_";
    s +=  city;
    s += ".txt";
    file_in= s;
    init(file_in.c_str(), T, year);
}

void environment::init(const char* ifile_in, const char isep, const int inTline, const int inTcol){
    file_in = ifile_in;
    sep = isep;
    nbTitleLines = inTline;
    nbTitleColumns = inTcol;
}

void environment::init(const char* ifile_in,vector<double>&data,const char *dataName,const char isep,const int inTline,const int inTcol){
    file_in=ifile_in;
    sep = isep;
    // Nombre de lignes / Colones de titres
    nbTitleLines = inTline;
    nbTitleColumns = inTcol;
    Getdata(file_in.c_str(), data,Time, dataName,nbTitleLines, nbTitleColumns, sep);
    T = data;
    Tmean=sum(T) / (int)T.size();
}


void environment::init(const char* ifile_in,vector<double>&data,const int index,const char isep,const int inTline,const int inTcol){
    file_in=ifile_in;
    sep = isep;
    // Nombre de lignes / Colones de titres
    nbTitleLines = inTline;
    nbTitleColumns = inTcol;
    Getdata(file_in.c_str(), data, Time, index, nbTitleLines, nbTitleColumns, sep);
    T = data;
    Tmean=sum(T) / (int)T.size();
}

void environment::Moredata(const char* ifile_in,vector<double>&data,const int index){
    Getdata(ifile_in,data,Time,index,nbTitleLines,nbTitleColumns,sep);
}

void environment::Moredata(vector<double>&data, const char* city, const char* year, const char* type){
    string s;
    s ="dataBase/Meteo_";
    s +=  city;
    s += ".txt";
    file_in= s;
    init(file_in.c_str(), data, year);
};


void environment::Moredata(vector<double>&data,const int index){
    Getdata(file_in.c_str(),data,Time,index,nbTitleLines,nbTitleColumns,sep);
}

void environment::Moredata(const char* ifile_in,vector<double>&data, const char* dataName){
    Getdata(ifile_in,data,Time, dataName,nbTitleLines,nbTitleColumns,sep);
}

void environment::Moredata(vector<double>&data, const char *dataName){
    Getdata(file_in.c_str(),data,Time,dataName,nbTitleLines,nbTitleColumns,sep);
}

void environment::Display(){
    cout << endl;
    cout <<"/****************************************************************************/" << endl;
    cout<<"                         Environment                           "<<endl<<endl;
    cout<<"Number of time-steps             [-] :       :"<< T.size()<<endl<<endl;
    cout<<"Mean temperature              [degC] :       :"<<Tmean-273.15<<endl;
    cout<<"Max temperature               [degC] :       :"<<max(T)-273.15<<endl;
    cout<<"Min temperature               [degC] :       :"<<min(T)-273.15<<endl<<endl;

    if(Tsol.empty()==0){
    cout<<"Mean ground temperature       [degC] :       :"<<mean(Tsol)-273.15<<endl;
    cout<<"Max ground temperature        [degC] :       :"<<max(Tsol)-273.15<<endl;
    cout<<"Min ground temperature        [degC] :       :"<<min(Tsol)-273.15<<endl<<endl;}

    if(Irr.empty()==0){
        cout<<"Mean radiation            [W/m2] :       :"<<mean(Irr)<<endl;
        cout<<"Max radiation             [W/m2] :       :"<<max(Irr)<<endl;
        cout<<"Min radiation             [W/m2] :       :"<<min(Irr)<<endl;}
    cout << endl;
    if(hum.empty()==0){
        cout<<"Mean humidity             [%] :          :"<<mean(hum)<<endl;
        cout<<"Max humidity              [%] :          :"<<max(hum)<<endl;
        cout<<"Min humidity              [%] :          :"<<min(hum)<<endl;}
    cout << endl;
    if(flow.empty()==0){
        cout<<"Mean flow                 [m3/s] :       :"<<mean(flow)<<endl;
        cout<<"Max flow                  [m3/s] :       :"<<max(flow)<<endl;
        cout<<"Min flow                  [m3/s] :       :"<<min(flow)<<endl;}

    cout<<endl;
    cout <<"/****************************************************************************/" << endl;
    cout << endl;
}