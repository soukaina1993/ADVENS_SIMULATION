//
// Created by lucile.schulthess on 13.08.2021.
//

#include <cmath>
#pragma once

#include <string>
#include "../iofile/iofile.h"
#include <vector>

/****************************************************************************/

class norm_energy {
public:
    explicit norm_energy(const char* ifile_in= "dataBase/EnergieBat_norm.txt",unsigned int inTline=1, unsigned int inTcol=1, char icolDelimiter='\t');

    void init(){
        energyconsum();
    };

    void heatconsum();
    void energyconsum();

public:
    vector<vector<double>>	Qh,Qc,Qv,Qe,R,Thref;
    vector<double>	Sur,AthAe,Tref,Tdem,Teref,Qw,Mw;
    int cat=4,cat_ren=4,nType=12;

private:
    char* file_in ;
    char sep;
    unsigned int nbTitleLines,nbTitleColumns;
};

/****************************************************************************/

class norm_ecs {
public:

    norm_ecs(const char* ifile_in= "dataBase/ECS.txt",unsigned int inTline = 1, unsigned int inTcol = 1, const char icolDelimiter='\t'){
        file_in=(char*)ifile_in;
        sep = icolDelimiter;
        nbTitleLines = inTline;
        nbTitleColumns = inTcol;
        init();}

    void init();

    // --------------------- DATA -----------------//
    vector<vector<double>>	profil_ECS;
    int hour=24, cat=12;

private:
    char* file_in;
    char sep;
    unsigned int nbTitleLines,nbTitleColumns;
};


/****************************************************************************/
const norm_energy Q_norm;    //global variable

class BuildingNorm{
public:
    BuildingNorm(){};
    BuildingNorm(double iQh);
    BuildingNorm(const char* file_in,const char* EGID, unsigned int nTline=1, unsigned int nTcol=1, char isep = '\t');
    BuildingNorm(const char* file_in,const int index=1, unsigned int nTline=1, unsigned int nTcol=1, char isep = '\t');
    void List(TSeqsArray& Records,const int index);
    void getB(TSequence* PtrSeq);
    void getBuildingName(TSeqsArray Fields,TSequence* PtrSeq,int nbFields);
    void init();
    void egid_init(const char* iegid);    //EGID = eidgenoessischer Gebaeudeindikator
    void Display();


// using norm to calculate the need
    void heat_norm();
    void consum_norm();
    void consum_norm(bool renov);
    void water_norm();
    unsigned int findID(TSeqsArray &Fields,char * IDname, int nFields);

    // avec correction


    //string sst_type;
    //double dT;

    double Qh=0,Qc=0,Qe=0,Qw=0,Thref=0,Teref=0,Tdem=0;
    double surf_sol=0,surf_pl=0,vol=0,SRE=0,Qv=0,heq=0,ren=0,Mw=0;
    int year=0, nb_stage=0, type=0,nb_user=0,ren_etat=0,cat=0;

private:
    char* file_in;
    char sep;
    unsigned int nbTitleLines,nbTitleColumns;
};