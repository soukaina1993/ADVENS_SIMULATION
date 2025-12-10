//
// Created by Lucile Schulthess on 16.08.2021.
//

#include "Norm_Building.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "../iofile/iofile.h"

using namespace std;

void BuildingNorm::getB(TSequence* PtrSeq){

    type=PtrSeq->Get(0)-1;
    year=PtrSeq->Get(1);
    surf_sol=PtrSeq->Get(2);
    surf_pl=PtrSeq->Get(3);
    //   SRE=PtrSeq->Get(findID(Fields,"SRE",InputFile.nFields()));
    nb_stage=PtrSeq->Get(4);
    vol=PtrSeq->Get(5);
    SRE=heq=Tdem=Teref=Thref=ren=ren_etat=nb_user=Mw=Qw=Qv=Qe=Qc=Qh=0;
/*
    Qh=PtrSeq->Get(findID(Fields,"Qh",nbFields));
    Qc=PtrSeq->Get(findID(Fields,"Qc",nbFields));
    Qv=PtrSeq->Get(findID(Fields,"Qv",nbFields));
    Qe=PtrSeq->Get(findID(Fields,"Qe",nbFields));
    Qw=PtrSeq->Get(findID(Fields,"Qw",nbFields));
    nb_user=PtrSeq->Get(findID(Fields,"nb_user",nbFields));
    ren_etat=PtrSeq->Get(findID(Fields,"ren_etat",nbFields));
    ren=PtrSeq->Get(findID(Fields,"ren",nbFields));
    Thref=PtrSeq->Get(findID(Fields,"Thref",nbFields));
    Teref=PtrSeq->Get(findID(Fields,"Teref",nbFields));
    Tdem=PtrSeq->Get(findID(Fields,"Tdem",nbFields));
    heq=PtrSeq->Get(findID(Fields,"heq",nbFields));
*/
}
void BuildingNorm::getBuildingName(TSeqsArray Fields,TSequence* PtrSeq,int nbFields){
    type=PtrSeq->Get(findID(Fields,(char*) "type",nbFields)) - 1;
    year=PtrSeq->Get(findID(Fields,(char*) "year",nbFields));
    surf_sol=PtrSeq->Get(findID(Fields,(char*) "surf_sol",nbFields));
    surf_pl=PtrSeq->Get(findID(Fields,(char*) "surf_pl",nbFields));
    //   SRE=PtrSeq->Get(findID(Fields,"SRE",nbFields));
    nb_stage=PtrSeq->Get(findID(Fields,(char*) "nb_stage",nbFields));
    vol=PtrSeq->Get(findID(Fields,(char*) "vol",nbFields));
/*
    Qh=PtrSeq->Get(findID(Fields,"Qh",nbFields));
    Qc=PtrSeq->Get(findID(Fields,"Qc",nbFields));
    Qv=PtrSeq->Get(findID(Fields,"Qv",nbFields));
    Qe=PtrSeq->Get(findID(Fields,"Qe",nbFields));
    Qw=PtrSeq->Get(findID(Fields,"Qw",nbFields));
    nb_user=PtrSeq->Get(findID(Fields,"nb_user",nbFields));
    ren_etat=PtrSeq->Get(findID(Fields,"ren_etat",nbFields));
    ren=PtrSeq->Get(findID(Fields,"ren",nbFields));
    Thref=PtrSeq->Get(findID(Fields,"Thref",nbFields));
    Teref=PtrSeq->Get(findID(Fields,"Teref",nbFields));
    Tdem=PtrSeq->Get(findID(Fields,"Tdem",nbFields));
    heq=PtrSeq->Get(findID(Fields,"heq",nbFields));
*/
}
void BuildingNorm::Display(){
    cout << endl;
    cout <<"/****************************************************************************/" << endl;
    cout<<"                         BuildingNorm                           "<<endl<<endl;
    cout<<"Heat energy total              [kwh] :       :"<<Qh<<endl;
    cout<<"ECS energy total               [kwh] :       :"<<Qw<<endl;
    cout<<"Cooling energy total           [kwh] :       :"<<Qc<<endl;
    cout<<"Ventil energy total (elec)     [kwh] :       :"<<Qv<<endl;
    cout<<"Electric energy total          [kwh] :       :"<<Qe<<endl;
    cout<<"Liter of water total           [L] :         :"<<Mw<<endl;
    cout<<"Reference Energetic Surface    [m2] :        :"<<SRE<<endl;
    cout<<"Consign temperature            [degC] :      :"<<Thref-273.15<<endl;
    cout<<"Demand temperature             [degC] :      :"<<Tdem-273.15<<endl;
    cout<<endl;
    cout <<"/****************************************************************************/" << endl;
    cout << endl;
}
/*
void BuildingNorm::List(TInputFile  &InputFile,const int index){
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    Records = *InputFile.GetRecords();
    TSequence* PtrSeq = Get(Records, index);
    getB(PtrSeq);
    init();
};*/

void BuildingNorm::List(TSeqsArray& Records, const int index) {
    TSequence *PtrSeq = Get(Records, index);
    getB(PtrSeq);
    init();
}

void  BuildingNorm::heat_norm(){
    if(Qh==0){
        Qh = Q_norm.Qh[cat][type] * SRE;}
    if (Thref == 0) {
        Thref = Q_norm.Thref[cat][type];}
}
void  BuildingNorm::consum_norm(bool renov){
    if (renov==false){
        consum_norm();
    }
    else {
        if (ren_etat == 0)
            ren = Q_norm.R[cat][type];
        if (Qh == 0)
            Qh = ((1 - ren) * Q_norm.Qh[cat][type] + ren * Q_norm.Qh[Q_norm.cat_ren][type]) * SRE;
        if (Qc == 0)
            Qc = ((1 - ren) * Q_norm.Qc[cat][type] + ren * Q_norm.Qc[Q_norm.cat_ren][type]) * SRE;
        if (Qv == 0)
            Qv = ((1 - ren) * Q_norm.Qv[cat][type] + ren * Q_norm.Qv[Q_norm.cat_ren][type]) * SRE;
        if (Qe == 0)
            Qe = ((1 - ren) * Q_norm.Qe[cat][type] + ren * Q_norm.Qe[Q_norm.cat_ren][type]) * SRE;
        if (Thref == 0)
            Thref =  Q_norm.Thref[cat][type]*(1-ren)+Q_norm.Thref[Q_norm.cat_ren][type]*ren ;
        if (Teref == 0)
            Teref = Q_norm.Teref[type];
        if (Tdem == 0)
            Tdem = Q_norm.Tdem[type];
    }
}

void  BuildingNorm::consum_norm(){
    if(Qh==0)
        Qh=Q_norm.Qh[cat][type]*SRE;
    if(Qc==0)
        Qc=Q_norm.Qc[cat][type]*SRE;
    if(Qv==0)
        Qv=Q_norm.Qv[cat][type]*SRE;
    if(Qe==0)
        Qe=Q_norm.Qe[cat][type]*SRE;
    if(Thref==0)
        Thref = Q_norm.Thref[cat][type];
    if(Teref==0)
        Teref=Q_norm.Teref[type];
    if (Tdem==0)
        Tdem=Q_norm.Tdem[type];
}

void  BuildingNorm::water_norm(){
    if(Qw==0)
        Qw=Q_norm.Qw[type]*SRE;
    if(Mw==0)
        Mw=Q_norm.Mw[type]*SRE;
}

void BuildingNorm::init(){
    if(SRE==0){
        if (surf_pl!=0)
            SRE=0.9*surf_pl;
        else if(surf_sol!=0&&nb_stage!=0){
            surf_pl=surf_sol*nb_stage;
            SRE=0.9*surf_pl;}
        else if(surf_sol!=0&&vol!=0){
            nb_stage=vol/surf_sol/2;
            surf_pl=surf_sol*nb_stage;
            SRE=0.9*surf_pl;}
        else if(nb_stage!=0&&vol!=0){
            surf_sol=vol/(nb_stage*2);
            surf_pl=surf_sol*nb_stage;
            SRE=0.9*surf_pl;}
        else if(surf_sol!=0)
            SRE=0.9*surf_sol;
        else{
            cout<<"Calcul par m2 de SRE"<<endl;
            SRE=1;}
    }
    if(year > 1000 && year <= 1980)
        cat=0;
    else if(year <=1990)
        cat=1;
    else if(year <=2000)
        cat=2;
    else if(year <=2019)
        cat=3;
    else{
        cout<<"Unknown category"<<endl;
        cat=1;}
    consum_norm();
    water_norm();
}
BuildingNorm::BuildingNorm(double iQh){
    Qh=iQh;}

BuildingNorm::BuildingNorm(const char* ifile_in,const int index, unsigned int nTline, unsigned int nTcol, char isep){
    file_in=(char*)ifile_in;
    sep = isep;
    nbTitleLines = nTline; //nombre de lignes avec des titres
    nbTitleColumns = nTcol; //nombre de colonnes avec des titres ie sans chiffre ie

    TInputFile  InputFile((char*)file_in, nbTitleLines, nbTitleColumns, sep);
    InputFile.open();
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    Records = *InputFile.GetRecords();
    TSequence* PtrSeq = Get(Records, index);
    getB(PtrSeq);
    init();
}

BuildingNorm::BuildingNorm(const char* ifile_in,const char* EGID, unsigned int nTline, unsigned int nTcol, char isep ){
    file_in=(char*)ifile_in;
    sep = isep;
    nbTitleLines = nTline; //nombre de lignes avec des titres
    nbTitleColumns = nTcol; //nombre de colonnes avec des titres ie sans chiffre ie

    egid_init(EGID);
}

void BuildingNorm::egid_init(const char* EGID){
    TInputFile  InputFile(file_in, nbTitleLines, nbTitleColumns, sep);
    InputFile.open();
    InputFile.config(Tiofile::onrcds);
    TSeqsArray Records(InputFile.nRecords());
    Records = *InputFile.GetRecords();

    TSequence* PtrSeq = Get(Records, (char*)EGID);
    getB(PtrSeq);
    init();
}

unsigned int BuildingNorm::findID(TSeqsArray &Fields,char * IDname, int nFields){
    unsigned int cpt;
    dynstring st = IDname;
    while (st != *Fields[cpt]->Title && cpt != Fields.size() - 1) {
        cpt++;}
    if (cpt>=nFields-1)
        cout<<"Not found "<<st<<endl;
    return cpt;
}


















norm_energy::norm_energy(const char* ifile_in, unsigned int nTline, unsigned int nTcol, const char colDelimiter){

    file_in=(char*)ifile_in;
    sep = colDelimiter;
    nbTitleLines = nTline;
    nbTitleColumns = nTcol;
    init();


}
void norm_energy::heatconsum(){
    TInputFile InputFile((char *) file_in, nbTitleLines, nbTitleColumns, sep);
    TSeqsArray Fields(InputFile.nFields());
    Fields = *InputFile.GetFields();
    Qh.clear();
    nType = InputFile.nRecords();
    multi_data(Fields,"Qh",cat,nType,Qh);
    multi_data(Fields,"Th",cat,nType,Thref);    //here: Thref = Th in Celsius
    one_dataT(Fields,"Tref",nType,Tref);
    one_data(Fields,"AthAe",nType,AthAe);
    for(int i(0);i<cat;i++){
        for (int j(0);j<Tref.size();j++){
            Thref[i][j] += AthAe[j]*2.5 + 0.8*(Tref[j]-293.15) + 273.15;}
    }
}

void norm_energy::energyconsum(){
    TInputFile InputFile((char *) file_in, nbTitleLines, nbTitleColumns, sep);
    TSeqsArray Fields(InputFile.nFields());
    Fields = *InputFile.GetFields();
    nType = InputFile.nRecords();

    multi_data(Fields,"Qh",cat,nType,Qh);
    multi_data(Fields,"Qc",cat,nType,Qc);
    multi_data(Fields,"Qv",cat,nType,Qv);
    multi_data(Fields,"Qe",cat,nType,Qe);

    one_data(Fields,"Qw",nType,Qw);

    multi_data(Fields,"R",cat,nType,R);

    multi_data(Fields,"Th",cat,nType,Thref);    //here: Thref = Th in Celsius
    one_dataT(Fields,"Tref",nType,Tref);
    one_data(Fields,"AthAe",nType,AthAe);
    for(int i(0);i<cat;i++){
        for (int j(0);j<nType;j++)
            Thref[i][j] += AthAe[j]*2.5 + 0.8*(Tref[j]-293.15) + 273.15;
    }
    one_dataT(Fields,"Teref",nType,Teref);
    one_dataT(Fields,"Tdem",nType,Tdem);
    one_data(Fields,"Mw",nType,Mw);
    one_data(Fields,"Sur",nType,Sur);
}
















void norm_ecs::init(){

    TInputFile InputFile((char *) file_in, nbTitleLines, nbTitleColumns, sep);
    hour = InputFile.nFields();
    cat = InputFile.nRecords();
    TSeqsArray Fields(hour);
    Fields = *InputFile.GetFields();

    vector<double> Qtemp(cat);
    vector<double> sum_Qtemp(cat,0.0);
    profil_ECS.clear();

    for (int j(1); j <= hour; ++j) {
        TSequence *PtrSeq = Get(Fields,j);

        if (PtrSeq == nullptr){
            cerr << "ECS norm cannot be read" << endl;
            exit(1);
        }
        for (int i(0); i < cat; ++i)
            Qtemp[i] = PtrSeq->Get(i);

        profil_ECS.push_back(Qtemp);
        sum_Qtemp = sum(sum_Qtemp, Qtemp);
    }

    // scale to 1
    for (int j(0); j < hour; ++j)
        for (int i(0); i < cat; ++i)
            profil_ECS[j][i] *= 1/sum_Qtemp[i];
}



/****************************************************************************/


