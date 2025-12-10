//
// Created by lucile.schulthe on 14.02.2022.
//

#ifndef MAIN_CPP_NETWORK_H
#define MAIN_CPP_NETWORK_H

#include <string>
#include "../iofile/iofile.h"
#include "../math/maths.h"

#include "../BasicEnergySystem/BasicES.h"
#include "../state/TState.h"
// #include "../BasicEnergySystem/TBranch.h"
#include "../consumer/sst.h"
#include "../fluids/FluidCore.h"
#include "TBranch.h"

class sst;


class Tnetwork {
public:
    Tnetwork(){};
    ~Tnetwork(){};
    Tnetwork(double &TDHN, double &dT, double &ResRatio,int &pressure,int &alt_plant,int &alt_max,int &alt_min, char &fluidin1);
    Tnetwork(double &TDHN, double &dT, double &ResRatio,int &pressure,int &alt_plant,int &alt_max,int &alt_min, string &fluidin1);
    void init_fluid( char &fluidin1);
    void init_fluid( string &fluidin1);
    void branch_init(char* &file_in);

    void Tcalcul(){};

    void getpower(){};


    void   init (environment & ext);


    // Fonction non implementée mais pour avoir l'héritage pour l'instant

 //   BasicES & performance();
 //   float objective(TSequence& seq);
 //   void display();



    FluidCore* fluid;

    double dP=0,dT=0,ResRatio=0;
    int alt_plant=0,alt_max=0, alt_min=0;

    double energy=0;
    double Power;
    int Z=0,A=0,N=0,B=0;

   std::vector<double>power;

    Flow_physic *flow;

   vector<sst> e; // sst des réseaux (tout type)
    vector<HPipe> pipe; // list des pipes


    vector<TBranch> branch; // all the data are put inside the branches
    string Ctype,pipetype,insultype;
};
#endif //MAIN_CPP_NETWORK_H
