//
// Created by lucile.schulthe on 19.01.2022.
//

#ifndef MAIN_CPP_SST_H
#define MAIN_CPP_SST_H
#include "../math/maths.h"

#include <cmath>
#include "../fluids/FluidCore.h"
#include <string>
#include "../iofile/iofile.h"
#include "../environment/environment.h"
#include "solar.h"

#include "../BasicEnergySystem/heatpump/THeatpump.h"
#include "../network/TConnec.h"

#include "../network/network.h"


class TConnec;
class Tnetwork;

class sst {
public:

    sst(){};

    virtual void calcul_power(environment &ext){calcul_power();};
    virtual void calcul_power(){};

    void init_connec();
    void calcul_flow(Tnetwork &network,double & power);
    void Display();

    void input_file(vector <string> &data_input);

//    void connexion avec le reseau avec sst_type

    vector<double>power,flow_vec;

    //mu
    //,deltaT,q,power_dem;
    // peut etre seulement flow ou power? mais pas tout ou meme pas...



    double energy,dT,power_max,L_racc,DN_racc,alt,Tmin,Tpinch,dP=0;
    string sst_type;

    TConnec connec;
    Flow *flow; // flux du coté source de la sst


};







class sstfile:public sst {
public:
    sstfile(vector <string> &data_input);

    void calcul_power(); //simplification avec de l'eau à 20°C
    void Display(){};
    char* file_in;
    string type;
};



#endif //MAIN_CPP_SST_H