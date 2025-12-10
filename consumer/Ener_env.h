//
// Created by lucil on 17.11.2021.
//

#ifndef MAIN_CPP_ENER_ENV_H
#define MAIN_CPP_ENER_ENV_H
#include "../environment/environment.h"
#include "../fluids/FluidCore.h"
#endif //MAIN_CPP_ENER_ENV_H
class Ener_env {
public:
    vector<double> E_Tfix(vector<double> &flow,double &T,double &Tmin,double &P,FluidCore &fluid);

    vector<double> E_Tvar(vector<double> &flow,vector<double> &T,double &Tmin,double &P,FluidCore &fluid);

    vector<double> Tmix(vector<double> &flow1,vector<double> &T1,vector<double> &flow2,double &T2);

    vector<double> E_Tmix(vector<double> &flow1,vector<double> &T1,vector<double> &flow2,double &T2,double &Tmin,double &P,FluidCore &fluid);

};