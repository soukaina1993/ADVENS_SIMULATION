//
// Created by lucil on 17.11.2021.
//
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "../iofile/iofile.h"
#include "Ener_env.h"
#include "../environment/environment.h"
#include "../fluids/FluidCore.h"

using namespace std;


vector<double> Ener_env::E_Tfix(vector<double> &flow,double &T,double &Tmin,double &P,FluidCore &fluid){
    double cv,cp,cpmin,a;
    vector<double> E(flow.size());

        double Rho=998;//1/fluid.v_P T_l(P,Tmin);
        fluid.CvCp_vT_l(fluid.v_PT_l(P,T),T,cv,cp);
        fluid.CvCp_vT_l(1/Rho,Tmin,cv,cpmin);
        a=(cp+cpmin)/2 *(T - Tmin)*Rho/3.6e6;


    for (int i(0); i < flow.size(); i++){
        E[i] =  a* flow[i];}
    return E;
}

vector<double> Ener_env::E_Tvar(vector<double> &flow,vector<double> &T,double &Tmin,double &P,FluidCore &fluid){
    double cv,cp,cpmin,Rho;
    vector<double> E(flow.size());


        Rho=998;//Rho=1/fluid.v_PT_l(P,Tmin);
        fluid.CvCp_vT_l(1/Rho,Tmin,cv,cpmin);


    for (int i(0); i < flow.size(); i++){

            fluid.CvCp_vT_l(1/Rho,T[i],cv,cp);
        E[i] =  (cp+cpmin)/2 *(T[i] - Tmin)*Rho/3.6e6* flow[i];}

    return E;
}

vector<double> Ener_env::Tmix(vector<double>& flow1,vector<double>& T1,vector<double> &flow2,double &T2){
    vector <double> flow(flow1.size());
    vector <double>T(flow.size());


    for (int i(0); i < flow.size(); i++){
        flow[i]=flow1[i]+flow2[i];
        if (flow[i]==0)
            T[i]=T2;
        else
            T[i]=(flow1[i] * T1[i] + flow2[i] * T2) / flow[i];}

    return T;
}

vector<double> Ener_env::E_Tmix(vector<double> &flow1,vector<double> &T1,vector<double> &flow2,double &T2,double &Tmin,double &P,FluidCore &fluid){
    vector <double> flow(flow1.size());
    vector <double>T(flow.size());

    T=Tmix(flow1,T1,flow2,T2);
    for (int i(0); i < flow.size(); i++)
        flow[i]=flow1[i]+flow2[i];

    double cv,cp,cpmin,Rho;
    vector<double> E(flow.size());


        Rho=998; //Rho=1/fluid.v_PT_l(P,Tmin);
        fluid.CvCp_vT_l(1/Rho,Tmin,cv,cpmin);


    for (int i(0); i < flow.size(); i++){

            fluid.CvCp_vT_l(1/Rho,T[i],cv,cp);
        E[i] =  (cp+cpmin)/2 *(T[i] - Tmin)*Rho/3.6e6* flow[i];}

    return E;
}