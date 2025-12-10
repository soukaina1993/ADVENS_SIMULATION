//
// Created by lucile.schulthe on 01.10.2021.
//
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "../environment/environment.h"
#include "solar.h"
#include "../math/maths.h"


vector <double> solar::ElecinTime(double Pelec, environment &Tempex,double albedo){

    double a=1.040356;
    double b=-0.00007;
    double c=-0.000382;
    double POA,Irr_nor;
    vector <double> E ;

    vector<double> coszen=calcul_coszen();
    vector<double> cosAOI=calcul_AOI(coszen);

    for(int i=0;i<8760;i++){
    Irr_nor=(Tempex.Irr[i]-Tempex.Irr_diff[i])/coszen[i];
    POA=Irr_nor*cosAOI[i]+albedo*Tempex.Irr[i]+Tempex.Irr_diff[i];
    E.push_back(a*POA+b*(POA/Tempex.T[i])*(POA/Tempex.T[i])+c*POA*Tempex.T[i]);
    }
    return vN(E,Pelec);

}

vector<double> solar::calcul_AOI(double azimuth,double tilt,double array_a){
//tilt Array Tilt Angle
// array_a Array Azimuth Definition
    vector<double> cosAOI;
    double phi=latitude/180*M_PI;
double theta_day,theta_hour,coszen,zen;
    for(int day=1;day<366;day++) {
        for(int hour=1;hour<25;hour++) {
            theta_day=23.45/180*M_PI*sin(2*M_PI*(284+day)/365);
            theta_hour = M_PI * (12 - hour)/12;
            coszen=(sin(phi)*sin(theta_day)+cos(phi)*cos(theta_day)*cos(theta_hour));
            zen=acos(coszen);
            cosAOI.push_back(coszen*cos(tilt)+sin(zen)*cos(azimuth-array_a));
        }

    }
            return cosAOI;
}

vector<double> solar::calcul_AOI(vector<double> coszen,double azimuth,double tilt,double array_a){
//tilt Array Tilt Angle
// array_a Array Azimuth Definition

    vector<double> cosAOI;
    double zen;
    for(int i=0;i<8760;i++) {
            zen=acos(coszen[i]);
            cosAOI.push_back(coszen[i]*cos(tilt)+sin(zen)*cos(azimuth-array_a));
    }
    return cosAOI;
}

vector<double> solar::calcul_coszen(){

    vector<double> coszen;
    for(int day=1;day<366;day++) {
        for(int hour=1;hour<25;hour++) {
            coszen.push_back(calcul_coszen(day,hour));
        }

    }
    return coszen;
}

double solar::calcul_coszen(double day,double hour){
    double phi=latitude/180*M_PI;
    double theta_day,theta_hour;
    theta_day=23.45/180*M_PI*sin(2*M_PI*(284+day)/365);
    theta_hour = M_PI * (12 - hour)/12;
    return (sin(phi)*sin(theta_day)+cos(phi)*cos(theta_day)*cos(theta_hour));
}


double solar::thersol(double&Eg, double&dt, double &eta0,double &k1, double &k2){
     return eta0-((k1*dt/1000)-k2*dt*dt/1000);

}


vector <double> solar::thersolinTime(double m2, environment &Tempex,double T_cap){


    double Irr_nor;
    vector <double> E ;

    vector<double> coszen=calcul_coszen();

    double eta0=0.8;
    double k1=2.5;
    double k2=0.01;
    double dt;
    for(int i=0;i<8760;i++){
        if(coszen[i]>0){
            //Irr_nor=(Tempex.Irr[i]-Tempex.Irr_diff[i])/coszen[i];
            //cout<<"coszen"<<coszen[i]<<endl;
            //cout<<"Irr"<<Tempex.Irr[i]<<endl;
            dt=T_cap-Tempex.T[i];
            E.push_back(Tempex.Irr[i]*thersol(Tempex.Irr[i],dt,eta0,k1,k2));
        }
        else
            E.push_back(0);

    }

    return E;

}