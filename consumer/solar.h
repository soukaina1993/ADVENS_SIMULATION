//
// Created by lucile.schulthe on 01.10.2021.
//

#ifndef MAIN_CPP_SOLAR_H
#define MAIN_CPP_SOLAR_H

#endif //MAIN_CPP_SOLAR_H


#include "../iofile/iofile.h"


class solar {
public:

    solar(){};
    vector <double> ElecinTime(double Pelec,environment &Tempex,double albedo=0.1);
    vector <double > calcul_AOI(double azimuth=0,double tilt=35/180*M_PI,double array_a=0);
    vector<double> calcul_coszen();
    double calcul_coszen(double day,double hour);
    vector<double> calcul_AOI(vector<double> coszen,double azimuth=0,double tilt=35/180*M_PI,double array_a=0);

    double thersol(double&Eg, double&dt, double &eta0,double &k1, double &k2);
    vector <double> thersolinTime(double m2, environment &Tempex,double T_cap);

    double latitude=23;



};