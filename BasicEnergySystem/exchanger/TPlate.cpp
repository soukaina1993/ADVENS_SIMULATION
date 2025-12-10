//
// Created by robin.lemaire on 11.08.2021.
//

#include "TPlate.h"

/****************************************************************************/

void TPlate::design() {
    getPower();   //calculate the required heat flux
    deltaTlm();   //logarithmic mean temp difference
    double kNew;

    kNew = 4500; //usually, the global transfer coefficient for plate heat exchanger is between 2-7 kW.m².K-1
    //this value is a first guess

    while(abs(k-kNew)>0.1)
    {
        k=kNew;

        Area = Q / k / DeltaTlog;

        N = Area / s;

        //Now we want to verify this rough value of the heat exchanger size

        double Dh = (4 * l * dPlate) / (2 * (1 + dPlate));  //hydraulic diameter

        double coldu = *Cold.Massflow * (Cold.pStateIn->GetV() + Cold.pStateOut->GetV()) / 2 / (3.14 * Dh * Dh);
        double hotu = *Hot.Massflow * ((Hot.pStateIn->GetV() + Hot.pStateOut->GetV()) / 2) / (3.14 * Dh * Dh);

        double coldRe = coldu * Dh / Cold.pStateIn->ViscCin;
        double hotRe = hotu * Dh / ((Hot.pStateIn->ViscCin + Hot.pStateOut->ViscCin) / 2);

        double coldA=1;
        double hotA=1;
        double coldB=1;
        double hotB=1;

        if (C < 30 && coldRe < 10) {
            coldA = 0.718;
            coldB = 0.349;
        }
        else if (C < 30 && coldRe >= 10) {
            coldA = 0.348;
            coldB = 0.663;
        }
        else if (C == 45 && coldRe < 10) {
            coldA = 0.718;
            coldB = 0.349;
        }
        else if (C == 45 && 10 <= coldRe && coldRe < 100) {
            coldA = 0.4;
            coldB = 0.598;
        }
        else if (C == 45 && coldRe >= 100) {
            coldA = 0.3;
            coldB = 0.663;
        }
        else if (C == 50 && coldRe < 20) {
            coldA = 0.63;
            coldB = 0.333;
        }
        else if (C == 50 && coldRe >= 20 && coldRe < 300) {
            coldA = 0.291;
            coldB = 0.591;
        }
        else if (C == 50 && coldRe >= 300) {
            coldA = 0.13;
            coldB = 0.732;
        }
        else if (C == 60 && coldRe < 20) {
            coldA = 0.562;
            coldB = 0.326;
        }
        else if (C == 60 && coldRe >= 20 && coldRe < 400) {
            coldA = 0.306;
            coldB = 0.529;
        }
        else if (C == 60 && coldRe >= 400) {
            coldA = 0.108;
            coldB = 0.703;
        }
        else if (C > 65 && coldRe < 20) {
            coldA = 0.562;
            coldB = 0.326;
        }
        else if (C > 65 && coldRe >= 20 && coldRe < 500) {
            coldA = 0.331;
            coldB = 0.503;
        }
        else if (C > 65 && coldRe > 500) {
            coldA = 0.087;
            coldB = 0.718;
        }
        else
            cerr << "coldA, coldB cannot be computed" << endl;

        if (C < 30 && hotRe < 10) {
            hotA = 0.718;
            hotB = 0.349;
        }
        else if (C < 30 && hotRe >= 10) {
            hotA = 0.348;
            hotB = 0.663;
        }
        else if (C == 45 && hotRe < 10) {
            hotA = 0.718;
            hotB = 0.349;
        }
        else if (C == 45 && 10 <= hotRe && hotRe < 100) {
            hotA = 0.4;
            hotB = 0.598;
        }
        else if (C == 45 && hotRe >= 100) {
            hotA = 0.3;
            hotB = 0.663;
        }
        else if (C == 50 && hotRe < 20) {
            hotA = 0.63;
            hotB = 0.333;
        }
        else if (C == 50 && hotRe >= 20 && hotRe < 300) {
            hotA = 0.291;
            hotB = 0.591;
        }
        else if (C == 50 && hotRe >= 300) {
            hotA = 0.13;
            hotB = 0.732;
        }
        else if (C == 60 && hotRe < 20) {
            hotA = 0.562;
            hotB = 0.326;
        }
        else if (C == 60 && hotRe >= 20 && hotRe < 400) {
            hotA = 0.306;
            hotB = 0.529;
        }
        else if (C == 60 && hotRe >= 400) {
            hotA = 0.108;
            hotB = 0.703;
        }
        else if (C > 65 && hotRe < 20) {
            hotA = 0.562;
            hotB = 0.326;
        }
        else if (C > 65 && hotRe >= 20 && hotRe < 500) {
            hotA = 0.331;
            hotB = 0.503;
        }
        else if (C > 65 && hotRe >= 500) {
            hotA = 0.087;
            hotB = 0.718;
        }
        else
            cerr << "hotA, hotB cannot be computed" << endl;

        double coldNu = coldA * pow(coldRe, coldB) * pow((Cold.pStateIn->Prandtl + Cold.pStateOut->Prandtl) / 2, 0.33) *
                 pow((Cold.pStateIn->Prandtl + Cold.pStateOut->Prandtl) / 2 / Cold.pStateWall->Prandtl, 0.13);
        double hotNu = hotA * pow(hotRe, hotB) * pow((Hot.pStateIn->Prandtl + Hot.pStateOut->Prandtl) / 2, 0.33) *
                pow((Hot.pStateIn->Prandtl + Hot.pStateOut->Prandtl) / 2 / Hot.pStateWall->Prandtl, 0.13);

        double hC = coldNu * ((Cold.pStateIn->ViscDyn + Cold.pStateOut->ViscDyn) / 2) / Dh;
        double hH = hotNu * ((Hot.pStateIn->ViscDyn + Hot.pStateOut->ViscDyn) / 2) / Dh;

        kNew = 1 / (1 / hC + Rfc + 1 / hH + Rfh + ep / lambda);
    }
    k=kNew;
    Area = Q / k / DeltaTlog;
    N = Area / s;
    cout <<"k= "<<k<<endl;
}

/****************************************************************************/

void TPlate::init()
{
    TExchanger::init();
    s=a;
    l=b;
    dPlate=c;
    C=d;
    Rfc=e;
    Rfh=f;
    ep=g;
    lambda=h;
}