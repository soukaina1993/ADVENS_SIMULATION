//
// Created by robin.lemaire on 30.07.2021.
//

#ifndef TTUBULAR_CPP_TDRYCOOLER_H
#define TTUBULAR_CPP_TDRYCOOLER_H

#include "TExchanger.h"

#include "../../Flow/FlowPhaseChange.h"

class TDryCooler : public TExchanger{
    TDryCooler(TExchanger::Exchanger tp	) : TExchanger(tp, TExchanger::COUNTER, TExchanger::COND) {};

    void kCalculation() override;
    void convExchCoeffCalculation();

    double hE, hI;
    double ng;              // rendement de la surface ailetée
    double Se, Si;          //surface extérieure et intérieure des tubes par mètres de longueur

    double DiamI, DiamE;    //Diametre Interieur et Diametre Exterieur
    double A;               // moitié de longueur de l'ailette
    double B;               // moitié de largeur de l'ailette
    double Rex;             // rayon trou ailette
    double eAil;            // epaisseur d'ailettes
    double n;               // nombre d'ailettes
    double Beta;            //coeff de dilatation

    void init() override;

    FlowPhaseChange HotChange;

};
#endif //TTUBULAR_CPP_TDRYCOOLER_H
