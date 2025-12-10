//
// Created by robin.lemaire on 28.07.2021.
//

#ifndef TEXCHANGER_CPP_TTUBULAR_H
#define TEXCHANGER_CPP_TTUBULAR_H

#include "TExchanger.h"



class TTubular : public TExchanger{
public:
    enum fluidIn {HOT=0, COLD=1}; //definir le fluide à l'intérieur des tubes

    TTubular(TExchanger::Exchanger tp,  TExchanger::Flowtype flw ,  TExchanger::Type typ   ,  TTubular::fluidIn fld );

    fluidIn fluidI (const fluidIn fluid)	{return fld = fluid;}
    fluidIn fluidI () 			{return fld;  }

    void kCalculation() override;
    void convExchCoeffCalculation();

    double hE, hI;
    double DiamI, DiamH, DiamE;

    fluidIn fld;


    FlowPhaseChange HotChange;
    FlowPhaseChange ColdChange;

	void init() override;
};
#endif //TEXCHANGER_CPP_TTUBULAR_H
