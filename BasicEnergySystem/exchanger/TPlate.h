//
// Created by robin.lemaire on 11.08.2021.
//

#ifndef TEXCHANGER_H_TPLATE_H
#define TEXCHANGER_H_TPLATE_H

#include "TExchanger.h"

class TPlate : public TExchanger{
public:
    TPlate(TExchanger::Exchanger tp	, const TExchanger::Flowtype flw, const TExchanger::Type typ) :  TExchanger(tp, flw, typ) {};
    void init() override;
    void design() override;

    double N;   //number of plates

    ////parameters

    double s;   //size of a single plate
    double l;   //width of the plates
    double dPlate;  //gap in between 2 plates
    double C;   //corrugation angle
    double Rfc; //fouling resistance on cold side
    double Rfh; //fouling resistance on hot side
    double ep; //thickness of the plates
    double lambda;  //thermal conductivity of the plate

   // ~TPlate();

};

#endif //TEXCHANGER_H_TPLATE_H
