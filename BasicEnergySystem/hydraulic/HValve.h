//
// Created by lucile.schulthe on 21.11.2022.
//

#ifndef FLUIDS_HVALVE_H
#define FLUIDS_HVALVE_H

#include "../pipe/Hydraulic.h"

class HValve: public Hydraulic {

public:

    HValve(Flow *flux=new Flow){};
    HValve(double& Diameter, double&Massflow) ;
    void init (){};



    void Initflow(FluidCore*&iFluid, double iTin, double iPin,double iMFlow);

    void PressureDropCalculation();
    void HeatlossCalculation();
    void LossCalculation();

    void        setMassflow(double& q);
    double      getFlow();
    double      getVelocity();
    double      getRe();
    double      getDiameter();
    void        setDiameter(double & diam);
    double      setHeadLoss(double q);
    double      getHeadLoss();




    double Diam, headloss, Speed, Re,DeltaT, DeltaPr, Kv=100;


public:
    TPhysicState *pstate = new TPhysicState;
    Flow *flow= new Flow; // peut être une référence?
};

#endif //FLUIDS_HVALVE_H
