//
// Created by lucile.schulthe on 21.11.2022.
//

#ifndef FLUIDS_HPUMP_H
#define FLUIDS_HPUMP_H
#include "../pipe/Hydraulic.h"
class HPump: public Hydraulic {

public:

    HPump(Flow *flux=new Flow);
    void init ();

    void PowerCalculation();


    void display();


    //-- GET SET ---//

    void        setFlow(double& q);
    double      getFlow();
    double      getVelocity();
    double      getRe();
    double      getDiameter();
    void        setDiameter(double & diam);




protected:
    double Diam, headloss, Speed, Re,DeltaPr;
public:
    TPhysicState *pstate = new TPhysicState;
// add for network from Matlab
    Flow *flow = new Flow; // peut être une référence?
};
#endif //FLUIDS_HPUMP_H
