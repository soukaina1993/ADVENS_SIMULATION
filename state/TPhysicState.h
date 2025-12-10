//
// Created by Yolaine Adihou on 09/03/2021.
//

#ifndef MAIN_CPP_TPHYSICSTATE_H
#define MAIN_CPP_TPHYSICSTATE_H
#include "TState.h"

class TPhysicState: public TState {
public:
    TPhysicState();
    TPhysicState(double& AViscDyn, double& AViscCin, double& ARo, double& AThermCond, double& APrandtl );
    void SetPhysicState (TPhysicState AphState);
    void GetPhysicState (TPhysicState AphState);
    void Display();
//private:
    double ViscDyn=0, ViscCin=0, Ro=0, ThermCond=0, Prandtl=0;


};



#endif //MAIN_CPP_TPHYSICSTATE_H
