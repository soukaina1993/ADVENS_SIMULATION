//
// Created by Yolaine Adihou on 09/03/2021.
//

#ifndef MAIN_CPP_TTHERMOSTATE_H
#define MAIN_CPP_TTHERMOSTATE_H
#include "TState.h"

class TThermoState: public TState {
public:
    TThermoState();
    TThermoState(double& ACv, double& ACp, double& ASv, double& AAv, double& ABp, double& AGt);
    void SetThermoState (TThermoState AthState);
    void GetThermoState (TThermoState AthState);
    void Display();
//private:
    double Cv, Cp, Sv, Av, Bp, Gt ;
};


#endif //MAIN_CPP_TTHERMOSTATE_H
