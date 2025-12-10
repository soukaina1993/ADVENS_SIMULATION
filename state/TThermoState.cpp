//
// Created by Yolaine Adihou on 09/03/2021.
//

#include "TThermoState.h"
#include <iostream>

TThermoState::TThermoState(double &ACv, double &ACp, double &ASv, double &AAv, double &ABp, double &AGt) : TState(), Cv(ACv), Cp(ACp), Sv(ASv), Av(AAv), Bp(ABp), Gt(AGt){

}

TThermoState::TThermoState() : TState(){

}

void TThermoState::SetThermoState(TThermoState AthState) {
    v = AthState.v ;
    P = AthState.P ;
    T = AthState.T ;
    h = AthState.h ;
    s = AthState.s ;
    x = AthState.x ;
    Cv = AthState.Cv ;
    Cp = AthState.Cp ;
    Sv = AthState.Sv ;
    Av = AthState.Av ;
    Bp = AthState.Bp ;
    Gt = AthState.Gt ;
}

void TThermoState::GetThermoState(TThermoState AthState) {
    AthState.v = v ;
    AthState.P = P ;
    AthState.T = T ;
    AthState.h = h ;
    AthState.s = s ;
    AthState.x = x ;
    AthState.Cv = Cv ;
    AthState.Cp = Cp ;
    AthState.Sv = Sv ;
    AthState.Av = Av ;
    AthState.Bp = Bp ;
    AthState.Gt = Gt ;
}

void TThermoState::Display(){
    std::cout << "Cv [J/kg K] : " << Cv << std::endl;
    std::cout << "Cp [J/kg K] : " << Cp << std::endl;
}
