//
// Created by Yolaine Adihou on 09/03/2021.
//

#include "TPhysicState.h"
#include <iostream>

TPhysicState::TPhysicState() : TState(){

}

TPhysicState::TPhysicState(double& AViscDyn, double& AViscCin, double& ARo, double& AThermCond, double& APrandtl) : ViscDyn(AViscDyn), ViscCin(AViscCin), Ro(ARo), ThermCond(AThermCond), Prandtl(APrandtl) {

}

void TPhysicState::SetPhysicState(TPhysicState AphState) {
    v = AphState.v ;
    P = AphState.P ;
    T = AphState.T ;
    h = AphState.h ;
    s = AphState.s ;
    x = AphState.x ;
    ViscDyn = AphState.ViscDyn ;
    ViscCin = AphState.ViscCin ;
    Ro = AphState.Ro ;
    ThermCond = AphState.ThermCond ;
    Prandtl = AphState.Prandtl ;
}

void TPhysicState::GetPhysicState(TPhysicState AphState) {
    AphState.v = v ;
    AphState.P = P ;
    AphState.T = T ;
    AphState.h = h ;
    AphState.s = s ;
    AphState.x = x ;
    AphState.ViscDyn = ViscDyn ;
    AphState.ViscCin = ViscCin ;
    AphState.Ro = Ro ;
    AphState.ThermCond = ThermCond ;
    AphState.Prandtl = Prandtl ;
}

void TPhysicState::Display() {
   // std::cout << "cv [J/kg K] : " << cv << std::endl;
   // std::cout << "cp [J/kg K] : " << cp << std::endl;
   TState::Display();
    std::cout << "Ro    [kg/m^3] : " << Ro << std::endl;
    std::cout << "ViscDyn   [kg/m s] : " << ViscDyn << std::endl;
    std::cout << "ViscCin   [m^2/s] : " << ViscCin << std::endl;
    std::cout << "ThermCond [W/m K] : " << ThermCond << std::endl ;
    std::cout << "Prandtl   [-] : " << Prandtl << std::endl;
}


