//
// Created by lucile.schulthe on 21.11.2022.
//



#include "HValve.h"
void HValve::Initflow(FluidCore*&iFluid, double iTin, double iPin,double iMFlow){
    flow->fluid=iFluid;
    flow->StateIn->T=iTin;
    flow->StateIn->P=iPin;
    *flow->Massflow=iMFlow;
}

void HValve::LossCalculation(){
    PressureDropCalculation();
    double Pout=flow->StateIn->P-DeltaPr;
    flow->setPout(Pout);

    HeatlossCalculation();
    double Tout=flow->StateIn->T-DeltaT;
    flow->setTout(Tout);
}

void HValve::HeatlossCalculation(){
    DeltaT=1; //TODO a calculer


}
void HValve::PressureDropCalculation(){
   // headloss=100; //TODO a calculer
    DeltaPr=sqr(*flow->Massflow)/sqr(Kv);

}

double  HValve::setHeadLoss (double q){
    *flow->Massflow=q;
    PressureDropCalculation();
    double Pout=flow->StateIn->P-DeltaPr;
    flow->setPout(Pout);
    return headloss;
}

double  HValve::getHeadLoss (){
    return headloss;
}

void HValve::setMassflow(double& q) {
    *flow->Massflow=q;
}

double HValve::getVelocity() {
    return Speed;
}

double HValve::getRe() {
    return Re;
}

double HValve::getFlow() {
    return *flow->Massflow;
}

double HValve::getDiameter() {
    return Diam;
}

void HValve::setDiameter(double & diam) {
    Diam=diam;

}

HValve::HValve(double& Diameter, double&Massflow):Diam(Diameter){
        *flow->Massflow=Massflow;
}



