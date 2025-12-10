//
// Created by Yolaine Adihou on 12/06/2021.
//

#include "EPANETValve.h"
EPANETValve::EPANETValve(Valve::ValveType ValveType,string ValveName, double& Diameter, double&Flow): HValve(Diameter,Flow) , valve(ValveName) { //initialisation de la vanne avec le constructeur de Link de EPANET
    valve.valveType = ValveType ;
     valve.diameter=diameter;
    valve.flow=flow;
     valve.hLoss=headloss;
    Re= valve.getRe(flow, VISCOSITY); //viscosité de l'eau constante déclarée dans constants.h
    velocity= valve.getVelocity();
}

void EPANETValve::init() {

}

BasicES &EPANETValve::performance() {
    return *this;
}

float EPANETValve::objective(TSequence &seq) {
    return 0;
}

std::string EPANETValve::getType() {
    return valve.typeStr();
}

void EPANETValve::setInitFlow() {
    //valve.setInitFlow();
    valve.setInitFlow();
    flow=valve.flow;
}

void EPANETValve::setVelocity() {
   //valve.getVelocity();
   velocity=valve.getVelocity();
}

double EPANETValve::setRe(const double q, const double viscos) {
    Re= valve.getRe(q,viscos);
    return Re;
}

double EPANETValve::getOpenHeadLoss(double q) {
    valve.findOpenHeadLoss(q);
    headloss=valve.hLoss;
    return headloss;

}

double EPANETValve::getPbvHeadLoss(double q) {
    valve.findPbvHeadLoss(q);
    headloss=valve.hLoss;
    return headloss;
}

double EPANETValve::getTcvHeadLoss(double q) {
    valve.findTcvHeadLoss(q);
    headloss=valve.hLoss;
    return headloss;
}

double EPANETValve::getFcvHeadLoss(double q) {
    valve.findFcvHeadLoss(q);
    headloss=valve.hLoss;
    return headloss;
}

double EPANETValve::getHeadLoss(double q){

    switch (valve.valveType)
    {
        case Valve::PBV:
            getPbvHeadLoss(q);
            break;
        case Valve::TCV:
            getTcvHeadLoss(q);
            break;
        case Valve::FCV:
            getFcvHeadLoss(q);
            break;
}
}