//
// Created by Yolaine Adihou on 07/06/2021.
//

#include "HValve.h"


void HValve::setFlow(double& q) {
    flow=q;
}

double HValve::getVelocity() {
    return velocity;
}

double HValve::getRe() {
    return Re;
}

double HValve::getFlow() {
    return flow;
}

double HValve::getDiameter() {
    return diameter;
}

void HValve::setDiameter(double & diam) {
    diameter=diam;

}

HValve::HValve(double& Diameter, double&Flow):diameter(Diameter),flow(Flow){

}



