//
// Created by Yolaine Adihou on 07/06/2021.
//

#ifndef FLUIDSMAIN_CPP_HVALVE_H
#define FLUIDSMAIN_CPP_HVALVE_H
#include "Hydraulic.h"
#include <string>

class HValve : public Hydraulic {

public:

    HValve(double& Diameter, double&Flow) ;
    void        setFlow(double& q);
    double      getFlow();
    double      getVelocity();
    double      getRe();
    double      getDiameter();
    void        setDiameter(double & diam);
    virtual double      getHeadLoss(double q)=0;



protected:
    double diameter, flow, headloss, velocity, Re;
};


#endif //FLUIDSMAIN_CPP_HVALVE_H
