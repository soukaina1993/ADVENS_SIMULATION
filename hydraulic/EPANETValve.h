//
// Created by Yolaine Adihou on 12/06/2021.
//

#ifndef FLUIDSMAIN_CPP_EPANETVALVE_H
#define FLUIDSMAIN_CPP_EPANETVALVE_H
#include "../EPANET/Elements/valve.h"
#include "../EPANET/Core/constants.h"
#include "HValve.h"
#include "../iofile/iofile.h"
#include "../state/TState.h"
#include "../fluids/PengRobinson.h"
#include <string>

class EPANETValve: public HValve {
public:
    EPANETValve(Valve::ValveType ValveType,string ValveName, double& Diameter, double&Flow);

    //EPANET

    std::string getType();
    void        setInitFlow();
    void        setVelocity();
    double      setRe(const double q, const double viscos);
    double      getOpenHeadLoss(double q); //  Find the head loss and its gradient for a fully open valve.
    double      getPbvHeadLoss(double q); //  Find the head loss and its gradient for a pressure breaker valve.
    double      getTcvHeadLoss(double q); //  Find the head loss and its gradient for a throttle control valve.
    double      getFcvHeadLoss(double q); //  Find the head loss and its gradient for a flow control valve.

    //ADVENS
    double      getHeadLoss( double q); //12.06.2021 Implémentation de la méthode générique pour le calcul des pertes de charge selon le type de vanne

    //BasicES
    void 	      init ();
    BasicES&  performance (); //sort tous les résulats
    float        objective(TSequence& seq);

private:
    Valve valve;

};


#endif //FLUIDSMAIN_CPP_EPANETVALVE_H
