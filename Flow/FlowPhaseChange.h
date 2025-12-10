//
// Created by lucile.schulthe on 09.05.2022.
//

#ifndef MAIN_CPP_FLUXPHASECHANGE_H
#define MAIN_CPP_FLUXPHASECHANGE_H


#include "Flow.h"
#include "../state/TThermoState.h"

class FlowPhaseChange: public Flow {        // TODO derive from Flow_physic?
public:
    explicit FlowPhaseChange(FluidCore* fluidin1, Type ityp=HOT) : Flow(fluidin1, ityp){
        StateIn=new TState;
        StateOut=new TState;
        Massflow=new double;
        States_created = Massflow_created = true;
    }

    void PSatCalculation(TState state);
    void PSatCalculation(double P);

    void TSatCalculation(TState state);
    void TSatCalculation(double T);

    void vCalculation(TPhysicState state);
    void vCalculation(TThermoState state);

    void PhysicCalculation(TPhysicState& pState);


    double Psat  (){return pSatState.P;}
    double Tsat  (){return pSatState.T;}
    double Psat  (const double Psat){
        tStategaz.P=Psat;
        tStateliquid.P=Psat;
        return pSatState.P=Psat;
    }
    double Tsat  (const double Tsat){
        tStategaz.T=Tsat;
        tStateliquid.T=Tsat;
        return pSatState.T=Tsat;
    }

    TState pSatState;
    TState tStateliquid, tStategaz;
};



#endif //MAIN_CPP_FLUXPHASECHANGE_H
