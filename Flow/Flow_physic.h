//
// Created by lucile.schulthe on 18.11.2022.
//

#ifndef FLUIDS_FLUX_PHYSIC_H
#define FLUIDS_FLUX_PHYSIC_H


#include "../state/TPhysicState.h"
#include "../fluids/FluidCore.h"
#include "Flow.h"
#include "FlowPhaseChange.h"

class Flow_physic: public Flow{
public:
    Flow_physic(bool withStateWall=false);
    Flow_physic(FluidCore* fluidin1, bool withStateWall=false);
    ~Flow_physic();
    void init(FluidCore* fluidin1) override;
    void copy(Flow_physic* in_flux);
    void PhysicCalculation(TPhysicState*& pState);  //from P,T --> not defined for wet steam
    void PhysicCalculation(TPhysicState& pState);
    void PhysicCalculation_vT(TPhysicState*& pState);   // from v,T
    void PhysicCalculation_vT(TPhysicState& pState);
    ////Set & Get


public:
    TPhysicState *pStateIn=nullptr, *pStateOut=nullptr;
    TPhysicState *pStateWall=nullptr;
};


#endif //FLUIDS_FLUX_PHYSIC_H
