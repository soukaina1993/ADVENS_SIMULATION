//
// Created by lucile.schulthe on 09.05.2022.
//

#include "FlowPhaseChange.h"



void FlowPhaseChange::PSatCalculation(TState state)
{
    PSatCalculation(state.P);
}

void FlowPhaseChange::PSatCalculation(const double P)
{
    double Tsat;
    fluid->Prop_P_s(P, tStateliquid.v, tStategaz.v, Tsat, tStateliquid.u,  tStategaz.u, tStateliquid.h,  tStategaz.h, tStateliquid.s,  tStategaz.s);
    tStateliquid.P=P;
    tStateliquid.T=Tsat;
    tStategaz.P=P;
    tStategaz.T=Tsat;
}

void FlowPhaseChange::TSatCalculation(TState state)
{
    TSatCalculation(state.T);
}

void FlowPhaseChange::TSatCalculation(const double T)
{
    double Psat;
    fluid->Prop_T_s(T, tStateliquid.v, tStategaz.v, Psat, tStateliquid.u,  tStategaz.u, tStateliquid.h,  tStategaz.h, tStateliquid.s,  tStategaz.s);
    tStateliquid.P=Psat;
    tStateliquid.T=T;
    tStategaz.P=Psat;
    tStategaz.T=T;
}

void FlowPhaseChange::vCalculation(TPhysicState state)
{
    fluid->v_PT(state.P,state.T, tStateliquid.v, tStategaz.v);
}

void FlowPhaseChange::vCalculation(TThermoState state)
{
    fluid->v_PT(state.P,state.T, tStateliquid.v, tStategaz.v);
}

void FlowPhaseChange::PhysicCalculation(TPhysicState& pState) {
    cerr<<"pas encore implementer"<<endl;
}

