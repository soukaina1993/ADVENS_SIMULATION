//
// Created by lucile.schulthe on 18.11.2022.
//

#include "Flow_physic.h"


// do not call default constructor of Flow because it would create new StateIn/StateOut !!
Flow_physic::Flow_physic(bool withStateWall) : Flow(new MWater) {
    pStateIn = new TPhysicState;
    pStateOut = new TPhysicState;
    if (withStateWall)
        pStateWall = new TPhysicState;
    StateIn = pStateIn;
    StateOut = pStateOut;
}
Flow_physic::Flow_physic(FluidCore* fluidin1, bool withStateWall) : Flow(fluidin1) {
    pStateIn = new TPhysicState;
    pStateOut = new TPhysicState;
    if (withStateWall)
        pStateWall = new TPhysicState;
    StateIn = pStateIn;
    StateOut = pStateOut;
}
Flow_physic::~Flow_physic() {
    delete pStateIn;
    delete pStateOut;
    delete pStateWall;
}

void Flow_physic::init(FluidCore* fluidin1) {
    fluid=fluidin1;
    StateCalculation(StateIn);
    StateCalculation(StateOut);
    PhysicCalculation(pStateIn);
    PhysicCalculation(pStateOut);
    if (pStateWall != nullptr){
        cerr << "Not yet implemented!" << endl;
        exit(1);
    }
}

void Flow_physic::copy(Flow_physic* in_flux){
    Flow::copy(in_flux);
    if (pStateIn != nullptr) {
        pStateIn->ViscDyn = in_flux->pStateIn->ViscDyn;
        pStateIn->ViscCin = in_flux->pStateIn->ViscCin;
        pStateIn->Ro = in_flux->pStateIn->Ro;
        pStateIn->ThermCond = in_flux->pStateIn->ThermCond;
        pStateIn->Prandtl = in_flux->pStateIn->Prandtl;
    }
    if (pStateOut != nullptr) {
        pStateOut->ViscDyn = in_flux->pStateOut->ViscDyn;
        pStateOut->ViscCin = in_flux->pStateOut->ViscCin;
        pStateOut->Ro = in_flux->pStateOut->Ro;
        pStateOut->ThermCond = in_flux->pStateOut->ThermCond;
        pStateOut->Prandtl = in_flux->pStateOut->Prandtl;
    }
    if (pStateWall != nullptr) {
        *pStateWall=*in_flux->pStateWall;
        pStateWall->ViscDyn=in_flux->pStateWall->ViscDyn;
        pStateWall->ViscCin=in_flux->pStateWall->ViscCin;
        pStateWall->Ro=in_flux->pStateWall->Ro;
        pStateWall->ThermCond=in_flux->pStateWall->ThermCond;
        pStateWall->Prandtl=in_flux->pStateWall->Prandtl;
    }
}

void Flow_physic::PhysicCalculation(TPhysicState & pstate ) {

    if(pstate.P!=0 && pstate.T!=0){
        double Psat=fluid->P_T_s(pstate.T);
        if (pstate.P < Psat)
            fluid->TransProp_PT_g(pstate.GetP(), pstate.GetT(), pstate.ViscDyn, pstate.ViscCin, pstate.Ro,
                                  pstate.ThermCond, pstate.Prandtl);
        else
            fluid->TransProp_PT_l(pstate.GetP(), pstate.GetT(),pstate.ViscDyn, pstate.ViscCin, pstate.Ro, pstate.ThermCond, pstate.Prandtl);
    }
    else{
        cerr << "P==0 or T==0" << endl;
        exit(1);
    }
}

void Flow_physic::PhysicCalculation(TPhysicState*& pstate ) {

    if(pstate->P!=0 && pstate->T!=0){
        double Psat=fluid->P_T_s(pstate->T);
        if (pstate->P < Psat)
            fluid->TransProp_PT_g(pstate->GetP(), pstate->GetT(), pstate->ViscDyn, pstate->ViscCin, pstate->Ro, pstate->ThermCond,pstate->Prandtl);
        else
            fluid->TransProp_PT_l(pstate->GetP(), pstate->GetT(),pstate->ViscDyn, pstate->ViscCin, pstate->Ro, pstate->ThermCond, pstate->Prandtl);
    }
    else{
        cerr << "P==0 or T==0" << endl;
        exit(1);
    }
}

void Flow_physic::PhysicCalculation_vT(TPhysicState*& pstate ) {

    if(pstate->v!=0 && pstate->T!=0){
        double Ts = fluid->T_P_s(pstate->P);
        double vl,vg;
        fluid->v_PT(pstate->P,Ts,vl,vg);
        if (pstate->v > vg){
            fluid->TransProp_vT_g(pstate->GetV(), pstate->GetT(), pstate->ViscDyn, pstate->ViscCin, pstate->Ro, pstate->ThermCond, pstate->Prandtl);
        }
        else if (pstate->v < vl){
            fluid->TransProp_vT_l(pstate->GetV(), pstate->GetT(), pstate->ViscDyn, pstate->ViscCin, pstate->Ro, pstate->ThermCond, pstate->Prandtl);
        }
        else{
            TPhysicState pstate_l, pstate_g;
            fluid->TransProp_vT_g(vg, pstate->GetT(), pstate_g.ViscDyn, pstate_g.ViscCin, pstate_g.Ro, pstate_g.ThermCond, pstate_g.Prandtl);
            fluid->TransProp_vT_l(vl, pstate->GetT(), pstate_l.ViscDyn, pstate_l.ViscCin, pstate_l.Ro, pstate_l.ThermCond, pstate_l.Prandtl);
            if (pstate->x == -1)    // not yet calculated
                pstate->x = Mfraction(vl,vg,pstate->GetV());
            pstate->ViscDyn = (1 - pstate->x) * pstate_l.ViscDyn + pstate->x * pstate_g.ViscDyn;
            pstate->ViscCin = (1 - pstate->x) * pstate_l.ViscCin + pstate->x * pstate_g.ViscCin;
            pstate->Ro = 1.0 / pstate->GetV();
            pstate->ThermCond = (1 - pstate->x) * pstate_l.ThermCond + pstate->x * pstate_g.ThermCond;
            pstate->Prandtl = (1 - pstate->x) * pstate_l.Prandtl + pstate->x * pstate_g.Prandtl;
        }
    }
    else{
        cerr << "v==0 or T==0" << endl;
        exit(1);
    }
}

void Flow_physic::PhysicCalculation_vT(TPhysicState& pstate ) {

    if(pstate.v!=0 && pstate.T!=0){
        double Ts = fluid->T_P_s(pstate.P);
        double vl,vg;
        fluid->v_PT(pstate.P,Ts,vl,vg);
        if (pstate.v > vg){
            fluid->TransProp_vT_g(pstate.GetV(), pstate.GetT(), pstate.ViscDyn, pstate.ViscCin, pstate.Ro, pstate.ThermCond, pstate.Prandtl);
        }
        else if (pstate.v < vl){
            fluid->TransProp_vT_l(pstate.GetV(), pstate.GetT(), pstate.ViscDyn, pstate.ViscCin, pstate.Ro, pstate.ThermCond, pstate.Prandtl);
        }
        else {
            TPhysicState pstate_l, pstate_g;
            fluid->TransProp_vT_g(vg, pstate.GetT(), pstate_g.ViscDyn, pstate_g.ViscCin, pstate_g.Ro,
                                  pstate_g.ThermCond, pstate_g.Prandtl);
            fluid->TransProp_vT_l(vl, pstate.GetT(), pstate_l.ViscDyn, pstate_l.ViscCin, pstate_l.Ro,
                                  pstate_l.ThermCond, pstate_l.Prandtl);
            if (pstate.x == -1)    // not yet calculated
                pstate.x = Mfraction(vl, vg, pstate.GetV());
            pstate.ViscDyn = (1 - pstate.x) * pstate_l.ViscDyn + pstate.x * pstate_g.ViscDyn;
            pstate.ViscCin = (1 - pstate.x) * pstate_l.ViscCin + pstate.x * pstate_g.ViscCin;
            pstate.Ro = 1.0 / pstate.GetV();
            pstate.ThermCond = (1 - pstate.x) * pstate_l.ThermCond + pstate.x * pstate_g.ThermCond;
            pstate.Prandtl = (1 - pstate.x) * pstate_l.Prandtl + pstate.x * pstate_g.Prandtl;
        }
    }
    else{
        cerr << "v==0 or T==0" << endl;
        exit(1);
    }
}