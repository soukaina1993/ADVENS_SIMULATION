//
// Created by lucile.schulthe on 28.04.2022.
//

#include "Flow.h"

Flow::Flow(FluidCore* fluidin1, TState *iStateIn, TState *iStateOut, Type ityp){
    typ=ityp;
    StateIn=iStateIn;
    StateOut=iStateOut;
    fluid=fluidin1;
    Massflow=new double;
    Massflow_created = true;
}

Flow::Flow() {
    typ=HOT;
    StateIn=new TState;
    StateOut=new TState;
    fluid=new MWater;
    Massflow=new double;
    States_created = Fluid_created = Massflow_created = true;
}

Flow::Flow(FluidCore* fluidin1, Type ityp) {
    typ=ityp;
    fluid=fluidin1;
    Massflow=new double;
    Massflow_created = true;
}

Flow::~Flow() {
    if (Massflow_created)
        delete Massflow;
    if (States_created) {
        delete StateIn;
        delete StateOut;
    }
    if (Fluid_created)
        delete fluid;
}

void Flow::init(){

}

void Flow::init(FluidCore* fluidin1) {
    fluid=fluidin1;
//    PhysicCalculation(StateIn, fluid);
//    PhysicCalculation(StateOut, fluid);
    StateCalculation(StateIn);
    StateCalculation(StateOut);
}
void Flow::exergyCalculation(double Ta){
    ex=(1-Ta/Tlog*(1-eta_s))*(*Massflow*abs(StateIn->h-StateOut->h));
    // TODO attention difference entre flux froid et chaud?
}

double Flow::getexergy(double Ta){
    return ex;
}

double Flow::getDeltah(){
    return StateOut->h-StateIn->h;
}

double Flow::getMassflow(){
    return *Massflow;
}

double Flow::getPower(){
    return Power;
}




void Flow::setStateIn_PT(double P,double T){
    StateIn->P=P;
    StateIn->T=T;
    StateIn->State_PT(fluid);
}

void Flow::setStateIn_Ps(double P,double s){
    StateIn->P=P;
    StateIn->s=s;
    StateIn->State_Ps(fluid);
}

void Flow::setStateOut_Ps(double P,double s){
    StateOut->P=P;
    StateOut->s=s;
    StateOut->State_Ps(fluid);
}

void Flow::setStateOut_PT(double P,double T){
    StateOut->P=P;
    StateOut->T=T;
    StateOut->State_PT(fluid);
}

void Flow::setStateOut_Ph(double P, double h){
    StateOut->State_Ph(fluid);
}

double Flow::getTlog(){
        return Tlog;
}


void Flow::etasCalculation(){
    double hs=fluid->h_Ps(StateOut->P,StateIn->s);

    if(typ==0)
        eta_s=(StateIn->h-StateOut->h)/(StateIn->h-hs);
    else
        eta_s=(StateIn->h-hs)/(StateIn->h-StateOut->h);
}

void Flow::connect(Flow* in_flux){
    StateIn=in_flux->StateOut;
    Massflow=in_flux->Massflow;
    fluid=in_flux->fluid;
}

void Flow::copy(Flow* in_flux){
    *StateIn=*in_flux->StateIn;
    *StateOut=*in_flux->StateOut;
    *Massflow=*in_flux->Massflow;
    *fluid=*in_flux->fluid;
    typ=in_flux->typ;
    SpecHeat=in_flux->SpecHeat;
    Tlog=in_flux->Tlog;
    eta_s=in_flux->eta_s;
    Cp=in_flux->Cp;
    Power=in_flux->Power;
    ex=in_flux->ex;
}


void Flow::setPout(double &P){
    StateOut->P=P;
    PressureCalculation();
}
void Flow::setPin(double &P)
{
    StateIn->P=P;
    PressureCalculation();
}
void Flow::setTin(double &T){
        StateIn->T=T;
}
void Flow::setTout(double &T){
    StateOut->T=T;
}

void Flow::PressureCalculation(){
    (StateIn->P+StateOut->P)/2;}

void Flow::TlogCalculation(){
    if (StateIn->T!=0 && StateOut->T!=0){
        if (StateIn->T==StateOut->T)
            Tlog=StateIn->T;
        else
            Tlog=(StateIn->T-StateOut->T)/log(StateIn->T/StateOut->T);
    }
    else
        Tlog=max(StateIn->T,StateOut->T);
}

double Flow::MassflowCalculation(){
    return Power/(StateIn->h-StateOut->h);
}

void Flow::StateCalculation(TState*& pState,FluidCore*fluidin1){
    fluid=fluidin1;
    StateCalculation(pState);
}
void Flow::StateCalculation(TState * &pState) {
    if (pState->T!=0){
        if (pState->v!=0){
            fluid->Prop_vT(pState->v,pState->T,pState->P,pState->u,pState->h,pState->s);
        }
        else if(pState->h!=0){
            double htemp;
            pState->v=fluid->v_Th_a(pState->T,pState->h);
            fluid->Prop_vT(pState->v,pState->T,pState->P,pState->u,htemp,pState->s);
        }
        else if(pState->P!=0){
            double ps=fluid->P_T_s(pState->T);
            if(pState->P>ps){                                   // cornelia.blanke: this is the liquid state
                //fluid->Prop_PT_g(pState->P,pState->T,pState->v,pState->u,pState->h,pState->s);
                fluid->Prop_PT_l(pState->P,pState->T,pState->v,pState->u,pState->h,pState->s);}
            else if (pState->P == ps){
                if (pState->x==-1) {
                    if (pState->h > 0){
                        pState->v=fluid->v_Th_a(pState->T,pState->h);
                        double vl,vg;
                        fluid->v_PT(pState->P,pState->T,vl,vg);
                        pState->x = Mfraction(vl, vg, pState->v);
                    }
                    else if (pState->s > 0){
                        pState->v=fluid->v_Ts_a(pState->T,pState->s);
                        double vl,vg;
                        fluid->v_PT(pState->P,pState->T,vl,vg);
                        pState->x = Mfraction(vl, vg, pState->v);
                    }
                    else {
                        pState->x = 0;
                        cerr << "Not enough data for biphasic state: x=0 is used" << endl;
                    }
                }
                fluid->Prop_Tx_s(pState->T,pState->x,pState->v,ps,pState->u,pState->h,pState->s);
            }
            else{
                //fluid->Prop_PT_l(pState->P,pState->T,pState->v,pState->u,pState->h,pState->s);
                fluid->Prop_PT_g(pState->P,pState->T,pState->v,pState->u,pState->h,pState->s);}
        }
        else if(pState->s!=0){
            double stemp;
            pState->v=fluid->v_Ts_a(pState->T,pState->s);
            fluid->Prop_vT(pState->v,pState->T,pState->P,pState->u,pState->h,stemp);
        }

        else if(pState->x!=-1){
            fluid->Prop_Tx_s(pState->T,pState->x,pState->v,pState->P,pState->u,pState->h,pState->s);
        }
        else {
            cerr << "Pas encore implementer"<<endl;
        }
    }
    else if(pState->P!=0){
        if (pState->v!=0){
            double p;
            pState->T=fluid->T_vP(pState->v,pState->P);
            fluid->Prop_vT(pState->v,pState->T,p,pState->u,pState->h,pState->s);
        }

        else if (pState->h!=0){
            double p,htemp;
            pState->v=fluid->v_Ph_a(pState->P,pState->h);
            pState->T=fluid->T_vP(pState->v,pState->P);
            fluid->Prop_vT(pState->v,pState->T,p,pState->u,htemp,pState->s);
        }
        else if(pState->s!=0){
            double p,stemp;
            pState->v=fluid->v_Ps_a(pState->P,pState->s);
            pState->T=fluid->T_vP(pState->v,pState->P);
            fluid->Prop_vT(pState->v,pState->T,p,pState->u,pState->h,stemp);
        }

        else if(pState->x!=-1){
            double p;
            pState->T=fluid->T_P_s(pState->P);
            fluid->Prop_Tx_s(pState->T,pState->x,pState->v,p,pState->u,pState->h,pState->s);
        }
        else {
            cout << "Pas encore implementer"<<endl;
        }
    }
    else {
        cerr << "need T or P to calculate" << endl;
        exit(1);
    }

}