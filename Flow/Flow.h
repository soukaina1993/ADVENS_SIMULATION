//
// Created by lucile.schulthe on 28.04.2022.
//

#ifndef MAIN_CPP_FLUX_H
#define MAIN_CPP_FLUX_H


#include "../state/TPhysicState.h"
#include "../fluids/FluidCore.h"
#include "../fluids/LeeKessler.h"
#include "../fluids/MWater.h"

class Flow {
public:
    enum Type{HOT=0, COLD=1};  // Flow chaud diminue d'enthalpie- Flow froid augmente son enthalpie


    Flow(Type ityp) {typ=ityp;};
    Flow(TState *&ipStateIn, TState *&ipStateOut, Type ityp=HOT) {StateIn=ipStateIn; StateOut=ipStateOut;};
    //Flow(FluidCore* fluidin1=new MWater,TState *iStateIn=new TState,TState *iStateOut=new TState, Type ityp=HOT);

    Flow();
    Flow(FluidCore* fluidin1, Type ityp=HOT);
    Flow(FluidCore* fluidin1, TState *iStateIn, TState *iStateOut, Type ityp=HOT);
    ~Flow();    //everything created with "new" needs to be deleted

    void init();
    virtual void init(FluidCore* fluidin1);
    void connect(Flow* in_flux);
    void copy(Flow* in_flux);


    void etasCalculation();
    void exergyCalculation(double Ta);
    void PressureCalculation();
    void TlogCalculation();
    void StateCalculation(TState* &pState);
    void StateCalculation(TState*& pState,FluidCore*fluid);
    double MassflowCalculation();

    ////Set & Get
    void setStateIn_PT(double P,double T);
    void setStateIn_Ps(double P,double s);
    void setStateOut_Ps(double P,double s);
    void setStateOut_PT(double P, double T);
    void setStateOut_Ph(double P, double h);
    double getMassflow(); // calcul la différence d'enthlapie


    double getTlog();
    double getPower(); // calcul la puissance
    double getexergy(double Ta); // calcul l'exergy a partir d'un temperature de ref
    double getDeltah(); // calcul la différence d'enthlapie

    double massflow 	    ()		{return *Massflow;}
    double massflow 	    (double &M)		{return *Massflow=M;}
    double TIn    ()		    {return StateIn->T;}
    double TIn (double Ti)		{return StateIn->T=Ti;}

    double T    ()		    {return StateIn->T;}
    double T (double Ti)		{return StateIn->T=Ti;}
    double TOut     ()		{return StateOut->T;}
    double TOut  (double To)		{return StateOut->T=To;}
    double xIn      ()          {return StateIn->x;}
    double xIn      (double xI)     {return StateIn->x=xI;}
    double xOut      ()          {return StateOut->x;}
    double xOut      (double xO)     {return StateOut->x=xO;}
    double vIn      ()          {return StateIn->v;}
    double vIn      (double vI)     {return StateIn->v=vI;}
    double vOut      ()          {return StateOut->v;}
    double vOut      (double vO)     {return StateOut->v=vO;}
    double PIn      ()          {return StateIn->P;}
    double PIn      (double Pi)     {return StateIn->P=Pi;}
    double POut      ()          {return StateOut->P;}
    double POut      (double PO)     {return StateOut->P=PO;}
    double specHeat       ()		{return SpecHeat;}
    double specHeat  (double iCp)		{return SpecHeat=iCp;}



    double getT(){return (StateIn->T+StateOut->T)/2;};
    double getP(){return (StateIn->P+StateOut->P)/2;};

    void setPout(double &P);
    void setPin(double &P);
    void setTin(double &T);
    void setTout(double &T);


public:
    FluidCore* fluid;
    double 	SpecHeat,  Tlog, eta_s=1, Cp=0, Power, ex;    // SpecHeat == Cp ?!
    double *Massflow;
    TState *StateIn, *StateOut;

    Type        typ;                        //Type

protected:
    bool Massflow_created=false, States_created=false, Fluid_created=false;
};



#endif //MAIN_CPP_FLUX_H
