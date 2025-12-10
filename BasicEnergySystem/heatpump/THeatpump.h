//
// Created by lucile.schulthe on 31.03.2022.
//

#ifndef MAIN_CPP_THEATPUMP_H
#define MAIN_CPP_THEATPUMP_H


#include "../BasicES.h"
//#include "../../fluids/FluidCore.h"
#include "../../Flow/FlowPhaseChange.h"



class THeatpump {
public:
    enum Type {DEFAULT=0, AIR=1, WATER=2, GEO=3};
    enum Heattype {HEAT=0, COOL=1, BOTH=2};
    enum Criteria  {RENDEXERG=0, POLL=1, COST=2};  //3 aspects à optimiser : rendement exergétique, pollution et coût --> peut être mis sur Basic energy system

    // Constructor, Destructor
    explicit THeatpump( FluidCore* &fluidin1, double ieta_comp=0.8, double ipinchCond=5, double ipinchEvap=5,
                        THeatpump::Type typ=WATER, THeatpump::Heattype heat=HEAT);
    explicit THeatpump( string fluidin1="R134a", double ieta_comp=0.8, double ipinchCond=5, double ipinchEvap=5,
                        THeatpump::Type typ=WATER, THeatpump::Heattype heat=HEAT);
    ~THeatpump();

    // herited
    void init ();
    BasicES & performance();
    float objective(TSequence& seq);
    void display();


    void COP_real(double &Thot,double &Tcold,double& Text);
    void COP_real(double &Thot,double &Tcold);
    double COP_max(const double &Thot, const double& Tcold);     // COP Carnot

    void compute_cycle(double &Tcond, double &Tevap);
    double DeltaS_Cp(double &Tcond, double &Tevap);

    void getexergy_old(double &Thot, double& Tcold, double& Text);
    void getexergy_old(double Text=273.15);
    //void getexergy(double Text=273.15);       // not implemented

    //void getMassFlow();
    double getCost();                         //TODO not implemented


    Heattype heattype (const Heattype HType)	{return heat = HType;}
    Heattype heattype () 			{return heat;  }

    Criteria crit(const Criteria Crit) {return obj=Crit;}
    Criteria crit ()            {return obj;  }

    Type type(const Type AType) {return typ=AType;}
    Type type ()             {return typ;   }

public:
    Criteria	obj;				//optimization criteria
    Type        typ;                //Default,Air,Water,Geo
    Heattype    heat;               //Heat,Cool,Both

    double eta_comp=0.8, pinchEvap=5, pinchCond=5;
    double T2_log, Cost=0;
    double eta_ex=-1;
    double COP=0;
    //double Massflow_in;

    FluidCore* fluidin= nullptr; // fluid inside Heat Pump (e.g. R134a...)

    // neighbouring fluids (e.g. Air, Water...)
    Flow *Hot= nullptr, *Cold= nullptr;       //Flow pointer to Flow Hot=Building, Cold=Network

private:
    // includes fluid pointer to fluid inside Heat Pump (e.g. R134a...)
    FlowPhaseChange *cond= nullptr, *evap= nullptr;
    bool Fluid_created=false;
};

#endif //MAIN_CPP_THEATPUMP_H
