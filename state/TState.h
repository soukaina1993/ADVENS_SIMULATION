//
// Created by Yolaine Adihou on 09/03/2021.
//

#ifndef MAIN_CPP_TSTATE_H
#define MAIN_CPP_TSTATE_H

#include "../fluids/FluidCore.h"

class TState {
public:
    TState() = default;
    TState(double iAv, double iAP, double iAT, double iAh, double iAs, double iAx);
    //void init();

    void StateCalculation(FluidCore* fluid);
    void State_PT(FluidCore* fluid);
    void State_Ps(FluidCore* fluid);
    void State_Ph(FluidCore* fluid);

    void SetV (const double iAv) { v = iAv; }
    void SetP (const double iAP) { P = iAP; }
    void SetT (const double iAT) { T = iAT; }
    void SetH (const double iAh) { h = iAh; }
    void SetS (const double iAs) { s = iAs; }
    void SetX (const double iAx) { x = iAx; }

    void SetState (const TState& AState);
    void SetState (TState *AState);

    double GetV() { return v; }
    double GetP() {return P; }
    double GetT() { return T; }
    double GetH() { return h; }
    double GetS() { return s; }
    double GetX() { return x; }

    void GetState (TState& AState); //pour fixer les valeurs du State en paramètre
    void GetState (TState* AState);
    void GetState(const char* file_in, int index=1);

    void Display();

    // pour module fluide
    void assign(string &var1, double &arg1);
    void calcul_state(FluidCore* &fluid);
    void init(FluidCore* &fluid);
    void init_all(FluidCore* &fluid);
    void init_all_x(FluidCore* &fluid);     // includes calculation of x
    void Display(string &var1);

public:
    double v=0, P=0, T=0, h=0, s=0, x=-1,u=0;
};


#endif //MAIN_CPP_TSTATE_H
