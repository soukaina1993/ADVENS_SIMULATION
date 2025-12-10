//
// Created by Yolaine Adihou on 09/03/2021.
//

#ifndef MAIN_CPP_TSONICSTATE_H
#define MAIN_CPP_TSONICSTATE_H

#include "TState.h"

class TSonicState: public TState {
public:
    TSonicState();
    TSonicState(double& AvLaval, double& APLaval, double& ATLaval, double& ACsonic);
    void SetSonicState (TSonicState AsoState);
    void GetSonicState (TSonicState AsoState);
    void SetVL (double ivl);
    void SetPL (double iPl);
    void SetTL (double iTl);
    void SetCs (double iCs);

    double GetVL();
    double GetPL();
    double GetTL();
    double GetCs();

private:
    double vLaval, PLaval, TLaval, Csonic ;
};


#endif //MAIN_CPP_TSONICSTATE_H
