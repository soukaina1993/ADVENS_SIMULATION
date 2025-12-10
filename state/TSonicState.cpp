//
// Created by Yolaine Adihou on 09/03/2021.
//

#include "TSonicState.h"

TSonicState::TSonicState() : TState(){

}

TSonicState::TSonicState(double &AvLaval, double &APLaval, double &ATLaval, double &ACsonic) : TState(), vLaval(AvLaval), PLaval(APLaval), TLaval(ATLaval), Csonic(ACsonic)  {

}

void TSonicState::SetSonicState(TSonicState AsoState) {
    v = AsoState.v ;
    P = AsoState.P ;
    T = AsoState.T ;
    h = AsoState.h ;
    s = AsoState.s ;
    x = AsoState.x ;
    vLaval = AsoState.vLaval ;
    PLaval = AsoState.PLaval ;
    TLaval = AsoState.TLaval ;
    Csonic = AsoState.Csonic ;
}

void TSonicState::GetSonicState(TSonicState AsoState) {
    AsoState.v = v ;
    AsoState.P = P ;
    AsoState.T = T ;
    AsoState.h = h ;
    AsoState.s = s ;
    AsoState.x = x ;
    AsoState.vLaval = vLaval ;
    AsoState.PLaval = PLaval ;
    AsoState.TLaval = TLaval ;
    AsoState.Csonic = Csonic ;
}

void TSonicState::SetVL(double ivl) {
vLaval=ivl;
}

void TSonicState::SetPL(double iPl) {
PLaval=iPl;
}

void TSonicState::SetTL(double iTl) {
TLaval=iTl;
}

void TSonicState::SetCs(double iCs) {
Csonic=iCs;
}

double TSonicState::GetVL() {
    return vLaval;
}

double TSonicState::GetPL() {
    return PLaval;
}

double TSonicState::GetTL() {
    return TLaval;
}

double TSonicState::GetCs() {
    return Csonic;
}

