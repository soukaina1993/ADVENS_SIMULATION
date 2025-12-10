//
// Created by lucile.schulthe on 13.06.2022.
//

#ifndef FLUIDS_LEEK1_H
#define FLUIDS_LEEK1_H

#include "../math/maths.h"
#include <string>
#include "../iofile/iofile.h"
#include "FluidCore.h"



/*
struct LKFluidRec{ // plus besoin, passer directement dans la classe sans struct
    string fluidName, fluidSymbole, fuildFormula, fluidSource;
    double Tfp, Tnb, Pc, Vc, Tc, MW, ALPHA, BETA, GAMMA, DELTA, roL, TL;
    double wPitzer, Dipm;
    double ADA0, ADA1, Tref, WP, AF, ZC, HREFO, SREFO;
    double SC, HC;
};*/



class LeeK1 {

public:
    //pourra être suppriment quand derive de fluidcore

// decalrer dans fluidcore
//    double Tref,href,sref,r;
//    char *FluidName;
//    double vc,Pc,Tc,uc,hc,sc,Zc, Tnb; //donn�es en conditions critiques
//    double Mw,DipM;
//    double w, k, Mur;
//    double E(int i);
//    double B(int i);
//






    // depuis: ReadAllFluids
 //   void init_LKConst(LKConstRec &LKConstRec_st);

// pas sur que ce soit utile!
 //   LKConstRec LKConst;
 //   LKFluidRec * LKFluid;
 //   vector<LKFluidRec> * ArrayOfFluid;



    // Cpo coefficients




};







#endif //FLUIDS_LEEK1_H
