//
// Created by lucile.schulthe on 03.06.2022.
//

#include "TConnec.h"



void TConnec::calcul_hotflow(){
    Hot.Massflow=0;
    getMassFlow();
};


void TConnec::calcul_coldflow(){
    Cold.Massflow=0;
    getMassFlow();
};


void  TConnec::calcul_hotflow(double*&Cold_Massflow){
    Hot.Massflow=0;
    Cold.Massflow=Cold_Massflow;
    getMassFlow();
}
void  TConnec::calcul_coldflow(double *&Hot_Massflow){
    Hot.Massflow=Hot_Massflow;
    Cold.Massflow=0;
    getMassFlow();
}