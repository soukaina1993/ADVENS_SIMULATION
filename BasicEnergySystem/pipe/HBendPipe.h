//
// Created by lucile.schulthe on 24.11.2022.
//

#ifndef TESTCOMPOSANT_CPP_HBENDPIPE2_H
#define TESTCOMPOSANT_CPP_HBENDPIPE2_H


#include "HPipe.h"

class HBendPipe: public HPipe {
public:
    HBendPipe();
    HBendPipe(const char *iPipeName, TInputFile &InputFile);
    HBendPipe(HPipe* pipe);
    HBendPipe(HPipe& pipe);
    HBendPipe(Flow* iflow);
    HBendPipe(Flow_physic* iflow);
    ~HBendPipe();

    double angle=0.0, bendradius=0.0;

    void PressureDropCalculation();
    double Heattransfer_conv();
};


#endif //TESTCOMPOSANT_CPP_HBENDPIPE2_H
