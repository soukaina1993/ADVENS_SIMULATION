//
// Created by lucile.schulthe on 03.06.2022.
//

#ifndef FLUIDS_TCONNEC_H
#define FLUIDS_TCONNEC_H


#include "../BasicEnergySystem/BasicES.h"
#include "../Flow/Flow.h"


enum Sort{HEAT=0, COOL=1, BOTH=2};
class TConnec: public BasicES {
public:


    TConnec(){};

    void calcul_hotflow();
    void calcul_coldflow();

    void calcul_hotflow(double *&Cold_Massflow);
    void calcul_coldflow(double *&Hot_Massflow);

    void sethotflow(Flow &flow){Hot=flow;};
    void setcoldflow(Flow &flow){Cold=flow;};

    Flow gethotflow(){return Hot;};
    Flow getcoldflow(){return Cold;};

    virtual double getloss(){return 0;};
    virtual double getMassFlow(){};

    virtual void   init (){};
    virtual BasicES & performance(){};
    virtual float objective(TSequence& seq){};
    virtual void display(){};

    Flow Hot;
    Flow Cold;

    double efficiency=0,th_efficiency=0.9;

   Sort    sor;
    double dP;
};


#endif //FLUIDS_TCONNEC_H
