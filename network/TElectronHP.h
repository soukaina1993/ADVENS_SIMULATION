//
// Created by cornelia.blanke on 25.10.2023.
//

#ifndef FLUIDS_TELECTRONHP_H
#define FLUIDS_TELECTRONHP_H

#include "TElectron.h"
#include "../BasicEnergySystem/heatpump/THeatpump.h"


class TElectronHP : public TElectron {
public:
    explicit TElectronHP(TElectron::Balance balance=TElectron::NEG) : TElectron(balance){}
    explicit TElectronHP(unsigned int ID, TElectron::Balance balance=TElectron::NEG) : TElectron(ID, balance){}
    TElectronHP (unsigned int ID, int l, unsigned int k, TElectron::Balance balance=TElectron::NEG) :
            TElectron(ID, l, k, balance){}

    ~TElectronHP() override {
        if (Heatpump_created)
            delete HP;
    }



    void create_heatpump(string fluidin1, double ieta_comp=0.8, double ipinchCond=5, double ipinchEvap=5);
    void create_heatpump(const char* fluidin1="R134a", double ieta_comp=0.8, double ipinchCond=5, double ipinchEvap=5);
    void set_COP(double COP);       // set a constant COP
    void compute_COP(double Thot);  // compute COP from temperatures
    void compute_COP();
    void compute_power_from_district() override;

    void compute_material_cost() override;
    void compute_engineering_cost(const char* region) override;


public:
    THeatpump* HP= nullptr;

private:
    bool Heatpump_created=false;
    bool useCOPconst=false;

};


#endif //FLUIDS_TELECTRONHP_H
