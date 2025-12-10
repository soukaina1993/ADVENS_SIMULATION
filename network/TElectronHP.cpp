//
// Created by cornelia.blanke on 25.10.2023.
//

#include "TElectronHP.h"


void TElectronHP::create_heatpump(string fluidin1, double ieta_comp, double ipinchCond, double ipinchEvap) {
    HP = new THeatpump(fluidin1, ieta_comp, ipinchCond, ipinchEvap);
    Heatpump_created=true;
}

void TElectronHP::create_heatpump(const char* fluidin1, double ieta_comp, double ipinchCond, double ipinchEvap) {
    string fluid;
    fluid = fluidin1;
    HP = new THeatpump(fluid, ieta_comp, ipinchCond, ipinchEvap);
    Heatpump_created=true;
}

void TElectronHP::compute_power_from_district(){
    double factor = 1;        //TODO heat loss

    // DeltaQ = lossfactor * (Tin - Text);   lossfactor in [W/K]
    // flow->Power = demand + DeltaQ;

    if(D.power.empty()){
        cerr << "District not defined/computed" << endl;
        exit(1);
    }
    demand = D.power[QNet->timestep];         // what is needed by district
    if (HP->COP == 0){
        cerr << "COP not calculated" << endl;
        exit(1);
    }
    flow->Power = factor * (1 - 1/HP->COP) * demand;
}

void TElectronHP::set_COP(const double COP) {
    if (COP <= 0.0){
        cerr << "A constant COP value of " << COP << " is not allowed" << endl;
        exit(1);
    }
    HP->COP = COP;
    useCOPconst = true;
}

void TElectronHP::compute_COP(double Thot) {
    if (!useCOPconst){
        double Tcold = *Tambient + (flow->TIn() - *Tambient) * exp(-tau);   // temperature at HP entry
        HP->COP_real(Thot, Tcold);
    }
}

void TElectronHP::compute_COP() {
    if (!useCOPconst)
        compute_COP(D.Tneed_max);
}

void TElectronHP::compute_material_cost(){
    // connecting pipes
    const char* fluidname = flow->fluid->FluidName;
    if (strcmp(fluidname, "H2O") == 0) {
        // project report of Jeremy Rolle page 25; good agreement with Briguet
        material_cost = 2 * length * (0.000471 * pow(DN, 2) + 0.8862 * DN + 57.412);
    }
    else if (strcmp(fluidname, "CO2") == 0){
        // master thesis Briguet page 24
        double alpha = 1.1;
        double DN_vapour=DN, DN_liquid=DN;  //TODO two different (!) pipes
        material_cost = length * 35 * (pow(DN_vapour, alpha) + pow(DN_liquid, alpha)) / pow(25, alpha);
    }
    else {
        cerr << "Costs of "<< fluidname << " pipes not implemented" << endl;
    }

    // heat pump (Jeremy Rolle page 29)
    double power = D.power_max/1000;    // in kW TODO (what/where is the definition of power of SST?)
    material_cost += (19.569 * log(power) - 35.018) * 1000;
}

void TElectronHP::compute_engineering_cost(const char* region) {
    // connecting pipes
    const char* fluidname = flow->fluid->FluidName;
    if (strcmp(fluidname, "H2O") == 0) {
        // master thesis Briguet page 24
        double excavation;
        if (strcmp(region, "Valais") == 0)
            excavation = length * (2.3414 * DN + 269.94);
        else if (strcmp(region, "Geneva") == 0)
            excavation = length * (4.8923 * DN + 645.48);
        else        // mean value
            excavation = length * (3.6169 * DN + 457.71);
        engineering_cost = excavation;
    }
    else if (strcmp(fluidname, "CO2") == 0){
        // master thesis Briguet page 24
        double excavation = length * (1.38 * DN + 200);          //valid for Valais
        engineering_cost = excavation;
    }
    else {
        cerr << "Costs of "<< fluidname << " pipes not implemented" << endl;
    }

    // heat pump (Jeremy Rolle page 29)
    double power = D.power_max/1000;    // in kW
    engineering_cost += (0.6376 * log(power) - 0.2742) * 1000;
}
