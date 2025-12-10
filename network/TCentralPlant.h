//
// Created by cornelia.blanke on 06.03.2023.
//

#ifndef FLUIDS_TCENTRALPLANT_H
#define FLUIDS_TCENTRALPLANT_H

//#include "Quantum.h"
#include "QNetwork.h"


/*** Boiler class ***/
class Boiler{
public:
    enum Type {OIL=0, GAS=1, WOOD=2, PELLET=3, NONE=4};
    // use Type=NONE for energy sources like industry heat or waste incineration
    // or for the definition of new types
    Boiler() = default;
    Boiler(Boiler::Type type, double max_power);
    Boiler(Boiler::Type type, double max_power, double price);

public:
    Boiler::Type type =NONE;
    double max_power=0.0, res_cost_per_kWh=0.0;
};


/*** Price list for resources (oil, gas, ...) ***/
static unordered_map<Boiler::Type, double> PRICE_PER_kWh =
        {{Boiler::OIL,  0.09}, {Boiler::GAS, 0.11},
         {Boiler::WOOD, 0.05}, {Boiler::PELLET, 0.05},
         {Boiler::NONE, 0.0}};


/*** Central Plant class ***/
class QNetwork;

class TCentralPlant: public Quantum {
public:
    TCentralPlant();
    void display() override;

public:
    void compute_m_from_table() override;
    void compute_m() override;
    void compute_mu_from_table() override;

    void compute_massflow_from_mu() override;
    void compute_mu_from_massflow() override;

    void modify_pricelist(vector<Boiler::Type>& types, vector<double>& prices);   // override default values
    void modify_pricelist(unordered_map<Boiler::Type, double>& prices);
    void init_boilers(vector<Boiler::Type>& types, vector<double>& max_powers, vector<double>& prices);
    void init_boilers(vector<Boiler::Type>& types, vector<double>& max_powers);     // use default prices
    void init_boilers();   // use with settings coming from read_settings(...)
    void add_boiler(Boiler::Type typ, double maxpower);
    void compute_maxvalues();       // if one of them is missing
    void compute_maxvalues_from_boilers();
    void init_pressure(double Pin, double Pout);
    void init_temperature(double Tin, double Tout);

    void connect_network() override;
    void connect_flow() override;
    void init_Tambient(bool below) override;
    void init_Tref() override;

    void compute_tau() override;

    void compute_material_cost() override;
    void compute_operating_cost() override;
    void compute_resources_cost() override;
    void compute_res_cost_per_kWh();    // weighted average of boilers



    QNetwork* QNet;
    vector<Boiler> boilers;
    double max_massflow=0.0, max_power=0.0; // nominal specification of power plant
    double res_cost_per_kWh=0.0;
};


#endif //FLUIDS_TCENTRALPLANT_H
