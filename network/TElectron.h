//
// Created by cornelia.blanke on 06.03.2023.
//

#ifndef FLUIDS_TELECTRON_H
#define FLUIDS_TELECTRON_H

//#include "Quantum.h"
#include "QNetwork.h"
#include "../BasicEnergySystem/pipe/HPipe.h"
#include "../consumer/District.h"


class QNetwork;

class TElectron: public Quantum {
public:

    enum Balance {NEUTRAL=0, POS=1, NEG=2};     //not yet used

    //TElectron();
    explicit TElectron(TElectron::Balance balance=TElectron::NEUTRAL);
    explicit TElectron(unsigned int ID, TElectron::Balance balance=TElectron::NEUTRAL);
    TElectron(unsigned int ID, int l, unsigned int k, TElectron::Balance balance=TElectron::NEUTRAL);
    ~TElectron() override;

    void display() override;

    // determine position on networktable
    void find_kl_from_ID();

public:
    void compute_m_from_table() override;
    void compute_m() override;
    void compute_mu_from_table() override;

    void compute_massflow_from_mu() override;
    void compute_mu_from_massflow() override;

    void connect_network() override;
    void connect_flow() override;

    // pipe functions
    void init_DN();     //inits DN and pipetype
//    void init_length();
//    void init_roughness(double irough=0.0);
    void init_pipeprop();
    void init_Tambient(bool below) override;
    void init_Tref() override;
    void update_dT();       // if dT time-dependent
    void compute_dP() override;
    void compute_tau() override;

    void compute_T_supply() override;
    void compute_T_return() override;

    void compute_material_cost() override;
    void compute_engineering_cost(const char* region) override;
    void compute_operating_cost() override;
    void compute_resources_cost() override;     // resources needed for that SST

    // district functions
    virtual void compute_power_from_district();     //overridden for TElectronHP

    void compute_mu_from_power();

    QNetwork* QNet;

    unsigned int ID=0, k;
    int l;

    vector<double> dTvar;           // holds time-dependent dT values, if needed

    HPipe Pipe = HPipe(flow);   // Raccordement
    double* pipetype=nullptr;      // pointer to a line in QNet->pipetable, call pipetype[i] to get variable values
    int DN=-1;
    double length=-1.0;
    //double roughness=0.0;

    District D;
};


#endif //FLUIDS_TELECTRON_H
