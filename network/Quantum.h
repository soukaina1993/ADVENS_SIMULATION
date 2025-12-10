//
// Created by cornelia.blanke on 06.03.2023.
//

#ifndef FLUIDS_QUANTUM_H
#define FLUIDS_QUANTUM_H

//#include <vector>
#include "../BasicEnergySystem/BasicES.h"
#include "../Flow/Flow_physic.h"
#include "../environment/environment.h"




class Quantum: public BasicES {
public:
    void init ();

    // Fonction non implementée mais pour avoir l'héritage pour l'instant
    BasicES & performance() override;
    float objective(TSequence& seq) override;
    void display() override;

public:
    unsigned int m=0;    // number of SST below
    bool below_ground=true;     // switch above/below ground
    double mu=0.0, mu_bar=0.0, delta=0.0, tau=0.0, nu_bar=0.0, epsilon=0.0;    // quantum numbers
    double speed=0.0;
    double demand=0.0;                  // total power demand of districts alimented by quantum
    double dT=0.0, dP=0.0, dP_m=0.0, dPower=0.0;  // differential temperature, pressure loss, pressure loss per meter, power loss
    double altitude=0.0;   //TODO
    double total_invest=0.0;
    double material_cost=0.0, engineering_cost=0.0, operating_cost=0.0,
            resources_cost=0.0, total_resources_cost=0.0;
    double *Tref=nullptr;
    double *Tambient=nullptr;
    double *Tsource=nullptr;


    Flow_physic* flow = new Flow_physic;         //TODO when to use Flow or Flow reference or Flow pointer?
    vector<Flow_physic*> nextflow;      // pointers to the flows of the next elements, size=2 for proton/neutron, 1 else
    vector<Quantum*> next;


    //FluidCore* fluid;       fluid is already a part of flow



// Constructors
    Quantum();
    ~Quantum() override;

    void init_fluid(FluidCore* &fluidin1);            // copied and modified (!) from Lucile
    void init_fluid(string fluidin1);

    virtual void init_Tambient(bool below)=0;
    virtual void init_Tref()=0;
    virtual void init_Tsource(Quantum* Source);

    virtual void compute_m_from_table()=0;     // get m from mtable, if exists
    virtual void compute_m()=0;                // compute m from networktable
    virtual void compute_mu_from_table()=0;    // get mu from muTable, if exists

    double massflow_from_mu(double mu, double max_source_massflow, unsigned int Z);
    virtual void compute_massflow_from_mu()=0;
    double mu_from_massflow(double Mdot, double max_source_massflow, unsigned int Z);
    virtual void compute_mu_from_massflow()=0;
    virtual void compute_mubar_from_mu();       // overridden for QNetwork (loop over all elements)

    virtual void connect_flow()=0;       //TODO probably not needed any more
    virtual void connect_network()=0;
    void set_nextTIn(double TIn);
    virtual double get_nextTIn();
    virtual double get_nextTOut();          // overridden for TProton (case of residual flow)

    void compute_dT_from_delta();
    void compute_delta_from_dT();

    virtual void compute_dP();
    virtual void compute_tau();

    //estimates dPower without computing the temperatures
    virtual void estimate_powerloss();      // overridden for QNetwork, TBranch

    virtual void compute_T_supply();        // overridden for TElectron
    virtual void compute_T_return();        // overridden for TElectron

    virtual void compute_nubar();           // overridden for QNetwork (loop over all elements)
    void compute_nubar_from_T(double k0);   // should only be used after T computation
    virtual void compute_epsilon();         // overridden for QNetwork (loop over all elements)

    void compute_massflow_from_power();
    virtual void compute_power_from_massflow();  // overridden for QNetwork (loop over all elements except TElectron)

    virtual void compute_material_cost();         // overridden for QNetwork (loop over all elements)
    virtual void compute_engineering_cost(const char* region);
    virtual void compute_total_invest();
    virtual void compute_operating_cost();
    virtual void compute_resources_cost();     // resources cost per time-step

};


#endif //FLUIDS_QUANTUM_H
