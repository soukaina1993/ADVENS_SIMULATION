//
// Created by lucile.schulthe on 16.05.2022.
// Modified by cornelia.blanke on 02.03.2023.
//

#ifndef FLUIDS_TSEGMENT_H
#define FLUIDS_TSEGMENT_H


//#include "../Flow/Flow.h"
//#include "Quantum.h"
#include "QNetwork.h"
#include "../BasicEnergySystem/pipe/HPipe.h"
#include "../BasicEnergySystem/pipe/HBendPipe.h"

class QNetwork;

class TSegment: public Quantum	  {
public:
    void init() override;

    // Fonction non implementée mais pour avoir l'héritage pour l'instant
    BasicES & performance();
    float objective(TSequence& seq);
    void display() override;

public:
    TSegment();
    TSegment(int l, unsigned int k);
    ~TSegment() override;

    unsigned int compute_m(int l_number, unsigned int k_number=1);
    void compute_m_from_table() override;    // get m from mtable, if exists
    void compute_m() override;
    void compute_mu_from_table() override;

    void compute_massflow_from_mu() override;
    void compute_mu_from_massflow() override;

    void connect_network() override;
    void connect_flow() override;

    void init_DN();     //inits DN and pipetype
    void init_length();
//    void init_roughness(double irough=0.0);
    void init_pipeprop();
    void init_Tambient(bool below) override;
    void init_Tref() override;
    void compute_dP() override;
    void compute_tau() override;

    void compute_material_cost() override;
    void compute_engineering_cost(const char* region) override;

    QNetwork* QNet= nullptr;

    //HPipe* Pipe = new HPipe(flow);
    HPipe Pipe = HPipe(flow);       // HPipe and TSegment share same flow
//    unsigned int nbElbows=0;      TODO
//    HBendPipe Elbow = HBendPipe(Pipe);    // Are all elbows the same? Then we need only one here!

    int l;              // quantum number l (constant)
    unsigned int k;     // quantum number k (constant)

    double* pipetype=nullptr;      // pointer to a line in QNet->pipetable, call pipetype[i] to get variable values
    int DN=0;
    double length=1.0;
    //double roughness=0.0;     constant for all pipes?

};


#endif //FLUIDS_TSEGMENT_H
