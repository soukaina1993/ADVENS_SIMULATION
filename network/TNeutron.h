//
// Created by cornelia.blanke on 06.03.2023.
//

#ifndef FLUIDS_TNEUTRON_H
#define FLUIDS_TNEUTRON_H

//#include "Quantum.h"
#include "QNetwork.h"

class QNetwork;

class TNeutron: public Quantum {
public:
    TNeutron();
    TNeutron(unsigned int n);
    ~TNeutron();
    void display() override;

public:
    unsigned int compute_m(unsigned int n_number);
    void compute_m_from_table() override;    // get m from mtable, if exists
    void compute_m() override;               // compute m from networktable
    void compute_mu_from_table() override;

    void compute_massflow_from_mu() override;
    void compute_mu_from_massflow() override;

    void connect_network() override;
    void connect_flow() override;
    void init_Tambient(bool below) override;
    void init_Tref() override;

    void compute_tau() override;


    QNetwork* QNet;      // TODO share network, pointer or reference?
    unsigned int n;

};


#endif //FLUIDS_TNEUTRON_H
