//
// Created by cornelia.blanke on 06.03.2023.
//

#ifndef FLUIDS_TPROTON_H
#define FLUIDS_TPROTON_H

//#include "Quantum.h"
#include "QNetwork.h"

class QNetwork;

class TProton: public Quantum {
public:
    TProton();
    TProton(unsigned int ID);
    TProton(unsigned int ID, int l, unsigned int k);
    ~TProton();

    void display() override;

public:
    unsigned int compute_m(int l_number, unsigned int k_number=1);   //starting from position (l_number,k_number)
    void compute_m_from_table() override;    // get m from mtable, if exists
    void compute_m() override;
    void compute_mu_from_table() override;

    void compute_massflow_from_mu() override;
    void compute_mu_from_massflow() override;

    void connect_network() override;
    void connect_flow() override;
    void init_Tambient(bool below) override;
    void init_Tref() override;

    void compute_tau() override;

    double get_nextTOut() override;


    QNetwork* QNet;
    unsigned int ID;
    int l;              // quantum number l (constant)
    unsigned int k;     // quantum number k (constant)

};


#endif //FLUIDS_TPROTON_H
