//
// Created by cornelia.blanke on 06.03.2023.
//

#include <stdio.h>
#include <algorithm>
#include "TNeutron.h"

TNeutron::TNeutron() {}
TNeutron::TNeutron(unsigned int n): n(n) {}
TNeutron::~TNeutron() {}

void TNeutron::display(){
    cout << "\n";
    cout <<"/****************************************************************************/\n";
    cout<<"TNeutron                              : \n\n";
    cout<<"Quantum number n                      :       "<< n <<"\n";
    cout<< "\n";
    cout <<"/****************************************************************************/\n";
    cout << endl;
}

void TNeutron::compute_m_from_table() {
    if (QNet->mtable.empty()){
        cerr << "mtable not defined" << endl;
        exit(1);
    }
    else {
        for (int j = 0; j < QNet->nbBranches; j++){
            if( QNet->networktable[QNet->nbRecords-1][j] == n){
                m = QNet->mtable[QNet->nbRecords-1][j];
                break;
            }
        }
    }
}

unsigned int TNeutron::compute_m(unsigned int n_number){
    unsigned int sum=0;
    int branch0 = (int)(QNet->nbBranches/2);

    for (int i = 1; i < QNet->nbRecords - 1; i++) {
        if (QNet->networktable[i][branch0 + n_number] > 0)     // positive branch
            sum++;
        if (QNet->networktable[i][branch0 - n_number] > 0)     // negative branch
            sum++;
    }
    unsigned int pos_next_branch = QNet->networktable[QNet->nbRecords - 1][branch0 + n_number];
    if (pos_next_branch > 0)
        sum += compute_m(pos_next_branch);
    unsigned int neg_next_branch = QNet->networktable[QNet->nbRecords - 1][branch0 - n_number];
    if (neg_next_branch > 0)
        sum += compute_m(neg_next_branch);

    return sum;
}

void TNeutron::compute_m() {
    m = compute_m(n);
}

void TNeutron::compute_mu_from_table() {
    if (QNet->muTable.empty()){
        cerr << "muTable not defined" << endl;
        exit(1);
    }
    else {
        for (int j = 0; j < QNet->nbBranches; j++){
            if( QNet->networktable[QNet->nbRecords-1][j] == n) {
                mu = QNet->muTable[QNet->nbRecords - 1][j];
                break;
            }
        }
    }
}

void TNeutron::compute_massflow_from_mu(){
    double Mdot = massflow_from_mu(mu, QNet->CentralPlant->max_massflow, QNet->Z);
    flow->massflow(Mdot);
}

void TNeutron::compute_mu_from_massflow(){
    mu = mu_from_massflow(flow->massflow(), QNet->CentralPlant->max_massflow, QNet->Z);
}

void TNeutron::connect_network() {
    if (QNet->Table_of_Segments.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    // negative segment is next[0], positive segment is next[1]
    next.push_back(QNet->Table_of_Segments[0][branch0 - n]);
    next.push_back(QNet->Table_of_Segments[0][branch0 + n]);
}

void TNeutron::connect_flow() {
    if (QNet->Table_of_Segments.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    // negative segment is nextflow[0], positive segment is nextflow[1]
    nextflow.push_back(QNet->Table_of_Segments[0][branch0 - n]->flow);
    nextflow.push_back(QNet->Table_of_Segments[0][branch0 + n]->flow);
}

void TNeutron::init_Tambient(const bool below) {
    if (below)
        Tambient = &(QNet->Tground);
    else
        Tambient = &(QNet->Toutside);
}

void TNeutron::init_Tref() {
    Tref = &(QNet->T0);
}

void TNeutron::compute_tau(){}