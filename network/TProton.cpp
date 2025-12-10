//
// Created by cornelia.blanke on 06.03.2023.
//

#include "TProton.h"

TProton::TProton() {}
TProton::TProton(unsigned int ID) : ID(ID) {}
TProton::TProton(unsigned int ID, int l, unsigned int k) : ID(ID), l(l), k(k) {}
TProton::~TProton() {}

void TProton::display(){
    cout << "\n";
    cout <<"/****************************************************************************/\n";
    cout<<"TProton                               :       ID "<<ID<<"\n\n";
    cout<<"Quantum numbers l, k                  :       "<<l<<", "<<k<<"\n";
    cout<< "\n";
    cout <<"/****************************************************************************/\n";
    cout << endl;
}

void TProton::compute_m_from_table() {
/*    if (QNet->mtable == NULL){
        cerr << "mtable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int)(QNet->nbBranches/2);
        m = QNet->mtable[k][branch0 + l];
    }*/
    m = 1;
}

unsigned int TProton::compute_m(int l_number, unsigned int k_number){
/*    unsigned int sum=0;
    int branch0 = (int)(QNet->nbBranches/2);
    for (int i = k_number; i < QNet->nbRecords - 1; i++) {
        if (QNet->networktable[i][branch0 + l_number] > 0)
            sum++;
        else
            break;
    }
    int next_branch = QNet->networktable[QNet->nbRecords - 1][branch0 + l_number];
    if (next_branch > 0)
        sum += compute_m(-next_branch) + compute_m(next_branch);
    return sum; */
    return 1;
}

void TProton::compute_m() {
//    m = compute_m(l, k);
    m = 1;
}

void TProton::compute_mu_from_table() {
    if (QNet->muTable.empty()){
        cerr << "muTable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int) (QNet->nbBranches / 2);
        mu = QNet->muTable[QNet->nbRecords - 1][branch0 + l];   // neutron below proton
        for (int i=QNet->nbRecords-2; i>=k; i--)                // add protons
            mu += QNet->muTable[i][branch0 + l];
    }
}

void TProton::compute_massflow_from_mu(){
    double Mdot = massflow_from_mu(mu, QNet->CentralPlant->max_massflow, QNet->Z);
    flow->massflow(Mdot);
}

void TProton::compute_mu_from_massflow(){
    mu = mu_from_massflow(flow->massflow(), QNet->CentralPlant->max_massflow, QNet->Z);
}

void TProton::connect_network() {
    if (QNet->Table_of_Segments.empty() || QNet->List_of_Electrons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    next.push_back(QNet->List_of_Electrons[ID-1]);     // connected electron is next[0]
    if (k<QNet->Table_of_Segments.size() && QNet->Table_of_Segments[k][branch0+l] != nullptr)        // following segment is next[1]
        next.push_back(QNet->Table_of_Segments[k][branch0 + l]);
//    else
//        next.push_back(nullptr);
}

void TProton::connect_flow() {
    if (QNet->Table_of_Segments.empty() || QNet->List_of_Electrons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    nextflow.push_back(QNet->List_of_Electrons[ID-1]->flow);     // connected electron is nextflow[0]
    if (k<QNet->Table_of_Segments.size() && QNet->Table_of_Segments[k][branch0+l] != nullptr)        // following segment is nextflow[1]
        nextflow.push_back(QNet->Table_of_Segments[k][branch0 + l]->flow);
//    else
//        nextflow.push_back(nullptr);
}

void TProton::init_Tambient(const bool below) {
    if (below)
        Tambient = &(QNet->Tground);
    else
        Tambient = &(QNet->Toutside);
}

void TProton::init_Tref() {
    Tref = &(QNet->T0);
}



void TProton::compute_tau(){}

double TProton::get_nextTOut(){
    if (next.empty()){
        cerr << "next-pointer not initialised" << endl;
        exit(1);
    }

    // TProton at open end of branch
    if (next.size()==1)
        return next[0]->flow->TOut();       //TODO add contribution from residual flow

    // TProton followed by TSegment
    else if (flow->massflow() != 0.0) {
        double sum=0.0;
        for (int i = 0; i < next.size(); i++)
            sum += next[i]->flow->TOut() * next[i]->flow->massflow();
        return sum/flow->massflow();       // mass-weighted average
    }
    else {
        //cout << "Warning: No massflow" << endl;
        double sum=0.0;
        for (int i = 0; i < next.size(); i++)
            sum += next[i]->flow->TOut();
        return sum/next.size();             // arithmetic average
    }
}