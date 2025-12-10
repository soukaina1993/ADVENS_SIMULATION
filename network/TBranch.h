//
// Created by lucile.schulthe on 16.05.2022.
// Modified by cornelia.blanke on 06.03.2023.
//

#ifndef FLUIDS_TBRANCH_H
#define FLUIDS_TBRANCH_H

#include "QNetwork.h"
//#include "../consumer/sst.h"

//class sst;
class QNetwork;
class TSegment;
class TProton;

class TBranch: public Quantum	  { // Branch is a part of quantum network
public:
    TBranch();
    explicit TBranch(int l);
    ~TBranch() override;
    void init () override;
    void display() override;

/*
    vector<double> getMassflow();

    double getMassflow(int Time);
*/


/*    // Fonction non implementée mais pour avoir l'héritage pour l'instant

    BasicES & performance();
    float objective(TSequence& seq);
    void display();
*/


public:
    unsigned int compute_m(int l_number);
    void compute_m() override;
    void compute_m_from_table() override;        // total number of sst alimented by branch
    void compute_sstnumber();           // number of sst on that branch
    void compute_sstnumber_from_table();
    void compute_segnumber();           // number of segments on that branch
    void compute_segnumber_from_table();
    void get_segments();
    void get_protons();
    void get_mu_from_segment();
    void compute_mu_from_table() override;       // total mu in branch

    void compute_massflow_from_mu() override;
    void compute_mu_from_massflow() override;

    void connect_network() override;
    void connect_flow() override;

    void init_Tambient(bool below) override;
    void init_Tref() override;


    void compute_dP() override;      // addition of dP of all segments
    void compute_tau() override;     // addition of tau of all segments
    void estimate_powerloss() override; // addition of dPower of all segments
    void get_TIn_TOut_dT();      // get TIn/TOut of first segment, dT of last segment on branch

    void compute_material_cost() override;      // addition of all segments
    void compute_engineering_cost(const char* region) override;

    QNetwork* QNet;

    int l;
    unsigned int sstnumber=0, segnumber=0;
    vector<TSegment*> Seg_on_Branch;    // list of pointers to all segments on branch, length is segnumber
    vector<TProton*> Prot_on_Branch;
    //double dP; // pressure loss


};


#endif //FLUIDS_TBRANCH_H
