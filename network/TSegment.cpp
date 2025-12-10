//
// Created by lucile.schulthe on 16.05.2022.
// Modified by cornelia.blanke on 02.03.2023.
//

#include "TSegment.h"


TSegment::TSegment()=default;
TSegment::TSegment(int l, unsigned int k): l(l), k(k) {}
TSegment::~TSegment()=default;

void 	      TSegment::init (){}
BasicES & TSegment::performance(){ }
float TSegment::objective(TSequence& seq){ return 0; }

void TSegment::display(){
    cout << "\n";
    cout <<"/****************************************************************************/\n";
    cout<<"TSegment                              : \n\n";
    cout<<"Quantum numbers l, k                  :       "<<l<<", "<<k<<"\n";
    cout<<"Pipe Properties                       :       "<<"DN"<<DN<<", Length "<<length<<"m\n";
    cout<< "\n";
    cout <<"/****************************************************************************/\n";
    cout << endl;
}

void TSegment::compute_m_from_table() {
    if (QNet->mtable.empty()){
        cerr << "mtable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int)(QNet->nbBranches/2);
        if (QNet->networktable[k][branch0 + l] > 0)   //segment followed by proton
            m = QNet->mtable[k][branch0 + l];
        else    //segment followed by neutron
            m = QNet->mtable[QNet->nbRecords - 1][branch0 + l];
    }
}

unsigned int TSegment::compute_m(int l_number, unsigned int k_number){
    unsigned int sum=0;
    int branch0 = (int)(QNet->nbBranches/2);
    int next_branch = QNet->networktable[QNet->nbRecords - 1][branch0 + l_number];
    for (unsigned int i = k_number; i < QNet->nbRecords - 1; i++) {
        if (QNet->networktable[i][branch0 + l_number] > 0)
            sum++;
        else if (next_branch > 0){
            sum += compute_m(-next_branch) + compute_m(next_branch);
            next_branch = 0;
        }
        else
            break;
    }
    return sum;
}

void TSegment::compute_m() {
    m = compute_m(l,k);
}

void TSegment::compute_mu_from_table() {
    if (QNet->muTable.empty()){
        cerr << "muTable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int)(QNet->nbBranches/2);
        mu = QNet->muTable[QNet->nbRecords - 1][branch0 + l];   // neutron below segment
        for (unsigned int i=QNet->nbRecords-2; i>=k; i--)                // add protons
            mu += QNet->muTable[i][branch0 + l];
    }
}

void TSegment::compute_massflow_from_mu(){
    double Mdot = massflow_from_mu(mu, QNet->CentralPlant->max_massflow, QNet->Z);
    flow->massflow(Mdot);
}

void TSegment::compute_mu_from_massflow(){
    mu = mu_from_massflow(flow->massflow(), QNet->CentralPlant->max_massflow, QNet->Z);
}

void TSegment::connect_network() {
    if (QNet->networktable.empty() || QNet->List_of_Protons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    int next_proton = QNet->networktable[k][branch0 + l];
    if (k < QNet->nbRecords-1 && next_proton > 0)    // proton below
        next.push_back(QNet->List_of_Protons[next_proton-1]);
    else    // neutron below
        next.push_back(QNet->List_of_Neutrons[QNet->networktable[QNet->nbRecords-1][branch0 + l] - 1]);
}

void TSegment::connect_flow() {
    if (QNet->networktable.empty() || QNet->List_of_Protons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    int next_proton = QNet->networktable[k][branch0 + l];
    if (k < QNet->nbRecords-1 && next_proton > 0)    // proton below
        nextflow.push_back(QNet->List_of_Protons[next_proton-1]->flow);
    else    // neutron below
        nextflow.push_back(QNet->List_of_Neutrons[QNet->networktable[QNet->nbRecords-1][branch0 + l] - 1]->flow);
}


void TSegment::init_DN(){
    if (k == 0){
        cerr << "Init k first " << endl;
        exit(1);
    }
    if (QNet->diamtable.empty()){
        cerr << "Read diamtable first " << endl;
        exit(1);
    }
    DN = QNet->diamtable[k-1][QNet->N + l];
    if (!QNet->pipetable.empty()){
        for (int i=0; i<QNet->pipetable.size(); i++) {
            if ((int)QNet->pipetable[i][0] == DN) {
                pipetype = &(QNet->pipetable[i][0]);
                break;
            }
        }
    }
}

void TSegment::init_length(){
    if (k == 0){
        cerr << "Init k first " << endl;
        exit(1);
    }
    if (QNet->lengthtable.empty()){
        cerr << "Read lengthtable first " << endl;
        exit(1);
    }
    length = QNet->lengthtable[k-1][QNet->N + l];
}

/*void TSegment::init_roughness(double irough){
    roughness = irough;
}*/

void TSegment::init_pipeprop(){
        Pipe.Length=length;
        if (length==0 || DN==0 || pipetype==nullptr){
            cerr << "Init length, DN and pipetype first" << endl;
            exit(1);
        }
        else if (QNet->pipetable[0].size() ==8){     // two-layer pipe
            Pipe.Diam = pipetype[1];
            Pipe.RugAbs = pipetype[7];
            if (pipetype[6] != 0)
                Pipe.U_pipe = pipetype[6];
            else {
                Pipe.U_pipe = Pipe.Heattransfer_mat(Pipe.Diam, pipetype[2], pipetype[3], pipetype[4], pipetype[5]);
                pipetype[6] = Pipe.U_pipe;  // update pipetable
            }
        }
        else {                          // three-layer pipe
            Pipe.Diam = pipetype[1];
            Pipe.RugAbs = pipetype[9];
            if (pipetype[8] != 0)
                Pipe.U_pipe = pipetype[8];
            else {
                Pipe.U_pipe = Pipe.Heattransfer_mat(Pipe.Diam, pipetype[2], pipetype[3], pipetype[4], pipetype[5], pipetype[6], pipetype[7]);
                pipetype[8] = Pipe.U_pipe;  // update pipetable
            }
        }
}

void TSegment::init_Tambient(const bool below) {
    if (below)
        Tambient = &(QNet->Tground);
    else
        Tambient = &(QNet->Toutside);
}

void TSegment::init_Tref() {
    Tref = &(QNet->T0);
}

void TSegment::compute_dP(){
    Pipe.PressureDropCalculation();
    dP_m = Pipe.DeltaPr_m;
    dP = 2 * Pipe.DeltaPr;          // two pipes !!
    speed = Pipe.Speed;
}

void TSegment::compute_tau(){
    Pipe.tauCalculation();
    tau = Pipe.tau;
}

void TSegment::compute_material_cost(){
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
}

void TSegment::compute_engineering_cost(const char* region) {
    const char* fluidname = flow->fluid->FluidName;
    if (strcmp(fluidname, "H2O") == 0) {
        // master thesis Briguet page 24
        double excavation;      // other types of engineering costs may be added
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
}