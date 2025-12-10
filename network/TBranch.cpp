//
// Created by lucile.schulthe on 16.05.2022.
// Modified by cornelia.blanke on 06.03.2023.
//

#include "TBranch.h"



TBranch::TBranch(){}
TBranch::TBranch(int l) : l(l) {}
TBranch::~TBranch(){}

void TBranch::init (){}


// BASICES
/*
BasicES & TBranch::performance(){}
float TBranch::objective(TSequence& seq){}
*/

void TBranch::display(){
    cout << "\n";
    cout <<"/****************************************************************************/\n";
    cout<<"TBranch                               : \n\n";
    cout<<"Quantum number l                      :       "<<l<<"\n";
    cout<<"Number of Segments                    :       "<< segnumber <<"\n";
    cout<<"Number of Substations                 :       "<< sstnumber <<"\n";
    cout<< "\n";
    cout <<"/****************************************************************************/\n";
    cout << endl;
}

void TBranch::compute_m_from_table() {
    if (QNet->mtable.empty()){
        cerr << "mtable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int)(QNet->nbBranches/2);
        if (QNet->networktable[1][branch0 + l] <= 0)   //first segment followed by neutron?!
            m = QNet->mtable[QNet->nbRecords - 1][branch0 + l];
        else    //first segment followed by proton
            m = QNet->mtable[1][branch0 + l];
    }
}

unsigned int TBranch::compute_m(int l_number){
    unsigned int sum=0;
    int branch0 = (int)(QNet->nbBranches/2);
    for (int i = 1; i < QNet->nbRecords - 1; i++) {
        if (QNet->networktable[i][branch0 + l_number] > 0)
            sum++;
        else if (QNet->networktable[i][branch0 + l_number] == -1){
            int next_branch = QNet->networktable[QNet->nbRecords - 1][branch0 + l_number];
            sum += compute_m(-next_branch) + compute_m(next_branch);
        }
        else
            break;
    }
    return sum;
}

void TBranch::compute_m() {
    m = compute_m(l);
}

void TBranch::compute_sstnumber(){
    int branch0 = (int)(QNet->nbBranches/2);
    sstnumber = 0;
    for(int i=1; i<QNet->nbRecords-1; i++){
        if(QNet->networktable[i][branch0 + l]> 0)
            sstnumber++;
        else
            break;
    }
}
void TBranch::compute_sstnumber_from_table(){
    if (QNet->mtable.empty()){
        cerr << "mtable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int)(QNet->nbBranches/2);
        sstnumber = QNet->mtable[0][branch0 + l];
    }
}

void TBranch::compute_segnumber(){
    compute_sstnumber();
    int branch0 = (int)(QNet->nbBranches/2);
    if (QNet->networktable[QNet->nbRecords - 1][branch0 + l] > 0)     // neutron below?
        segnumber = sstnumber + 1;
    else
        segnumber = sstnumber;
}
void TBranch::compute_segnumber_from_table(){
    compute_sstnumber_from_table();
    int branch0 = (int)(QNet->nbBranches/2);
    if (QNet->networktable[QNet->nbRecords - 1][branch0 + l] > 0)     // neutron below?
        segnumber = sstnumber + 1;
    else
        segnumber = sstnumber;
}

void TBranch::get_segments(){
    if (QNet->Table_of_Segments.empty()){
        cerr << "Table_of_Segments not defined" << endl;
        exit(1);
    }
    Seg_on_Branch.clear();
    int branch0 = (int)(QNet->nbBranches/2);
    for (int i=0; i<QNet->Table_of_Segments.size(); i++){
        if (QNet->Table_of_Segments[i][l + branch0] != nullptr)
            Seg_on_Branch.push_back( QNet->Table_of_Segments[i][l + branch0] );
    }
    segnumber = Seg_on_Branch.size();
}

void TBranch::get_protons(){
    if (QNet->List_of_Protons.empty()){
        cerr << "List_of_Protons not defined" << endl;
        exit(1);
    }
    Prot_on_Branch.clear();
    int branch0 = (int)(QNet->nbBranches/2);
    for (int i=1; i<QNet->nbRecords-1; i++) {
        int id = QNet->networktable[i][l + branch0];
        if (id > 0) {
            if (id != QNet->List_of_Protons[id - 1]->ID) {
                cerr << "Mismatch in List_of_Protons" << endl;
                exit(1);
            }
            else
                Prot_on_Branch.push_back(QNet->List_of_Protons[id - 1]);
        }
        else
            break;
    }
    sstnumber = Prot_on_Branch.size();
}

void TBranch::get_mu_from_segment(){
    if (QNet->Table_of_Segments.empty()){
        cerr << "Table_of_Segments not defined" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    mu = QNet->Table_of_Segments[0][branch0 + l]->mu;
}

void TBranch::compute_mu_from_table() {
    if (QNet->muTable.empty()){
        cerr << "muTable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int)(QNet->nbBranches/2);
        mu = QNet->muTable[0][branch0 + l];
    }
}

void TBranch::compute_massflow_from_mu(){
    double Mdot = massflow_from_mu(mu, QNet->CentralPlant->max_massflow, QNet->Z);
    flow->massflow(Mdot);
}

void TBranch::compute_mu_from_massflow(){
    mu = mu_from_massflow(flow->massflow(), QNet->CentralPlant->max_massflow, QNet->Z);
}

void TBranch::connect_network() {
    if (QNet->networktable.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    int next_neutron = QNet->networktable[QNet->nbRecords-1][branch0 + l];
    if (next_neutron > 0)
        next.push_back(QNet->List_of_Neutrons[next_neutron-1]);
    else
        next.push_back(nullptr);
}

void TBranch::connect_flow() {
    if (QNet->networktable.empty() || QNet->List_of_Neutrons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    int next_neutron = QNet->networktable[QNet->nbRecords-1][branch0 + l];
    if (next_neutron > 0)
        nextflow.push_back(QNet->List_of_Neutrons[next_neutron-1]->flow);
    else
        nextflow.push_back(nullptr);
}

void TBranch::init_Tambient(const bool below) {
    if (below)
        Tambient = &(QNet->Tground);
    else
        Tambient = &(QNet->Toutside);
}

void TBranch::init_Tref() {
    Tref = &(QNet->T0);
}

void TBranch::compute_dP(){
    dP = 0.0;
    double length = 0.0;
    for (int i=0; i<segnumber; i++){
        dP += Seg_on_Branch[i]->dP;
        length += Seg_on_Branch[i]->Pipe.Length;
    }
    for (int i=0; i<sstnumber; i++){
        dP += Prot_on_Branch[i]->dP;
    }
    dP_m = dP / (2 * length);       // two pipes !
}

void TBranch::compute_tau(){
    tau = 0.0;
    for (int i=0; i<segnumber; i++){
        tau += Seg_on_Branch[i]->tau;
    }
    for (int i=0; i<sstnumber; i++){
        tau += Prot_on_Branch[i]->tau;
    }
}

void TBranch::estimate_powerloss() {
    dPower = 0.0;
    for (int i=0; i<segnumber; i++){
        dPower += Seg_on_Branch[i]->dPower;
    }
    for (int i=0; i<sstnumber; i++){
        dPower += Prot_on_Branch[i]->dPower;
    }
}

void TBranch::get_TIn_TOut_dT() {
        flow->pStateIn->T = Seg_on_Branch[0]->flow->TIn();
        flow->pStateOut->T = Seg_on_Branch[0]->flow->TOut();
        dT = Seg_on_Branch[segnumber-1]->dT;
}

void TBranch::compute_material_cost() {
    material_cost = 0.0;
    for (int i=0; i<segnumber; i++){
        material_cost += Seg_on_Branch[i]->material_cost;
    }
    for (int i=0; i<sstnumber; i++){
        material_cost += Prot_on_Branch[i]->material_cost;
    }
}

void TBranch::compute_engineering_cost(const char* region) {
    engineering_cost = 0.0;
    for (int i=0; i<segnumber; i++){
        engineering_cost += Seg_on_Branch[i]->engineering_cost;
    }
    for (int i=0; i<sstnumber; i++){
        engineering_cost += Prot_on_Branch[i]->engineering_cost;
    }
}