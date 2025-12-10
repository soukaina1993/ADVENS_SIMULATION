//
// Created by cornelia.blanke on 06.03.2023.
//

#include "TElectron.h"


//TElectron::TElectron() {}
TElectron::TElectron(TElectron::Balance balance) {}
TElectron::TElectron(unsigned int ID, TElectron::Balance balance) : ID(ID) {}
TElectron::TElectron(unsigned int ID, int l, unsigned int k, TElectron::Balance balance) : ID(ID), l(l), k(k) {}
TElectron::~TElectron() = default;

void TElectron::display(){
    cout << "\n";
    cout <<"/****************************************************************************/\n";
    cout<<"TElectron                             :       ID "<<ID<<"\n\n";
    cout<<"Quantum numbers l, k                  :       "<<l<<", "<<k<<"\n";
    cout<<"Connecting Pipe                       :       "<<"DN"<<DN<<", Length "<<length<<"m\n";
    cout<< "\n";
    cout <<"/****************************************************************************/\n";
    cout << endl;
}

void TElectron::compute_m() {
    m = 1;
}

void TElectron::compute_m_from_table() {
    m = 1;
}

void TElectron::find_kl_from_ID(){
    int branch0 = (int)(QNet->nbBranches/2);
    for (int i=1; i<QNet->nbRecords-1; i++){
        for (int j=0; j<QNet->nbBranches; j++){
            if (QNet->networktable[i][j] == ID){
                k = i;
                l = j - branch0;
                break;
            }
        }
    }
}

void TElectron::compute_mu_from_table() {
    if (QNet->muTable.empty()){
        cerr << "muTable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int) (QNet->nbBranches / 2);
        mu = QNet->muTable[k][branch0 + l];
    }
}

void TElectron::compute_massflow_from_mu(){
    double Mdot = massflow_from_mu(mu, QNet->CentralPlant->max_massflow, QNet->Z);
    flow->massflow(Mdot);
}

void TElectron::compute_mu_from_massflow(){
    mu = mu_from_massflow(flow->massflow(), QNet->CentralPlant->max_massflow, QNet->Z);
}

void TElectron::connect_network() {
    if (QNet->List_of_Protons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    //next.push_back(QNet->List_of_Protons[ID-1]);
}

void TElectron::connect_flow() {
    if (QNet->List_of_Protons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    //nextflow.push_back(QNet->List_of_Protons[ID-1]->flow);
}

void TElectron::init_DN(){
    if (DN == -1){
        cerr << "Init DN first " << endl;
        exit(1);
    }
    if (DN == 0)
        return;
    if (QNet->pipetable.empty())
        cerr << "Read pipetable first" << endl;

    for (int i=0; i<QNet->pipetable.size(); i++) {
        if ((int)QNet->pipetable[i][0] == DN) {
            pipetype = &(QNet->pipetable[i][0]);
            break;
        }
    }
    if (pipetype==nullptr){
        cerr << "Pipe DN " << DN << " is not in the database" << endl;
        exit(1);
    }
}

void TElectron::init_pipeprop(){
    Pipe.Length=length;
    Pipe.Diam = DN;
    if (length==0)
        cout << "Warning: Pipe of substation "<< ID << " has length 0" << endl;
    if (DN == 0)
        cout << "Warning: Pipe of substation "<< ID << " has diameter 0" << endl;
    else if (DN==-1 || pipetype==nullptr){
        cerr << "Init DN and pipetype first" << endl;
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

void TElectron::init_Tambient(const bool below) {
    if (below)
        Tambient = &(QNet->Tground);
    else
        Tambient = &(QNet->Toutside);
}

void TElectron::init_Tref() {
    Tref = &(QNet->T0);
}

void TElectron::update_dT(){
    if (!dTvar.empty())
        dT = dTvar[QNet->timestep];
    }

void TElectron::compute_dP(){
    if (DN > 0.0 && length > 0.0) {
        Pipe.PressureDropCalculation();
        dP_m = Pipe.DeltaPr_m;
        dP = 2 * Pipe.DeltaPr;          // two pipes !!
        speed = Pipe.Speed;
    }
}

void TElectron::compute_tau(){  // currently only the tau of pipes
    if (DN > 0.0 && length > 0.0) {
        Pipe.tauCalculation();
        tau = Pipe.tau;
    }
}


void TElectron::compute_T_supply() {
    // does nothing
}

void TElectron::compute_T_return() {
    if (flow->Cp == 0){
        cerr << "Cp is zero!" << endl;
        exit(1);
    }
    if (tau == 0.0) {
        flow->TOut(flow->TIn() - dT);
        dPower = 0.0;
    }
    else {
        flow->TOut(*Tambient + ((flow->TIn() - *Tambient) * exp(-tau) - dT) * exp(-tau));
        dPower = flow->massflow() * flow->Cp * (flow->TIn() - flow->TOut() - dT);
    }
}

void TElectron::compute_material_cost(){
    // connecting pipes
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

    // heat exchanger (Jeremy Rolle page 32)
    double power = D.power_max/1000;    // in kW TODO (what/where is the definition of power of SST?)
    material_cost += (5.56e-10 * pow(power, 3) - 6.73e-06 * pow(power, 2) + 0.0262 * power + 2.617) * 1000;
}

void TElectron::compute_engineering_cost(const char* region) {
    // connecting pipes
    const char* fluidname = flow->fluid->FluidName;
    if (strcmp(fluidname, "H2O") == 0) {
        // master thesis Briguet page 24
        double excavation;
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

    // heat exchanger (Jeremy Rolle page 32): substract blue points from orange points
    double power = D.power_max/1000;    // in kW
    engineering_cost += (2.35E-10 * pow(power, 3) + 2.26E-07 * pow(power, 2) + 6.36E-03 * power + 6.63) * 1000;
}

void TElectron::compute_operating_cost(){
    // maintenance cost for SST (Jeremy Rolle page 34)
    operating_cost = 0.12915 * pow(D.power_max / 1000, -0.813) * D.power_max;
}

void TElectron::compute_resources_cost(){
    if (QNet->timestep == 0)
        total_resources_cost = 0.0;
    // all time-steps
    resources_cost = QNet->res_cost_per_kWh * (flow->Power/1000);
    total_resources_cost += resources_cost;
}

void TElectron::compute_power_from_district(){
    double factor = 1;        //TODO heat loss (this is NOT efficiency of heat exchanger!)

    // DeltaQ = lossfactor * (Tin - Text);   lossfactor in [W/K]
    // flow->Power = demand + DeltaQ;

    if(D.power.empty()){
        cerr << "District not defined/computed" << endl;
        exit(1);
    }
    demand = D.power[QNet->timestep];         // what is needed by district
    flow->Power = factor * demand;            // what is yielded by SST
}

void TElectron::compute_mu_from_power(){
    compute_massflow_from_power();
    compute_mu_from_massflow();
}

