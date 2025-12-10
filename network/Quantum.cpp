//
// Created by cornelia.blanke on 06.03.2023.
//

#include "Quantum.h"


Quantum::Quantum() {
    init();
}
Quantum::~Quantum() {
    delete flow;
}

void Quantum::init (){}

BasicES & Quantum::performance(){}
float Quantum::objective(TSequence& seq){return 0.0;}
void Quantum::display(){}

void Quantum::init_fluid(FluidCore* &fluidin1){
    char* fluidname= (char*)fluidin1;
    fluid_define(flow->fluid,fluidname);
    flow->init(flow->fluid);
    double Cv;  // not used
    flow->fluid->CvCp_vT(flow->vIn(), flow->TIn(),Cv,flow->Cp);
};

void Quantum::init_fluid(string fluidin1){
    fluid_define(flow->fluid, fluidin1);
    flow->init(flow->fluid);
    double Cv;  // not used
    flow->fluid->CvCp_vT((flow->vIn()+flow->vOut())/2, (flow->TIn()+flow->TOut())/2,Cv,flow->Cp);
}

void Quantum::init_Tsource(Quantum* Source){
    Tsource = &(Source->flow->pStateIn->T);
}

void Quantum::compute_mubar_from_mu(){
    if (m==0) {
        cerr << "m not defined" << endl;
        exit(1);
    }
    else
        mu_bar = mu/m;
}

double Quantum::massflow_from_mu(double mu, double max_source_massflow, unsigned int Z){
    if (Z==0){
        cerr << "Z not defined" << endl;
        exit(1);
    }
    return mu * max_source_massflow/Z;
}
double Quantum::mu_from_massflow(double Mdot, double max_source_massflow, unsigned int Z){
    if (max_source_massflow==0){
        cerr << "Max Massflow not defined" << endl;
        exit(1);
    }
    return Mdot * Z/max_source_massflow;
}

void Quantum::set_nextTIn(double TIn){
    if (next.empty()){      // okay if TElectron, else error
        cerr << "next-pointer not initialised" << endl;
        exit(1);
    }
    for (int i=0; i<next.size(); i++)
            next[i]->flow->TIn(TIn);
}

double Quantum::get_nextTIn(){
    if (next.empty()){      // okay if TElectron, else error
        cerr << "next-pointer not initialised" << endl;
        exit(1);
    }
    return next[0]->flow->TIn();
}

double Quantum::get_nextTOut(){
    if (next.empty()){
        cerr << "next-pointer not initialised" << endl;
        exit(1);
    }
    if (next.size()==1)
        return next[0]->flow->TOut();
    else if (flow->massflow() != 0.0) {
        double sum=0.0;
        for (int i = 0; i < next.size(); i++)
            sum += next[i]->flow->TOut() * next[i]->flow->massflow();
        return sum/flow->massflow();
    }
    else {
        cout << "Warning: No massflow" << endl;
        double sum=0.0;
        for (int i = 0; i < next.size(); i++)
            sum += next[i]->flow->TOut();
        return sum/next.size();
    }
}

void Quantum::compute_dT_from_delta(){
    dT = delta * (*Tsource - *Tref);
}
void Quantum::compute_delta_from_dT(){
    if (*Tsource == *Tref){
        cerr << "Division by zero: Tsource == Tref" << endl;
        exit(1);
    }
    delta = dT / (*Tsource - *Tref);
}

void Quantum::compute_nubar() {
    if (mu_bar == 0.0)
        compute_mubar_from_mu();
    nu_bar = mu_bar * (1 - exp(-2*tau)*(1-delta)) * (flow->TIn() - *Tambient);
}

void Quantum::compute_nubar_from_T(double k0) {
    if (m*k0 == 0.0){
        cerr << "Division by zero" << endl;
        exit(1);
    }
    nu_bar = flow->massflow() * flow->Cp * (flow->TIn() - flow->TOut()) / (m * k0);
}

void Quantum::compute_epsilon() {
    if (tau == 0.0)
        epsilon = 1.0;
    else {
        double denominator = 1 - exp(-2 * tau) * (1 - delta);
        if (denominator == 0) {
            cerr << "Division by zero" << endl;
            exit(1);
        }
        epsilon = (delta * exp(-tau)) / denominator;
    }
}

void Quantum::compute_tau() {}

void Quantum::compute_dP() {}

void Quantum::estimate_powerloss() {
    if (flow->Cp == 0){
        cerr << "Cp is zero!" << endl;
        exit(1);
    }
    if (flow->massflow() == 0)
        dPower = 0.0;
    else {
        double factor = exp(-tau);
        dPower = flow->massflow() * flow->Cp *
                 ((flow->TIn() - *Tambient) * (1-pow(factor,2)) - dT*(1-factor));
    }
}


void Quantum::compute_T_supply() {
    if (tau == 0.0)
        set_nextTIn(flow->TIn());
    else
        set_nextTIn(*Tambient + (flow->TIn() - *Tambient) * exp(-tau) );
}

void Quantum::compute_T_return() {
    dT = get_nextTIn() - get_nextTOut();
    compute_delta_from_dT();
    if (tau == 0.0) {
        flow->TOut(get_nextTOut());
        dPower = 0.0;
    }
    else {
        flow->TOut(*Tambient + (get_nextTOut() - *Tambient) * exp(-tau));
        dPower = flow->massflow() * flow->Cp * (flow->TIn() - flow->TOut() - dT);   //sum of supply + return
    }
}

void Quantum::compute_massflow_from_power(){
    if(dT==0){
        cerr << "dT = 0" << endl;
        exit(1);
    }
    double Mdot = flow->Power / (flow->Cp * dT);
    flow->massflow(Mdot);
}

void Quantum::compute_power_from_massflow(){
    flow->Power = flow->massflow() * flow->Cp * dT;
}

void Quantum::compute_material_cost(){}
void Quantum::compute_engineering_cost(const char* region){}
void Quantum::compute_operating_cost(){}
void Quantum::compute_resources_cost(){}

void Quantum::compute_total_invest() {
    total_invest = material_cost + engineering_cost;
}