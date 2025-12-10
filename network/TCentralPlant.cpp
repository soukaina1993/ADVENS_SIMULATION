//
// Created by cornelia.blanke on 06.03.2023.
//

#include "TCentralPlant.h"
#include <iostream>

using namespace std;

Boiler::Boiler(Boiler::Type type, double max_power) : type(type), max_power(max_power) {
    if (PRICE_PER_kWh.count(type) == 1)
        res_cost_per_kWh = PRICE_PER_kWh[type];
}
Boiler::Boiler(Boiler::Type type, double max_power, double price): type(type), max_power(max_power),
                                                           res_cost_per_kWh(price) {}


TCentralPlant::TCentralPlant()=default;

void TCentralPlant::display(){
    cout << "\n";
    cout <<"/****************************************************************************/\n";
    cout<<"TCentralPlant                         :  \n\n";
    cout<<"Maximum Power                         :       " << max_power/1e6 << " [MW]\n";
    cout<<"Maximum Mass Flow                     :       " << max_massflow  << " [kg/s]\n";
    cout<< "\n";
    cout <<"/****************************************************************************/\n";
    cout << endl;
}

void TCentralPlant::compute_m_from_table() {
    if (QNet->mtable.empty()){
        cerr << "mtable not defined" << endl;
        exit(1);
    }
    else
        m = QNet->Z;
}
void TCentralPlant::compute_m() {
    m = QNet->Z;
}

void TCentralPlant::compute_mu_from_table(){
    if (QNet->muTable.empty()){
        cerr << "muTable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int) (QNet->nbBranches / 2);
        mu = QNet->muTable[0][branch0];
    }
}

void TCentralPlant::compute_massflow_from_mu(){
    double Mdot = massflow_from_mu(mu, max_massflow, QNet->Z);
    flow->massflow(Mdot);
}

void TCentralPlant::compute_mu_from_massflow(){
    mu = mu_from_massflow(flow->massflow(), max_massflow, QNet->Z);
}

void TCentralPlant::modify_pricelist(vector<Boiler::Type>& types, vector<double>& prices) {
    if (types.size() != prices.size()) {
        cerr << "Mismatch of dimensions" << endl;
        exit(1);
    }
    for (int i=0; i<types.size(); i++)
        PRICE_PER_kWh[types[i]] = prices[i];
}

void TCentralPlant::modify_pricelist(unordered_map<Boiler::Type, double>& prices) {
    for (auto& price : prices)
        PRICE_PER_kWh[price.first] = price.second;
}

void TCentralPlant::init_boilers(vector<Boiler::Type>& types, vector<double>& max_powers, vector<double>& prices) {
    if (types.size() != max_powers.size() || types.size() != prices.size()) {
        cerr << "Mismatch of dimensions" << endl;
        exit(1);
    }
    boilers.resize(types.size());
    for (int i=0; i<types.size(); i++) {
        Boiler boiler(types[i], max_powers[i]);
        if (prices[i] > 0.0)    // override default value
            boiler.res_cost_per_kWh = prices[i];
        boilers[i] = boiler;
    }
}

void TCentralPlant::init_boilers(vector<Boiler::Type>& types, vector<double>& max_powers) {
    vector<double> prices (types.size(),-1.0);
    init_boilers(types, max_powers, prices);
}

void TCentralPlant::init_boilers(){
    stringstream boilertypes;
    boilertypes << QNet->get_setting("Boiler Type");
    stringstream boilerpowers;
    boilerpowers << QNet->get_setting("Boiler Power");
    stringstream boilerprices;
    boilerprices << QNet->get_setting("Boiler Price/kWh");
    string str;
    vector<Boiler::Type> types;
    vector<double> max_powers, prices;
    int i;
    while (getline(boilertypes, str, '/')){
        const vector<const char*> boilertypes = {"OIL", "GAS", "WOOD", "PELLET", "NONE"};
        for(i=0; i < boilertypes.size(); i++)
            if (strcmp(str.c_str(), boilertypes[i]) == 0)
                break;
        if (i == boilertypes.size() + 1){
            cerr << "Boiler type "<< str << " not available" << endl;
            exit(1);
        }
        else
            types.push_back((Boiler::Type) i);
    }
    while (getline(boilerpowers, str, '/'))
        max_powers.push_back(stod(str));

    while (getline(boilerprices, str, '/')) {
        try {
            prices.push_back(stod(str));
        }
        catch(const std::invalid_argument&) {
            prices.push_back(-1.0);
        }
    }
    init_boilers(types, max_powers, prices);
}

void TCentralPlant::add_boiler(Boiler::Type typ, double maxpower) {
    Boiler boiler(typ, maxpower);
    boilers.push_back(boiler);
}

void TCentralPlant::compute_maxvalues(){
    if (max_massflow == 0.0 && max_power == 0.0){
        cerr << "Specify at least max_massflow or max_power" << endl;
        exit(1);
    }
    else if (max_massflow == 0.0)
        max_massflow = max_power / (flow->Cp * dT);
    else if (max_power == 0.0)
        max_power = max_massflow * flow->Cp * dT;
}

void TCentralPlant::compute_maxvalues_from_boilers(){
    double sum_power=0.0;
    for (auto & boiler : boilers)
        sum_power += boiler.max_power;

    if (max_power==0.0)                     // nothing specified
        max_power = sum_power;
    else if (max_power > sum_power) {       // missing boiler definition --> add boiler
        cout << "Adding a boiler of type NONE" << "\n";
        add_boiler(Boiler::NONE, max_power - sum_power);
    }
    else if (max_power < sum_power)
        cerr << "Warning: Max Power of Central Station is smaller than sum of all boilers." << endl;
    compute_maxvalues();
}

void TCentralPlant::init_pressure(const double Pin, const double Pout){
    flow->PIn (Pin);
    flow->POut (Pout);
    dP = Pin - Pout;
}

void TCentralPlant::init_temperature(double Tin, double Tout){
    if (Tin < 273)
        Tin += 273.15;
    if (Tout < 273)
        Tout += 273.15;
    flow->TIn (Tin);
    flow->TOut (Tout);
    dT = Tin - Tout;
}

void TCentralPlant::connect_network() {
    if (QNet->Table_of_Segments.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    next.push_back(QNet->Table_of_Segments[0][branch0]);
}

void TCentralPlant::connect_flow() {
    if (QNet->Table_of_Segments.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    int branch0 = (int)(QNet->nbBranches/2);
    nextflow.push_back(QNet->Table_of_Segments[0][branch0]->flow);
}

void TCentralPlant::init_Tambient(const bool below) {
    if (below)
        Tambient = &(QNet->Tground);
    else
        Tambient = &(QNet->Toutside);
}

void TCentralPlant::init_Tref() {
    Tref = &(QNet->T0);
}

void TCentralPlant::compute_tau(){}

void TCentralPlant::compute_material_cost(){
    // boilers (Jeremy Rolle page 25/26)
    for (auto & boiler : boilers){
        double power = boiler.max_power/1000;       //power in kW
        if (boiler.type == Boiler::OIL || boiler.type == Boiler::GAS) {
            material_cost = (29.324 * ln(power) - 134.1) * 1000;
            // add a chimney (Jeremy Rolle page 28)
            material_cost += (779.53 * ln(power) - 3597.8) * 1000;
        }
        else if (boiler.type == Boiler::WOOD) {
            material_cost = (-3.77e-8 * pow(power, 3) + 0.000303 * pow(power, 2) -
                             0.0933 * power + 322.54) * 1000;
            material_cost += (757.97 * ln(power) - 3230.1) * 1000;
        }
        else if (boiler.type == Boiler::PELLET) {
            material_cost = (70.327 * ln(power) - 206.12) * 1000;
            material_cost += (477.54 * ln(power) - 1555.5) * 1000;
        }
        else if (boiler.type == Boiler::NONE)
            material_cost = 0.0;
        else {
            cerr << "Boiler type not implemented" << endl;
            exit(1);
        }
        // for each boiler, add an electrofilter (Jeremy Rolle page 27)
        if (boiler.type != Boiler::NONE)
            material_cost += -(5.31e-6 * pow(power, 2) + 0.1015 * power + 11.745) * 1000;
    }

    // pumps (Jeremy Rolle page 32)      TODO needs to be evaluated after simulation from the results (in Python?)
    /*
    if (flow->pStateIn->Ro == 0){
        cerr << "Init flow first" << endl;
        exit(1);
    }
    if (max_massflow == 0)
        cerr << "Warning: Max massflow not defined" << endl;

    double V = max_massflow / flow->pStateIn->Ro / 3600;    // in m^3/h
    double P = flow->pStateIn->GetP()/1.0e5;                // in bar
    material_cost += 1051 - 304.9 * P + 9.934 * V + 72.66 * pow(P, 2) + 5.996 * P * V + 0.01878 * pow(V, 2);
     */
}

void TCentralPlant::compute_operating_cost(){
    // maintenance cost for the boilers (Jeremy Rolle page 33)
    operating_cost = 0.0;
    for (auto & boiler : boilers) {
        if (boiler.type == Boiler::OIL || boiler.type == Boiler::GAS)
            operating_cost += 1.355 * pow(boiler.max_power / 1000, -0.729) * boiler.max_power;
        else if (boiler.type == Boiler::WOOD)
            operating_cost += 0.82887 * pow(boiler.max_power / 1000, -0.579) * boiler.max_power;
        else if (boiler.type == Boiler::PELLET)
            operating_cost += 5.9796 * pow(boiler.max_power / 1000, -0.943) * boiler.max_power;
        else if (boiler.type == Boiler::NONE)
            operating_cost += 0.0;
        else {
            cerr << "Boiler type not implemented" << endl;
            exit(1);
        }
    }
}

void TCentralPlant::compute_res_cost_per_kWh(){
    double sum_power = 0.0;
    for (auto & boiler : boilers)
        sum_power += boiler.max_power;
    if (sum_power == 0.0)
        cerr << "Operating costs cannot be computed: Division by zero!" << endl;
    else {
        res_cost_per_kWh = 0.0;
        for (auto &boiler: boilers) {      // weighted average
            //TODO find out at which percentage each boiler is used
            double percentage = boiler.max_power / sum_power;
            res_cost_per_kWh += percentage * boiler.res_cost_per_kWh;
        }
    }
}

void TCentralPlant::compute_resources_cost(){   // evaluate after each time step!
    if (QNet->timestep == 0) {
        total_resources_cost = 0.0;
        compute_res_cost_per_kWh();
    }
    // all time-steps
    resources_cost = res_cost_per_kWh * (flow->Power/1000);
    total_resources_cost += resources_cost;
}