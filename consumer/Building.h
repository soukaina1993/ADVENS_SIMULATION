//
// Created by lucile.schulthess on 13.08.2021.
//

#pragma once

#include <iostream>
#include <unordered_map>
#include "Norm_Building.h"
//#include "sst.h"
#include "../environment/environment.h"
#include <cstring>
#include <cstdlib>


//const vector<string> type_heatmode =  {"Heat","HeatECS","ECS","Cool"};

//const vector<string> type_building = { "Habitat_collectif", "Habitat_individuel",
//                       "Administration", "Ecole","Commerce","Restauration","Lieu_de_rassemblement",
//                       "Hopital","Industrie","Depot","Installation_sportive","Piscine_couverte" };
//static unordered_map<string, int> type_building =
//        {{"Habitat_collectif", 0}, {"Habitat_individuel", 1},
//         {"Administration", 2}, {"Ecole", 3}, {"Commerce", 4}, {"Restauration", 5},
//         {"Lieu_de_rassemblement", 6}, {"Hopital", 7}, {"Industrie", 8},
//         {"Depot", 9}, {"Installation_sportive", 10}, {"Piscine_couverte", 11} };



//vector<double> partECS2={0.25,0.25,0.15,0.15,0.12,0.56,0.044,0.44,0.27,0.014,0.47,0.92};
//const vector<double> partECS_norm={0.17,0.077,0.039,0.033,0.12,0.56,0.044,0.44,0.27,0.014,0.47,0.92}; // faux
const vector<double> partECS_norm={0.147,0.072,0.038,0.032,0.029,0.360,0.042,0.304,0.026,0.014,0.317,0.480};   // sigma

const norm_ecs ECS;

/****************************************************************************/

//class Building:public sst {
class Building {        //TODO derive from?
public:
    Building();
    void Display();

// pour l'héritage sst
//    Building(vector <string> &data_input);

    void calcul_power(environment &ext);

    void define_type(const char* type_name);
    void define_heatmode(const char* heatmode);

    static int get_type(const char* type_name);         // static function, to be used outside of Building class
    static string get_type(int type_int);

    void init_norm(const char* ifile_in, const char* egid);
    void init_energy_from_norm(BuildingNorm &Bref);
    void init_energy_from_data(double ienergy, double iTref, double iTdem);

    vector<double> HourlyCoolingPower(environment &Tempex,double Teref);
    vector<double> HourlyHeatingPower(environment &Tempex,double Thref);
    vector<double> HeatingPower(environment &Tempex, bool correction=false);
    vector<double> CoolingPower(environment &Tempex, bool correction=false);
    vector<double> ECSinTime(double &Q);

    double Tref=273.15;   // Tref in Kelvin, when Heating or Cooling is switched on
    double Tdem=273.15;   // Tdem in Kelvin, minimum temperature needed for radiators
    int type=1;
    string EGID, type_heat;     // default value set in constructor type_heat="HeatECS"
    vector<double>power;
    double energy=-1.0, power_max=-1.0;       // energy [kWh] in one year!
};