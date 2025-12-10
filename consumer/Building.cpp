//
// Created by Lucile Schulthess on 16.08.2021.
// Modified by Cornelia Blanke on 09.05.2023.
//

#include <algorithm>    // std::max, std::min
#include "Building.h"


using namespace std;

Building::Building(){
    type_heat = "HeatECS";
}

vector<double> Building::CoolingPower(environment &Tempex, bool correction){
    vector<double> P(Tempex.Time);
    vector<double> DDC=HourlyCoolingPower(Tempex,Tref);
    double corr=1.0;
    double corr_factor=Tempex.Time/8760.0;

    if (correction && (corr_factor == round(corr_factor)) ) {  // only if full year
        int countdays = 0;
        for (int i = 0; i < Tempex.Time; i++) {     // count number of days with cooling needs
            if (DDC[i] != 0.0) {
                countdays++;
                i = (i / 24) * 24 + 23;     // jump to next day
            }
        }
        // Nombre de cooling degree day à météo de référence
        double DDC_ref = -2.2373 * pow(Tref+1-273.15,3) + 188.59 * pow(Tref+1-273.15,2) - 5309.2 * (Tref+1-273.15)+49941;
        corr = countdays / DDC_ref;
    }
    corr *= corr_factor;

    for(int i(0); i<Tempex.Time; ++i)
        P[i] = DDC[i]*energy*corr;
    return P;
}

vector<double> Building::HeatingPower(environment &Tempex, bool correction){
    vector<double> P(Tempex.Time);
    vector<double> DDH=HourlyHeatingPower(Tempex,Tref);
    double corr=1.0;
    double corr_factor=Tempex.Time/8760.0;

    if (correction && (corr_factor == round(corr_factor)) )  // only if full year
        corr+=(9.4+273.15-mean(Tempex.T))*0.06;
    corr *= corr_factor;

    for(int i(0); i<Tempex.Time; ++i)
        P[i] = DDH[i]*energy*corr;
    return P;
}


vector<double> Building::HourlyHeatingPower(environment &Tempex,double Thref){
    vector<double> DH(Tempex.Time);
    for(int i(0); i<Tempex.Time; ++i){
        DH[i] = std::max(Thref-Tempex.T[i], 0.0);}
    return vN(DH);
}

vector<double> Building::HourlyCoolingPower(environment &Tempex,double Teref){
    vector<double> DC(Tempex.Time);
    for(int i(0); i<Tempex.Time; ++i){
        DC[i] = std::min(Teref-Tempex.T[i], 0.0);}
    return vN(DC);
}


vector<double> Building::ECSinTime(double &Q){
    vector<double> P(ECS.hour);
    for (int i(0);i<ECS.hour;++i){
        P[i] = ECS.profil_ECS[i][type]*Q/365;
    }
    return P;
}

void Building::Display(){
    // find building type from int type
    string buildingtype;
    buildingtype = get_type(type);

    cout << endl;
    cout <<"/****************************************************************************/" << endl;
    cout<<"Building EGID                         :       "<<EGID<<endl<<endl;
    cout << "Building Type                         :       " << buildingtype << endl;
    cout<<"Demand energy total             [kWh] :       "<<energy<<endl;
    cout<<"Switch-on temperature          [degC] :       "<<Tref-273.15<<endl;
    cout<<"Demand temperature             [degC] :       "<<Tdem-273.15<<endl;
    cout<<endl;
    cout <<"/****************************************************************************/" << endl;
    cout << endl;
}
//// pour l'héritage sst
/*Building::Building(vector <string> &data_input){

    Tref=stod(data_input[3])+273.15;

    sst_type=data_input[4];
    dT=stod(data_input[5]);
    type=define_type(data_input[8]);
    type_heat=data_input[9];
    Tmin=stod(data_input[6]);
    Tpinch=stod(data_input[7]);
    L_racc=stod(data_input[10]);
    DN_racc=stod(data_input[11]);
    alt=stod(data_input[12]);

    if(data_input[1]=="Energy"){
        energy= stod(data_input[2]);
        if (type_heat=="Cool" && energy>0)
            energy=-energy;}
    else if(data_input[1]=="Power")
        power_max= stod(data_input[2]);
    else if(data_input[1]=="EGID"){
        int n = data_input[2].length();
        char egid[n+1];
        strcpy(egid,data_input[2].c_str());
        init_norm((const char *) egid);
    }
};*/
static bool calcul_power_output=true;   // print the warning only once
void Building::calcul_power(environment &ext){
    double corr_factor=ext.Time/8760.0;
    if (corr_factor != round(corr_factor) && calcul_power_output ) {
        cout << "Attention: calculation is not for a full year!" << endl;
        calcul_power_output = false;
    }
    if(energy==0.0 || power_max==0.0)   // no energy demand
        power = vector<double>(ext.T.size(), 0);

    else if(type_heat=="Heat")
        power=HeatingPower(ext);

    else if(type_heat=="HeatECS") {
        vector<double>DDH=HourlyHeatingPower(ext,Tref);
        double Qwater=energy*partECS_norm[type];        //TODO add correction for colder/hotter years

        power= multi(DDH,(energy-Qwater)*corr_factor);  // heating
        //norm_ecs ECS;
        vector<double>power_ECS=ECSinTime(Qwater);

        for(int t=0; t<ext.T.size(); t++)        // add water
            power[t]+=(power_ECS[t%24]);
    }
    else if(type_heat=="ECS"){
        //norm_ecs ECS;
        vector<double>power_ECS=ECSinTime(energy);  // 24 hours
        power = vector<double>(ext.T.size(), 0);
        for(int t=0; t<ext.T.size(); t++)
            power[t] = power_ECS[t%24];
    }
    else if(type_heat=="Cool"){
        power=CoolingPower(ext);
    }
    else {
        cerr << "Heat mode not defined" << endl;
        exit(1);
    }

    if(energy==-1.0 && power_max!=-1.0){        // energy not defined, power_max defined
        double corr=power_max/min(power);   // power is negative
        power=multi(power,corr);
        energy=sum(power);}
    else power_max=max(power);
}

void Building::define_type(const char* type_name){
    type = get_type(type_name);
}

void Building::define_heatmode(const char* heatmode){
    const vector<const char*> type_heatmode =  {"Heat", "HeatECS" ,"ECS", "Cool"};
    type_heat = "HeatECS";    // default
    for(const auto& mode : type_heatmode) {
        if (strcmp(heatmode, mode) == 0)
            type_heat = mode;
    }
}

int Building::get_type(const char* type_name){
    int i=0;
    string s;
    while (s != "Unknown"){
        s = get_type(i);
        if (strcmp(s.c_str(), type_name)==0)
            return i;
        i++;
    }
    return 1;       // default
}

string Building::get_type(const int type_int){
    const vector<const char*> type_building = { "Habitat_collectif", "Habitat_individuel",
                                          "Administration", "Ecole","Commerce","Restauration","Lieu_de_rassemblement",
                                          "Hopital","Industrie","Depot","Installation_sportive","Piscine_couverte" };
    string str;
    if (type_int >= 0 && type_int < type_building.size())
        str = type_building[type_int];
    else
        str = "Unknown";
    return str;
}

void Building::init_norm(const char* ifile_in, const char* egid){
    EGID = egid;
    BuildingNorm Bref(ifile_in, egid);
    init_energy_from_norm(Bref);
}

void Building::init_energy_from_norm(BuildingNorm &Bref){
    type = Bref.type;
    Tdem = Bref.Tdem;
    if (type_heat=="Heat"){
        energy=Bref.Qh;
        Tref=Bref.Thref;}
    else if (type_heat=="HeatECS"){
        energy=Bref.Qh+Bref.Qw;
        Tref=Bref.Thref;}
    else if (type_heat=="ECS"){
        energy=Bref.Qw;
        Tref=Bref.Thref;}
    else if (type_heat=="Cool"){
        energy=Bref.Qc;
        Tref=Bref.Teref;}
}

void Building::init_energy_from_data(double ienergy, double iTref, double iTdem){
    energy=ienergy;
    Tref=iTref;
    Tdem = iTdem;
}