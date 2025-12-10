//
// Created by Lucile Schulthess on 16.08.2021.
//

#include <algorithm>
#include "District.h"

using namespace std;

void District::iterate_list(const char* ifile_in, unsigned int inbTitleLines, unsigned int inbTitleColumns,char isep,bool iren) {
    TInputFile InputFile((char *) ifile_in, inbTitleLines, inbTitleColumns, isep);
    InputFile.open();
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    Records = *InputFile.GetRecords();
    Building Building;
    BuildingNorm Building_list; //default constructor does nothing!

    for (int i=1; i <= InputFile.nRecords(); ++i) {
        Building.EGID = ((Records.Get(i-1))->Title)->liste();
        Building.type_heat = "HeatECS";  //TODO where to find this information?

        Building_list.List(Records, i);   // get i-th building from InputFile and inits the norm
        Qh.push_back(Building_list.Qh);     //TODO do we need these vectors?
        Qc.push_back(Building_list.Qc);
        Qv.push_back(Building_list.Qv);
        Qe.push_back(Building_list.Qe);
        Qw.push_back(Building_list.Qw);
        Mw.push_back(Building_list.Mw);
        type.push_back(Building_list.type);
        Thref.push_back(Building_list.Thref);
        Teref.push_back(Building_list.Teref);
        Tdem.push_back(Building_list.Tdem);
        heq.push_back(Building_list.heq);

        Building.init_energy_from_norm(Building_list);
        B.push_back(Building);
    }
}

void District::init_buildings(vector<double> &energy, vector<double> &Tref, vector<int> &types, vector<const char*> &type_heat){
    if (energy.size() != Tref.size() || energy.size() != type_heat.size() || energy.size() != types.size()){
        cerr << "Input dimensions do not match!" << endl;
        exit(1);
    }
    Building Building;
    for (int i=0; i<energy.size(); i++){
        Building.energy = energy[i];
        if (Tref[i] < 100)  // in Celsius
            Building.Tref = Tref[i] + 273.15;
        else                // in Kelvin
            Building.Tref = Tref[i];
        Building.type = types[i];
        Building.type_heat = type_heat[i];
        B.push_back(Building);
    }
}

void District::init_buildings(vector<double> &energy, vector<double> &Tref, vector<double> &Tneed, vector<int> &types, vector<const char*> &type_heat){
    if (energy.size() != Tref.size() || energy.size() != type_heat.size() || energy.size() != types.size()){
        cerr << "Input dimensions do not match!" << endl;
        exit(1);
    }
    Building Building;
    for (int i=0; i<energy.size(); i++){
        Building.energy = energy[i];
        if (Tref[i] < 100)  // in Celsius
            Building.Tref = Tref[i] + 273.15;
        else                // in Kelvin
            Building.Tref = Tref[i];
        if (Tneed[i] < 100)  // in Celsius
            Building.Tdem = Tneed[i] + 273.15;
        else                // in Kelvin
            Building.Tdem = Tneed[i];
        Building.type = types[i];
        Building.type_heat = type_heat[i];
        B.push_back(Building);
    }
}


void District::calcul_power(environment &ext){  // sum of all buildings B
    power=vector<double>(ext.Time, 0);
    for(auto & building : B) {
        building.calcul_power(ext);
        power = sum(building.power, power);
    }
    Q_total = sum(power);   //TODO round-off errors in power: Q_total != SST-Energy
    power_max = max(power);
}

/*
vector<double> District::HeatinginTime(environment &Tempex){
   vector<double> P(Tempex.Time);
   vector<double> DDH(Tempex.Time);
   double corr=1.0;
   double corr_factor=Tempex.Time/(365.0*24.0);
   if (corr_factor == round(corr_factor) )  // only if full year
       corr+=(9.4+273.15-mean(Tempex.T))*0.06;
   corr *= corr_factor;
   for (int j(0);j<Qh.size();++j) {
       DDH = HourlyHeatingPower(Tempex,Thref[j]);
       P = sum(multi(DDH,Qh[j]*corr), P);
   }
   return P;
}

vector<double> District::CoolinginTime(environment &Tempex){
    vector<double> P(Tempex.Time);
    vector<double> DDC(Tempex.Time);
    double corr = 1.0;
    double corr_factor=Tempex.Time/(365.0*24.0);

    for (int j(0);j<Qc.size();++j) {
        DDC = HourlyCoolingPower(Tempex, Teref[j]+1);
        if (corr_factor == round(corr_factor) ) {  // only if full year
            int countdays = 0;
            for (int i = 0; i < Tempex.Time; i++) {     // count number of days with cooling needs
                if (DDC[i] != 0.0) {
                    countdays++;
                    i = (i / 24) * 24 + 23;     // jump to next day
                }
            }
            // Nombre de cooling degree day à météo de référence
            double DDC_ref = -2.2373 * pow(Teref[j]+1-273.15,3) + 188.59 * pow(Teref[j]+1-273.15,2) - 5309.2 * (Teref[j]+1-273.15)+49941;
            corr = countdays / DDC_ref;
        }
        corr *= corr_factor;
        P = sum(multi(DDC,Qc[j]*corr), P);
    }
    return P;
}


vector<double> District::ECSinTime(vector <double> &Q){
    vector<double> P(ECS.hour);
    for (int j(0);j<Q.size();++j) {
        for (int i(0);i<ECS.hour;++i){
            P[i]+=ECS.profil_ECS[i][type[j]]*Q[j]/365.0;      //sum of all P[i]
        }
    }
    return P;
}
*/

void District::Display(){
    cout << endl;
    cout <<"/****************************************************************************/" << endl;
    cout<<"                         District                           "<<endl<<endl;
    cout<<"Heat energy total              [kwh] :       :"<<sum(Qh)<<endl;
    cout<<"ECS energy total               [kwh] :       :"<<sum(Qw)<<endl;
    cout<<"Cooling energy total           [kwh] :       :"<<sum(Qc)<<endl;
    cout<<"Ventil energy total (elec)     [kwh] :       :"<<sum(Qv)<<endl;
    cout<<"Electric energy total          [kwh] :       :"<<sum(Qe)<<endl;
    cout<<"Liter of water total           [L] :         :"<<sum(Mw)<<endl;
    cout<<endl;
    cout <<"/****************************************************************************/" << endl;
    cout << endl;
}

/*
vector<double> District::HourlyHeatingPower(environment &Tempex, double Tref){
    vector<double> DH(Tempex.Time);
    for(int i(0); i<Tempex.Time; ++i)
        DH[i] = std::max(Tref-Tempex.T[i], 0.0);
    return vN(DH);
}


vector<double> District::HourlyCoolingPower(environment &Tempex, double Tref){
    vector<double> DC(Tempex.Time);
    for(int i(0); i<Tempex.Time; ++i)
        DC[i] = std::min(Tref-Tempex.T[i], 0.0);
    return vN(DC);
}
*/
/****************************************************************************/




