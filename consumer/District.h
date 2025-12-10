//
// Created by lucile.schulthess on 13.08.2021.
//


#pragma once
#include <utility>

#include "Building.h"


/****************************************************************************/


/****************************************************************************/

//class District: public sst{
class District{
public:
    District(){};
    ~District(){};

    void calcul_power(environment &ext);
    void calcul_flow(double Tinlet){};


    District(const char* ifile_in, unsigned int nTline=1, unsigned int nTcol=1, char isep = '\t',bool iren=false){
        file_in=(char*)ifile_in;
        sep = isep;
        nbTitleLines = nTline;
        nbTitleColumns = nTcol;
        ren=iren;
        init();
    };

    District(bool iren, const char* ifile_in, unsigned int nTline=1, unsigned int nTcol=1, char isep = '\t'){
        ren=iren;
        file_in=(char*)ifile_in;
        sep = isep;
        nbTitleLines = nTline;
        nbTitleColumns = nTcol;
        init();
    };

    void iterate_list(const char* ifile_in, unsigned int inbTitleLines = 1,unsigned int inbTitleColumns = 1,char isep = '\t',bool iren=false);

    void init(){
        iterate_list(file_in, nbTitleLines, nbTitleColumns, sep, ren);
        Tneed_max = max(Tdem);
    };

    // without Tneed
    District(vector<double> &energy, vector<double> &Tref, vector<int> &type,
             vector<const char*> &type_heat, bool iren=false){  // List of Buildings
        ren = iren;
        init_buildings(energy, Tref, type, type_heat);
    };

    // without Tneed
    District(vector<double> &energy, vector<double> &Tref, vector<const char*> &typestring,
             vector<const char*> &type_heat, bool iren=false){  // List of Buildings
        ren = iren;
        type.reserve(typestring.size());
        for (int i=0; i<typestring.size(); i++)
            type[i] = Building::get_type(typestring[i]);

        init_buildings(energy, Tref, type, type_heat);
    };

    // with Tneed
    District(vector<double> &energy, vector<double> &Tref, vector<double> &Tneed, vector<int> &type,
             vector<const char*> &type_heat, bool iren=false){  // List of Buildings
        ren = iren;
        for (double & temp : Tneed) {
            if (temp < 100)
                temp += 273.15;
        }
        init_buildings(energy, Tref, Tneed, type, type_heat);
        Tneed_max = max(Tneed);
    };

    // with Tneed
    District(vector<double> &energy, vector<double> &Tref, vector<double> &Tneed, vector<const char*> &typestring,
             vector<const char*> &type_heat, bool iren=false){  // List of Buildings
        ren = iren;
        type.reserve(typestring.size());
        for (int i=0; i<typestring.size(); i++)
            type[i] = Building::get_type(typestring[i]);
        for (double & temp : Tneed) {
            if (temp < 100)
                temp += 273.15;
        }
        init_buildings(energy, Tref, Tneed, type, type_heat);
        Tneed_max = max(Tneed);
    };

    // without Tneed
    District(double energy, double Tref, int type, const char* type_heat, bool iren=false){  // one Building / overall District
        ren = iren;
        vector<double> vec_energy = {energy};
        vector<double> vec_Tref = {Tref};
        vector<int> vec_type = {type};
        vector<const char*> vec_type_heat = {type_heat};
        init_buildings(vec_energy, vec_Tref, vec_type, vec_type_heat);
    };

    // without Tneed
    District(double energy, double Tref, const char* typestring, const char* type_heat, bool iren=false){  // one Building / overall District
        ren = iren;
        vector<double> vec_energy = {energy};
        vector<double> vec_Tref = {Tref};
        vector<int> vec_type = { Building::get_type(typestring) };
        vector<const char*> vec_type_heat = {type_heat};
        init_buildings(vec_energy, vec_Tref, vec_type, vec_type_heat);
    };

    // with Tneed
    District(double energy, double Tref, double Tneed, int type, const char* type_heat, bool iren=false){  // one Building / overall District
        ren = iren;
        vector<double> vec_energy = {energy};
        vector<double> vec_Tref = {Tref};
        if (Tneed < 100)
            Tneed += 273.15;
        vector<double> vec_Tneed = {Tneed};
        vector<int> vec_type = {type};
        vector<const char*> vec_type_heat = {type_heat};
        init_buildings(vec_energy, vec_Tref, vec_Tneed, vec_type, vec_type_heat);
        Tneed_max = Tneed;
    };

    //with Tneed
    District(double energy, double Tref, double Tneed, const char* typestring, const char* type_heat, bool iren=false){  // one Building / overall District
        ren = iren;
        vector<double> vec_energy = {energy};
        vector<double> vec_Tref = {Tref};
        if (Tneed < 100)
            Tneed += 273.15;
        vector<double> vec_Tneed = {Tneed};
        vector<int> vec_type = { Building::get_type(typestring) };
        vector<const char*> vec_type_heat = {type_heat};
        init_buildings(vec_energy, vec_Tref, vec_Tneed, vec_type, vec_type_heat);
        Tneed_max = Tneed;
    };

    void init_buildings(vector<double> &energy, vector<double> &Tref, vector<int> &type, vector<const char*> &type_heat);
    void init_buildings(vector<double> &energy, vector<double> &Tref, vector<double> &Tneed, vector<int> &type, vector<const char*> &type_heat);

    void Display();

    /*
    vector<double> HeatinginTime(environment &Tempex);
    vector<double> HourlyHeatingPower(environment &Tempex,double Tref);
    vector<double> CoolinginTime(environment &Tempex);
    vector<double> HourlyCoolingPower(environment &Tempex,double Teref);
    vector<double> ECSinTime(vector <double> &Q);
     */

    vector<Building> B;
    vector<double> Qh,Qc,Qv,Qe,Qw,heq,Mw;   // data from buildings
    vector<double> Thref,Teref,Tdem;
    vector<int> type;
    vector<double> power;
    bool ren=false;     //TODO implement renovation status
    double Q_total, power_max, Tneed_max;

private:
    char sep;
    unsigned int nbTitleLines, nbTitleColumns;
    char* file_in;
    };

