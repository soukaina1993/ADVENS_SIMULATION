//
// Created by lucile.schulthe on 14.02.2022.
//

#include "network.h"



Tnetwork::Tnetwork(double &TDHN, double &dT, double &ResRatio,int &pressure,int &alt_plant,int &alt_max,int &alt_min, char &fluidin1):
    dT(dT), ResRatio(ResRatio),alt_plant(alt_plant),alt_max(alt_max),alt_min(alt_min) {
    flow->pStateIn.T=TDHN;
    flow->StateIn->P=pressure;
    flow->pStateOut.P= flow->StateIn->P-( alt_max-alt_plant)*9.81;
    char* fludiname= (char*)fluidin1;
    fluid_define(fluid,fludiname);
    double cv;
    flow->StateCalculation(flow->StateIn,fluid);
    flow->PhysicCalculation(flow->pStateIn);
};
Tnetwork::Tnetwork(double &TDHN, double &dT, double &ResRatio,int &pressure,int &alt_plant,int &alt_max,int &alt_min, string &fluidin1):
        dT(dT), ResRatio(ResRatio),alt_plant(alt_plant),alt_max(alt_max),alt_min(alt_min) {
    flow->pStateIn.T=TDHN;
    flow->StateIn->P=pressure;
    flow->pStateOut.P=flow->StateIn->P-( alt_max-alt_plant)*9.81;
    flow->pStateOut.P=flow->pStateIn.P-dP;
    int n = fluidin1.length();
    char char_array[n + 1];
    strcpy(char_array, fluidin1.c_str());
    char* fludiname= (char*)char_array;
    fluid_define(fluid,fludiname);
    double cv;
    flow->StateCalculation(flow->StateIn,fluid);
    flow->PhysicCalculation(flow->pStateIn);
};

void Tnetwork::init_fluid(char &fluidin1){
    char* fludiname= (char*)fluidin1;
    fluid_define(fluid,fludiname);
    flow->init(fluid);
};

void Tnetwork::init_fluid(string &fluidin1){
    int n = fluidin1.length();
    char char_array[n + 1];
    strcpy(char_array, fluidin1.c_str());
    char* fludiname= (char*)char_array;
    fluid_define(fluid,fludiname);
    flow->init(fluid);
};

void Tnetwork::branch_init(char* &file_in){
    vector<vector<string>> table;
    table=read_record(file_in);

    int nb_br=table[0].size()-3; // nombre de branche
    int nb_lin=table.size(); // nombre de ligne
    int nb_couche=(nb_lin-3)/4; // nombre de sous-couches d'electron

    for(int i=0; i<nb_br; i++){
        TBranch temp_branch;
        for(int j=0; j<nb_couche;j++) {

            if (table[3+j][i + 2].size() > 0 && table[3+j][i + 2]!="0") {
                sst* e_temp;
                e_temp=&e[stoi(table[3 +j][i + 2]) - 1];
                temp_branch.e.push_back(e_temp);
            }
            if (table[3+j+nb_couche][i + 2].size() > 0 && table[3+j][i + 2]!="0") {
                TSegment seg_temp;
                seg_temp.iPipe = seg_temp.oPipe = &pipe[stoi(table[3+j+nb_couche][i + 2]) - 1];
                temp_branch.seg.push_back(seg_temp);
            }
        }



        branch.push_back(temp_branch);
    }

}





void   Tnetwork::init (environment & ext){
    N=(branch.size()-1)/2;
    B=branch.size();
    Z=e.size();
    A=N+Z;

    for(int i=0; i<Z;i++)
        e[i].calcul_power(ext);


    for(int i=0; i<N;i++){
        if (branch[i].N_num!=0)
        branch[i].K=branch[i].e.size()+branch[0].K-1;

    }


}


