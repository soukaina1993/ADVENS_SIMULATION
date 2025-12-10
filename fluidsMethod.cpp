//
// Created by lucile.schulthe on 20.06.2022.
//


#include <cmath>
#include <string>


#include "math/maths.h"
#include "iofile/iofile.h"
#include"iofile/dynstring.h"
#include "fluids/FluidCore.h"
#include "state/TState.h"


#include <iostream>
void GetState(FluidCore *&fluid);

int main(int argc,char *argv[]){
    if (argc<4) {
        cout << "error not valide parameter. 0: Method, 1: Fluid name, 2: Required  parameter, 3: Input name variable1, 4: Input value variable1, 5: Input name variable2, 6: Input value variable2 "<<endl;
    }else {
        string Method = argv[1];
        string fluidname = argv[2];
        string varout = argv[3];

        cout<<"Methode demande :"<<Method<<endl;
        FluidCore *fluid;
        fluid_define(fluid, fluidname, Method);

        if (varout=="Pc")
            cout<<fluid->getPc()<<endl;
        else if(varout=="Tc")
            cout<<fluid->getTc()<<endl;
        else if(varout=="curve"){
            vector<double> result;


            string varout2 = argv[4];
            string par1 = argv[5];
            double arg1 = stod(argv[6]);


            double arg2 = stod(argv[8]);
            double arg3 = stod(argv[9]);

            string par2 = argv[7];
            for(double i=arg2; i<arg3; i+=(arg3-arg2)/20){
                TState state;
                state.assign(par1, arg1);
                state.assign(par2, i);
                state.init_all(fluid);
                state.Display(varout2);
            }

        }

        else{
            if (argc == 8) {
                TState state;
                string par1 = argv[4];
                string par2 = argv[6];
                double arg1 = stod(argv[5]);
                double arg2 = stod(argv[7]);
                state.assign(par1, arg1);
                state.assign(par2, arg2);
                state.init_all(fluid);
                if(varout=="Cp"){
                    double Cp,Cv;
                    fluid->CvCp_PT(state.P,state.T,Cv,Cp);
                    cout<<Cp<<endl;
                }
                else if(varout=="TransProp"){
                    double mu,nu,Ro,L,Pr;
                    fluid->TransProp_vT_a(state.v,state.T,mu,nu,Ro,L,Pr);
                    cout<<mu<<" "<<nu<<" "<<Ro<<" "<<L<<" "<<Pr<<endl;
                }
                else if(varout=="Sonic"){
                    cout<< fluid->sv_vT(state.v,state.T)<<endl;
                }
                else
                    state.Display(varout);

            }  else if (argc ==4){
                vector <TState> state;
                GetState(fluid); //From File
            }else {
                cout << "error not valide parameter.0: Method, 1: Fluid name, 2: Required  parameter, 3: Input name variable1, 4: Input value variable1, 5: Input name variable2, 6: Input value variable2 "<<endl;
            }
        }

    }
    return 0;
};


void GetState(FluidCore *&fluid) {
    vector<vector<string>> table;
    vector<vector<double>> Results;

    table=read_record("data_temp");
    TState state_temp;
    for(int i=1;i<table.size(); i++) {

        state_temp.T = stod(table[i][0]);
        state_temp.P = stod(table[i][1]);
        state_temp.v = stod(table[i][2]);
        state_temp.h = stod(table[i][3]);
        state_temp.s = stod(table[i][4]);
        state_temp.x = stod(table[i][5]);

        state_temp.init_all(fluid);

        Results[i][0]=state_temp.T;
        Results[i][1]=state_temp.P;
        Results[i][2]=state_temp.v;
        Results[i][3]=state_temp.h;
        Results[i][4]=state_temp.s;
        Results[i][5]=state_temp.x;

    }
    SetOutput(Results, (char *) "results.txt");
}