//
// Created by cornelia.blanke on 04.03.2024.
//

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "./state/TState.h"
#include "./fluids/PengRobinson.h"
#include "./fluids/MWater.h"
#include "./fluids/LeeKessler.h"
#include "./fluids/ExtCoolProp.h"

using namespace std;

// declarations, implementation below
void setValue(TState& pState, const string& prop, const double& arg);
string read_datafluids(const char * inputfile, char sep, int titlelines);
vector<string> split(const string &s, char sep);

/***
 * call:
 * ADVENS_FLUID.exe FluidsList method --> get all available fluids
 *
 * ADVENS_FLUID.exe inprop1 invalue1 inprop2 invalue2 fluid method --> get properties
 * (invalue1 and invalue2 may be comma-separated strings)
 * ***/

int main(int argc, char *argv[]){
    const char* base_path = getenv("MY_BASE_PATH");
    string totalpath;
    if (base_path != nullptr)
        totalpath = string(base_path) + "/dataBase/DataFluids.txt";
    else
        totalpath = "dataBase/DataFluids.txt";
    const char* file_in = totalpath.c_str();


    // get all fluids
    if (argc==3 && strcmp(argv[1], "FluidsList")==0){
        string fluidlist;
        if (strcmp(argv[2], "CoolProp")==0)
            fluidlist = get_global_param_string("FluidsList");
        else if ((strcmp(argv[2], "PRFluid")==0) || (strcmp(argv[2], "LeeKessler")==0)){
            // read DataFluids.txt
            fluidlist = read_datafluids(file_in, '\t', 2);
        }
        cout << fluidlist << endl;
    }

    // compute properties
    else if (argc==7){
        FluidCore* fluid = nullptr;

        // create fluid
        if (strcmp(argv[6], "PRFluid")==0)
            fluid = new PRfluid(argv[5], file_in);
        else if (strcmp(argv[6], "MWater")==0)
            fluid = new MWater(argv[5], file_in);
        else if (strcmp(argv[6], "LeeKessler")==0)
            fluid = new LeeKessler(argv[5], file_in);
        else if (strcmp(argv[6], "CoolProp")==0)
            fluid = new ExtCoolProp(argv[5]);
        else
            cerr << "Error: Invalid method" << endl;


        vector<string> arg2 = split(argv[2], ',');
        vector<string> arg4 = split(argv[4], ',');

        for (const auto & i : arg2) {
            for (const auto &j: arg4) {
                TState pState;
                // set input
                setValue(pState, argv[1], stod(i));
                setValue(pState, argv[3], stod(j));

                // compute all properties
                pState.init_all_x(fluid);

                //output
                cout << pState.P << " ";       //Pa
                cout << pState.T << " ";       //K
                cout << pState.v << " ";       //m^3/kg
                cout << pState.h << " ";       //J/kg
                cout << pState.s << " ";       //J/kg.K
                cout << pState.u << " ";       //J/kg
                cout << pState.x << " ";       //-
                cout << fluid->Pc << " ";       //Pa
                cout << fluid->Tc << " ";       //K
                cout << fluid->vc << " ";       //m^3/kg
                cout << fluid->hc << " ";       //J/kg
                cout << fluid->sc << " ";       //J/kg.K
                cout << fluid->uc << "\n";       //J/kg
            }
        }

        delete fluid;
}
    else
        cerr << "Error: Wrong input" << endl;

    return 0;
}

/*** ------------------------------------------------------------------------------ ***/


void setValue(TState& pState, const string& prop, const double& arg){
    if (prop == "v")
        pState.SetV(arg);
    else if (prop == "P")
        pState.SetP(arg);
    else if (prop == "T")
        pState.SetT(arg);
    else if (prop == "h")
        pState.SetH(arg);
    else if (prop == "s")
        pState.SetS(arg);
    else if (prop == "x")
        pState.SetX(arg);
//        else if (prop == "u")
//            pState->u = arg;
}

string read_datafluids(const char * inputfile, char sep, int titlelines) {
    string content;
    fstream file(inputfile, ios::in);
    if (file.is_open()) {
        string line, word;
        int count=0;
        while (getline(file, line)) {
            if (count == titlelines) {
                stringstream str(line);
                getline(str, word, sep);
                content = word;
            }
            else if (count > titlelines) {
                stringstream str(line);
                getline(str, word, sep);
                content += "," + word;
            }
            count++;
        }
        file.close();
    }
    return content;
}

vector<string> split(const string &s, char sep) {
    vector<string> result;
    stringstream ss (s);
    string item;
    while (getline (ss, item, sep))
        result.push_back(item);
    return result;
}