//
// Created by cornelia.blanke on 29.08.2023.
//

/*
 * Compile in MATLAB with:
 * mex ADVENS_Props.cpp ../fluids/*.cpp ../math/*.cpp ../iofile/*.cpp ../state/*.cpp ../Flow/*.cpp
 * mex ADVENS_Props.cpp ../fluids/*.cpp ../math/*.cpp ../iofile/*.cpp ../state/*.cpp ../Flow/*.cpp -llibCoolProp -L"C:\Users\cornelia.blanke\Documents\adv\mex"
 *
 * Usage: ADVENS_Props('Ro','P',101325,'T',80+273.15,'H2O','CP')
 */

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <iostream>
#include <fstream>
#include <vector>

#include "../fluids/PengRobinson.h"
#include "../fluids/ExtCoolProp.h"
#include "../Flow/Flow_physic.h"


using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;
    shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    ostringstream stream;
    const vector<string> constants{"Pc", "Tc"};
    const vector<string> props{"v", "P", "T", "h", "s", "x", "u"};
    const vector<string> props_out{"ViscDyn", "ViscCin", "Ro", "ThermCond", "Prandtl", "Cp"};

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {

        // Validate and get arguments
        string out;
        double arg1, arg2;
        string prop1, prop2, fluidname;
        string fluidlib;
        vector<double> fluidvar;
        fluidvar.clear();
        vector<double> results(props.size() + props_out.size());

        getArguments(outputs, inputs, out, prop1, arg1, prop2, arg2, fluidname, fluidlib, fluidvar);

        // Define fluid and flow
        FluidCore *fluid;

        if (fluidlib == "CP")                       // use CoolProp
            fluid = new ExtCoolProp(fluidname.c_str());
        else if (fluidlib == "MW")                  // use MWater
            fluid = new MWater;
        else if (fluidlib == "PR")                // use PengRobinson
            fluid = new PRfluid((char *) fluidname.c_str(), fluidvar[0], fluidvar[1], fluidvar[2], fluidvar[3],
                                fluidvar[4], fluidvar[5], fluidvar[6], fluidvar[7], fluidvar[8],
                                fluidvar[9], fluidvar[10], fluidvar[11]);
        else                                      // use LeeKessler
            fluid = new LeeKessler((char *) fluidname.c_str(), fluidvar[0], fluidvar[1], fluidvar[2], fluidvar[3],
                                   fluidvar[4], fluidvar[5], fluidvar[6], fluidvar[7], fluidvar[8],
                                   fluidvar[9], fluidvar[10], fluidvar[11], fluidvar[12], fluidvar[13],
                                   fluidvar[14]);


        if (out == "Pc")
            outputs[0] = factory.createScalar(fluid->getPc());
        else if (out == "Tc")
            outputs[0] = factory.createScalar(fluid->getTc());
        else {
            double Psat;
            if ((prop1 == "P" && prop2 == "T" && (abs(((Psat=fluid->P_T_s(arg2)) - arg1)/arg1) < 1e-3)) ||
                    (prop1 == "T" && prop2 == "P" && (abs(((Psat=fluid->P_T_s(arg1)) - arg2)/arg2) < 1e-3))) {
                displayWarningOnMATLAB("P is close to saturation pressure " + to_string(Psat) +
                                       "; remember that (P,T) cannot be used for biphasic state");
            }

            Flow_physic *flow = new Flow_physic(fluid);
//        getProperties(flow, results);
//        displayFluidProperties(results, "old values");

            // Implement function
            setValue(flow->pStateIn, prop1, arg1);
            setValue(flow->pStateIn, prop2, arg2);

            if (!(prop1 == "P" || prop1 == "T" || prop2 == "P" || prop2 == "T"))
                displayErrorOnMATLAB("Always enter P or T");

            flow->StateCalculation(flow->StateIn);
            if (flow->pStateIn->GetV() == 0 || flow->pStateIn->GetT() == 0)
                displayErrorOnMATLAB("State calculation failed");
            else {
                flow->PhysicCalculation_vT(flow->pStateIn);
                double Cv;
                flow->fluid->CvCp_vT(flow->vIn(), flow->TIn(), Cv, flow->Cp);
            }

            getProperties(flow, results);
            //displayFluidProperties(results, "new values");

            // Assign outputs
            if (out == "all")
                // returns an array
                outputs[0] = factory.createArray<double>({1, results.size()}, results.data(),
                                                         results.data() + results.size());
            else {
                // returns a scalar
                int idx = 0;
                while (out != props[idx] && idx < props.size())
                    idx++;
                if (idx == props.size()) {      // not yet found
                    while (out != props_out[idx - props.size()] && idx < props.size() + props_out.size())
                        idx++;
                }
                outputs[0] = factory.createScalar(results[idx]);
            }

            // Destructors
            delete flow;
        }
        delete fluid;
    }

    void getArguments(ArgumentList outputs, ArgumentList inputs,
                      string &out, string &prop1, double &arg1, string &prop2, double &arg2,
                      string &fluidname, string &fluidlib, vector<double> &fluidvar) {

        // check number of inputs
        if (inputs.size() == 3) {
            // check each input
            if (!matlab_to_string(inputs[0], out))
                displayErrorOnMATLAB("First argument not a char array or string");
            if (find(constants.begin(), constants.end(), out) == constants.end())
                displayErrorOnMATLAB("First argument is not a valid property");

            if (!matlab_to_string(inputs[1], fluidname))
                displayErrorOnMATLAB("Second argument not a char array or string");

            if (!matlab_to_string(inputs[2], fluidlib))
                displayErrorOnMATLAB("Third argument not a char array or string");
            else if (fluidlib != "CP" && fluidlib != "MW" && fluidlib != "PR" && fluidlib != "LK")
                displayErrorOnMATLAB("Third argument is not a valid fluid library");
            else if (fluidlib == "CP" && !is_in_CoolProp(fluidname))
                displayErrorOnMATLAB("Fluid not found in CoolProp database");
            else if (fluidlib == "MW" && fluidname != "H2O" && fluidname != "Water")
                displayErrorOnMATLAB("MWater is only for H2O");
            else if (fluidlib == "PR" || fluidlib == "LK") {
                if (fluidname == "H2O" || fluidname == "Water")
                    displayErrorOnMATLAB("Do not use PengRobinson or LeeKessler model for H2O");
                else {
                    fluidvar = get_DataFluids(fluidname);
                    if (fluidvar.empty())
                        displayErrorOnMATLAB("Fluid not found in DataFluids.txt");
                }
            }
        }

        else if (inputs.size() == 7) {
            // check each input
            if (!matlab_to_string(inputs[0], out))
                displayErrorOnMATLAB("First argument not a char array or string");
            if (find(props.begin(), props.end(), out) == props.end() &&
                find(props_out.begin(), props_out.end(), out) == props_out.end() && (out != "all"))
                displayErrorOnMATLAB("First argument is not a valid property");

            if (!matlab_to_string(inputs[1], prop1))
                displayErrorOnMATLAB("Second argument not a char array or string");
            if (find(props.begin(), props.end(), prop1) == props.end())     // prop1 not found in props
                displayErrorOnMATLAB("Second argument is not a valid property");

            if (!matlab_to_double(inputs[2], arg1))
                displayErrorOnMATLAB("Third argument not a double");

            if (!matlab_to_string(inputs[3], prop2))
                displayErrorOnMATLAB("Forth argument not a char array or string");
            if (find(props.begin(), props.end(), prop2) == props.end())     // prop2 not found in props
                displayErrorOnMATLAB("Forth argument is not a valid property");

            if (!matlab_to_double(inputs[4], arg2))
                displayErrorOnMATLAB("Fifth argument not a double");

            if (!matlab_to_string(inputs[5], fluidname))
                displayErrorOnMATLAB("Sixth argument not a char array or string");

            if (!matlab_to_string(inputs[6], fluidlib))
                displayErrorOnMATLAB("Seventh argument not a char array or string");
            if (fluidlib != "CP" && fluidlib != "MW" && fluidlib != "PR" && fluidlib != "LK")
                displayErrorOnMATLAB("Seventh argument is not a valid fluid library");
            else if (fluidlib == "CP" && !is_in_CoolProp(fluidname))
                displayErrorOnMATLAB("Fluid not found in CoolProp database");
            else if (fluidlib == "MW" && fluidname != "H2O" && fluidname != "Water")
                displayErrorOnMATLAB("MWater is only for H2O");
            else if (fluidlib == "PR" || fluidlib == "LK") {
                if (fluidname == "H2O" || fluidname == "Water")
                    displayErrorOnMATLAB("Do not use PengRobinson or LeeKessler model for H2O !");
                else {
                    fluidvar = get_DataFluids(fluidname);
                    if (fluidvar.empty())
                        displayErrorOnMATLAB("Fluid not found in DataFluids.txt");
                }
            }
        } else
            displayErrorOnMATLAB("Wrong number of inputs");

        // check size of output
        if (outputs.size() > 1)
            displayErrorOnMATLAB("Only one output is returned");

    }

    bool matlab_to_string(matlab::data::Array in, string& out){
        if (in.getType() == ArrayType::CHAR) {
            const CharArray arr = in;
            out = arr.toAscii();
            return true;
        } else if (in.getType() == ArrayType::MATLAB_STRING) {
            string tmp = in[0];
            out = tmp;
            return true;
        }
        else
            return false;
    }

    bool matlab_to_double(matlab::data::Array in, double& out) {
        if (in.getType() != ArrayType::DOUBLE)
            return false;
        else {
            out = in[0];
            return true;
        }
    }

    bool is_in_CoolProp(string &fluidname) {
        if (fluidname == "H2O")     // rename fluid
            fluidname = "Water";
        stringstream allfluid(get_global_param_string("FluidsList"));
        string fluid;
        while (getline(allfluid, fluid, ',')) {
            if (fluid == fluidname)
                return true;
        }
        return false;
    }

    vector<double> get_DataFluids(const string& fluidname){
        vector<vector<string>> fluidtable = read_csv((char *) "DataFluids.txt", '\t');
        vector<double> fluidvar;
        for (auto fluidrow: fluidtable) {
            if (fluidname == fluidrow[0]) {
                for (int j = 1; j < fluidrow.size(); j++)
                    fluidvar.push_back(stod(fluidrow[j]));
                break;
            }
        }
        return fluidvar;
    }

    void setValue(TPhysicState* pState, const string& prop, const double& arg){
        if (prop == "v")
            pState->SetV(arg);
        else if (prop == "P")
            pState->SetP(arg);
        else if (prop == "T")
            pState->SetT(arg);
        else if (prop == "h")
            pState->SetH(arg);
        else if (prop == "s")
            pState->SetS(arg);
        else if (prop == "x")
            pState->SetX(arg);
//        else if (prop == "u")
//            pState->u = arg;
    }

    void getProperties(Flow_physic* flow, vector<double>& results)
    {
        results[0] = flow->pStateIn->v;
        results[1] = flow->pStateIn->P;
        results[2] = flow->pStateIn->T;
        results[3] = flow->pStateIn->h;
        results[4] = flow->pStateIn->s;
        results[5] = flow->pStateIn->x;
        results[6] = flow->pStateIn->u;
        results[7] = flow->pStateIn->ViscDyn;
        results[8] = flow->pStateIn->ViscCin;
        results[9] = flow->pStateIn->Ro;
        results[10] = flow->pStateIn->ThermCond;
        results[11] = flow->pStateIn->Prandtl;
        results[12] = flow->Cp;
    }

    void displayFluidProperties(const vector<double>& results, string title=""){
        stream << title << '\n';
        for(int i=0; i<props.size(); i++)
            stream << props[i] << ":\t" << results[i] << '\n';
        for(int i=0; i<props_out.size(); i++)
            stream << props_out[i] << ":\t" << results[props.size() + i] << '\n';
        displayOnMATLAB(stream);
    }


    /*** helper functions --> could be summarized in a #include-file if needed in several mex-files ***/
    void displayOnMATLAB(ostringstream& ostream) {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
                         vector<matlab::data::Array>({ factory.createScalar(ostream.str()) }));
        // Clear stream buffer
        ostream.str("");
    }

    void displayOnMATLAB(string str) {
        // Pass string to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
                         vector<matlab::data::Array>({ factory.createScalar(str + "\n") }));
    }

    void displayOnMATLAB(double d) {
        displayOnMATLAB(to_string(d));
    }

    void displayErrorOnMATLAB(string str){
        matlabPtr->feval(u"error", 0,
                         vector<matlab::data::Array>({ factory.createScalar(str + "\n") }));
    }

    void displayWarningOnMATLAB(string str){
        matlabPtr->feval(u"warning", 0,
                         vector<matlab::data::Array>({ factory.createScalar(str + "\n") }));
    }

    vector<vector<string>> read_csv(const char * inputfile, char sep) {
        vector<vector<string>> content;
        vector<string> row;
        string line, word;

        fstream file(inputfile, ios::in);
        if (file.is_open()) {
            while (getline(file, line)) {
                row.clear();
                stringstream str(line);
                while (getline(str, word, sep))
                    row.push_back(word);
                content.push_back(row);
            }
        } else
            displayErrorOnMATLAB("Could not open the file " + string(inputfile));
        file.close();
        return content;
    }
};
