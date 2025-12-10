//
// Created by lucile.schulthe on 11.10.2021.
//

#include <cmath>


#include "fluids/MWater.h"

#include "fluids/PengRobinson.h"
#include "math/maths.h"
#include"iofile/iofile.h"
#include"iofile/dynstring.h"
#include "state/TPhysicState.h"
#include "state/TSonicState.h"
#include "state/TThermoState.h"
#include "environment/environment.h"
//#include "consumer/Building.h"
//#include "consumer/District.h"
//#include "consumer/Norm_Building.h"
#include <iostream>

#include "BasicEnergySystem/exchanger/TExchanger.h"
#include "BasicEnergySystem/exchanger/TPlate.h"
#include "BasicEnergySystem/exchanger/TTubular.h"
#include "BasicEnergySystem/heatpump/THeatpump.h"
#include "iofile/iofilecsv.h"
#include "fluids/LeeKessler.h"

//#include "composant/PAC.h"


void displayVector(const std::vector<string> v){for (int i(0); i != v.size(); ++i)
        cout << "\n" << v[i];}

int main() {


    //-------------- Input file test-----------------//
    /*
    string line,data;
    ifstream infile;
    infile.open("Simu_para.txt");
    int Cols;
    infile >> data >> data>>Cols;

    cout <<data<<endl;

    vector <vector<string>> a;

    while (getline(infile,line)) {
        vector <string> z;
        for (int col = 0; col < Cols; ++col)
        {
            string n;
            infile >> n;
            cout <<n<<endl;
            z.push_back(n);
        }
        a.push_back(z);
    }

    cout<<endl<<a[0][0]<<endl; */
//-------------- Input file test-----------------//
/*
    const char* file_irrd = "dataBase/Posieux_irr_diff.txt";
    const char* file_irrw = "dataBase/rayonnement_v2.txt";
    const char* file_air = "dataBase/Posieux_kelvin.txt";
    const char* data = "mean";

    const char* file_out0 = "Results/moy/Env.txt";

    vector<vector <double>> ResultsEn;
    ResultsEn.clear();

    environment Tex((const char*) file_air,Tex.T,3);
    Tex.Moredata((const char*)file_irrd,Tex.Irr,data);
    //   Tex.Display();
     ResultsEn.push_back(Tex.T);
     ResultsEn.push_back(Tex.Irr);
    // Tmonth=Myear_month(Tex.T);
    //   Irrmonth=Myear_month(Tex.Irr);
     SetOutput2(ResultsEn,(char*)file_out0);
   cout<<Tex.T[0]<<endl;

*/


//-------------- Energy solaire -----------------//
/*
    const char* file_air = "dataBase/TempMoy.txt";
    const char* file_irr = "dataBase/Posieux_irr_glob.txt";
    const char* file_irrdiff = "dataBase/Posieux_irr_diff.txt";

    const char* data = "mean";
    const char* file_out = "Results/moy/Solar.txt";

    environment Tex((const char*) file_air,Tex.T,"temp");
    Tex.Moredata((const char*)file_irr,Tex.Irr,data);
    Tex.Moredata((const char*)file_irrdiff,Tex.Irr_diff,data);
    //   Tex.Display();

    solar thermo;
    vector<double>E(8760);
    E=thermo.thersolinTime(1,Tex,323.15);

    cout<<E[1]<<endl;
    SetOutput( E,(char*)file_out);


 */

    //-------------- Fluid test-----------------//

  FluidCore* fluid1;
    string fludiname ="R134a";
    //string fludiname ="Water";
    string method="CP";
    fluid_define(fluid1,fludiname);
    LeeKessler fluid2;
//    fluid_define(fluid3,fludiname);
    double P;
    double T=300;
//    double s=197096;
//    double v=0.0854;
    double u,h,x,v,s,ul,ug,hl,hg,sl,sg,vl,vg;
//    state.h=422865;
//    state.s=1747.89;
//    state.u=405201;
//    state.x=0.7;
//    state.P=fluid1->P_T_s(T);
//    state.P=P;
//    state.u=s;
//    state.v=v;
double Cv,Cp,mu,nu,Ro,L,Pr;
//    P=fluid1->P_T_s(T);
//    cout<<fluid1->P_T_s(T)<<endl;
//    cout<<fluid1->T_P_s(P)<<endl;

cout<<"fin init"<<endl<<endl;
    P=7.08e5;
//    v=fluid1->v_PT_l(P,T);



//    cout<<fluid1->P_T_s(300);
//    cout<<fluid1->T_P_s(8e5);

//    cout<<fluid2.Psat_T(300);
//    cout<<"T = "<<fluid1->T_P_s(P)<<endl;

//    fluid1->Prop_P_s(P,vl,vg,T,ul,ug,hl,hg,sl,sg);
//    fluid1->Prop_Tx_s(T,x,v,P,u,h,s);

//    fluid1->Prop_PT_l(P,T,v,u,h,s);
    // fluid1->CvCp_PT(P, T, Cv, Cp);

    // cout<<Cp<<endl;
    // P=1e5;
    //fluid1->CvCp_PT(P, T, Cv, Cp);
//    cout<<Cp<<endl;
//    fluid1->TransProp_PT_l(P, T,mu,nu,Ro,L,Pr);
//    cout<<mu<<' '<<nu<<' '<<Ro<<' '<<L<<' '<<Pr<<endl;


    v=fluid1->v_Tx_a(T,0);

    cout<<"v ="<<v<<endl;
cout<<"P ici "<<fluid1->P_vT(v,T)<<endl;

    cout<< " v est calcule "<<endl<<endl;

    fluid1->TransProp_PT_l(P, T,mu,nu,Ro,L,Pr);
    cout<<mu<<' '<<nu<<' '<<Ro<<' '<<L<<' '<<Pr<<endl;

    fluid1->TransProp_vT_l(v, 353.15,mu,nu,Ro,L,Pr);
    cout<<mu<<' '<<nu<<' '<<Ro<<' '<<L<<' '<<Pr<<endl;

    fluid1->CvCp_PT(P, 353.15,Cv,Cp);
    cout<<"test "<<Cp<<endl;

    cout<<P<<endl;

// ------------------Test central -----------------------//


/*
    Ccentral simu;


//    cout<<simu.connex[0][0]<<simu.connex[0][1]<<simu.connex[0][2]<<endl;
//    cout<<simu.connex[1][0]<<simu.connex[1][1]<<simu.connex[1][2]<<endl;

//    cout<<simu.countC<<endl;

//    cout<<simu.sst_connex[0][0]<<simu.sst_connex[0][1]<<simu.sst_connex[0][2]<<endl;


    simu.calcul();

//    cout<<"energy= "<<simu.network[0].energy<<endl;
//    cout<<"energy= "<<simu.network[1].energy<<endl;
//    cout<<"energy= "<<simu.network[2].energy<<endl;

//    simu.comp[0]->Display();
//    simu.comp[1]->Display();

    simu.Display();
*/
    //-------------- PAC test-----------------//
/*
    FluidCore* fluid;

    fluid=new PRfluid("R134a");

    double eta_com=0.9;
    double Thot=30+273.15;
    double Tcold=-20+273.15;
    double Tex=10+273.15;


    THeatpump PAC(fluid,eta_com, 0, 0);

    double COPreal =PAC.COP_real(Thot, Tcold,Tex);
    cout <<"COP real : "<<COPreal << endl;
    cout <<"COP real : "<<PAC.COP_real(Thot,Tcold) << endl;

    //   cout <<"COP max : "<<PAC.DeltaS_Cp(Thot,Tcold) << endl;

*/
//--------------Turbine-----------------//






    //-------------- Network test-----------------//

    //---------Network init
/*
    Csimunetwork simu;
    cout<<simu.network.pipe[0].Diam<<endl;

    cout<<simu.network.e[0].energy<<endl;

    cout<<simu.network.branch[0].e[0]->energy<<endl;

*/
// ---- Heat exchanger ------------//

/*
    const char* file_inFLUID = "fluids/LKDataFluidsNEW5.txt";
    //DECLARATION DES FLUIDES
    MWater fluidMain("H2O");
    MWater fluidPlateHot("H2O");

    //DECLARATION DE L ECHANGEUR
    TPlate plateOne("P1",TExchanger::COUNTER, TExchanger::DEFAULT);
    plateOne.init();

    //On initialise les fluides chauds et froids (récuperer les données des états)
    plateOne.Hot.init(fluidPlateHot);
    plateOne.Cold.init(fluidMain);

    plateOne.performance();

    cout<<"---------------------Plates Heat Exchanger n°1-------------------------------"<<endl;
    plateOne.display();


    /////Evaporator
    MWater fluidEvapHot("H2O");
    TTubular evap("Cond",TExchanger::COUNTER, TExchanger::EVAP, TTubular::HOT);

    evap.Cold.pStateIn.SetT(plateOne.Cold.pStateOut.T+273.15);
    evap.Cold.pStateIn.SetP(1E5);
    evap.Cold.Massflow=plateOne.Cold.Massflow;

    evap.Hot.pStateIn.SetT(423.15);
    evap.Hot.pStateOut.SetT(323.15);
    evap.Hot.pStateIn.SetP(1E5);



    double Ta;
    Ta=((evap.Cold.pStateIn.GetT()+evap.Cold.pStateOut.GetT())/2+(evap.Hot.pStateIn.GetT()+evap.Hot.pStateOut.GetT())/2)/2;
    evap.Cold.pStateWall.SetT(Ta);
    evap.Hot.pStateWall.SetT(Ta);

    evap.Hot.init(fluidEvapHot);
    evap.Cold.init(fluidMain);

    evap.performance();


    evap.display();

    evap.Q=1;

    evap.performance();
    cout<<evap.Q<<endl;
    evap.display();

*/



    //-------------- END-----------------//
    return 0;
}
