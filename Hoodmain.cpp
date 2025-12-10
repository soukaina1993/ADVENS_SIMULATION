#include <cmath>
#include <string>
#include "fluids/PengRobinson.h"
#include "fluids/PRR.h"
#include "fluids/FluidCore.h"
#include "math/maths.h"
#include"iofile/iofile.h"
#include"iofile/dynstring.h"
#include "state/TPhysicState.h"
#include "state/TSonicState.h"
#include "state/TThermoState.h"
#include "fluids/MWater.h"
#include "environment/environment.h"
#include "consumer/Building.h"
#include "consumer/District.h"
#include "consumer/Norm_Building.h"
#include "consumer/solar.h"
#include <iostream>
#include "iofile/iofilecsv.h"


void input_analyse(vector<vector<string>> table,int counttype[5]);

int main(){

    const char* file_air = "dataBase/Posieux_kelvin.txt";

    const char* file_out0 = "Results/moy/Env.txt";
    const char* file_out1 = "Results/moy/PRuiss.txt";
    const char* file_out2 = "Results/moy/PEauxUsees.txt";
    const char* file_out3 = "Results/moy/PhQuartier.txt";
    const char* file_out4 = "Results/moy/PwQuartier.txt";
    const char* file_out5 = "Results/moy/PcQuartier.txt";
    const char* file_out6 = "Results/moy/Ev.txt";
    const char* file_out7 = "Results/moy/COPWater.txt";
    const char* file_out8 = "Results/moy/TWater.txt";
    const char* file_irr = "dataBase/Posieux_irr_glob.txt";
    const char* file_irrdiff = "dataBase/Posieux_irr_diff.txt";
    const char* file_rui = "InputData/Ruisselement_Q.txt";
    const char* file_usee = "InputData/EauxUsees_Quartier1Copie.txt";
    const char*file_hood1="InputData/Fr_100.txt";
    const char*file_hood2="InputData/Fr_200.txt";
    const char*file_hood3="InputData/Fr_300.txt";
    const char*file_hood4="InputData/Fr_400.txt";
    const char*file_hood5="InputData/Fr_BF_PAC.txt";
    const char*file_hood6="InputData/Fr_BF_PAC.txt"; // Ideal meme chose que le PAC

    const char* data = "mean";


//    const char* Q = "Perolles";
//

    vector<vector <double>> ResultsQR,ResultsEn,ResultsPh,ResultsPc,ResultsQT,  ResultsPw,mr, mt,Twater,COPWater;
    ResultsEn.clear();
    ResultsPh.clear();
    ResultsQT.clear();
    ResultsQR.clear();
    ResultsPc.clear();
    ResultsPw.clear();

    vector<double> P,W,Waterday,Tmonth,QR(8760),QT(24),Irrmonth;
    vector<double> EauxUsees,Ruiss,mr_s,mt_s,Twater_s,COPWater_s;
double Cp,Cv;

    int Time;
    int nb_Q=6;

//    double P_Qrmoy= 0.018,B_Qrmoy=0.011,E_Qrmoy=0.007;
//    double P_Qtmoy= 0.012,B_Qtmoy=0.006,E_Qtmoy=0.004;

    double Qrmoy= 0.018;
    double Qtmoy= 0.012;

    cout<<"Debut"<<endl;

//------------- SET ENVIRONNEMENT -----------------------///

    environment Tex((const char*) file_air,Tex.T,data);
    Tex.Moredata((const char*)file_irr,Tex.Irr,data);
    Tex.Moredata((const char*)file_irrdiff,Tex.Irr_diff,data);
    //   Tex.Display();
    ResultsEn.push_back(Tex.T);
    //   ResultsEn.push_back(Tex.Irr);
    // Tmonth=Myear_month(Tex.T);
    //   Irrmonth=Myear_month(Tex.Irr);
     SetOutput(ResultsEn,(char*)file_out0);
    cout<<"Env"<<endl;

//------------- SET WATER -----------------------///
    MWater water;
    double Cp0;
    double T0=273.15;
    water.CvCp_vT_l(273.15,Cp0,Cv);
    double T_EU=273.15+11, Cp11;
    water.CvCp_vT_l(T_EU,Cp11,Cv);
    // double Rho_EU=water.Rho_T_l(T_EU);
    double Rho_0=water.Rho_T_l(T0);
    cout<<"Water"<<endl;
//------------- CALCUL DE L ENERGIE RUISSELEMENT -----------------------///
/*
    ResultsQR.push_back(Tmonth);
    Getdata((const char *) file_rui, Waterday, Time, 1);
    ResultsQR.push_back(Waterday);

    for(int Q_index(0);Q_index<3;Q_index++) {
        Getdata((const char *) file_rui, Ruiss, Time, Q_index+2);
        for (int i(0); i < Tmonth.size(); i++) {
            if (Tmonth[i] > 273.15) {
                water.CvCp_vT_l(Tmonth[i], Cp[i], Cv[i]);
                QR[i] = ((Cp[i]+Cp0)/2 * (Tmonth[i] - 273.15)*Rho_0/3.6e6) * Ruiss[i];
            }
        }
        mr_s=Mmonth_year(Ruiss);
        mr.push_back(mr_s);
        ResultsQR.push_back(QR);
    }

    SetOutput(ResultsQR,(char*)file_out1);
*/
    for(int Q_index(0);Q_index<nb_Q;Q_index++) {
        Getdata((const char *) file_rui, Ruiss, Time, Q_index+1);
        for (int i(0); i < Tex.T.size(); i++) {
            if (Tex.T[i] > 273.15) {
                water.CvCp_vT_l(Tex.T[i], Cp, Cv);
                QR[i] = ((Cp+Cp0)/2 * (Tex.T[i] - 273.15)*Rho_0/3.6e6) * Ruiss[i];
            }
        }
        mr.push_back(Ruiss);
        ResultsQR.push_back(QR);
    }

    SetOutput2(ResultsQR,(char*)file_out1);
    cout<<"RUISSELEMENT"<<endl;


    //------------- CALCUL DE L ENERGIE EAUX USEES -----------------------///


    double a=(Cp11+Cp0)/2 *(T_EU - 273.15)*Rho_0/3.6e6;

    for(int Q_index(0);Q_index<nb_Q;Q_index++) {
        Getdata((const char *) file_usee, EauxUsees, Time, Q_index+1);
        for (int i(0); i < EauxUsees.size(); i++) {
            QT[i] =  a* EauxUsees[i];
        }
        mt_s=Mhour_year(EauxUsees);
        mt.push_back(mt_s);
        ResultsQT.push_back(QT);
    }
    SetOutput2(ResultsQT,(char*)file_out2);


    cout<<"Eaux usees"<<endl;

    //---------Temperature EAU---------//
    for(int j(0); j<nb_Q;j++) {
        Twater_s.clear();
        for (int i(0); i < 8760; i++) {
            if (Tex.T[i] > 273.15)
                Twater_s.push_back((mt[j][i] * 284.15 + mr[j][i] * Tex.T[i]) / (mr[j][i] + mt[j][i]));
            else
                Twater_s.push_back((mt[j][i] * 284.15 + mr[j][i] * 273.15) / (mr[j][i] + mt[j][i]));
        }
        Twater.push_back(Twater_s);
    }
    SetOutput2(Twater,(char*)file_out8);
    cout<<"Temperature eau"<<endl;
//------------- DECLARATION QUARTIER -----------------------///

    norm_ecs normEcs;
    norm_energy E;
    bool ren=0;
    District Q1(ren,(char*) file_hood1);
    District Q2(ren,(char*) file_hood2);
    District Q3(ren,(char*) file_hood3);
    District Q4(ren,(char*) file_hood4);
    District Q5(ren,(char*) file_hood5);
    District Q6(ren,(char*) file_hood6);


//    Q1.Display();
//    Q2.Display();
//    Q3.Display();
    cout<<"Declaration quartier"<<endl;
//------------- CALCUL DE BESOIN CHAUFFAGE -----------------------///

    P=Q1.HeatinginTime(Tex);
    ResultsPh.push_back(P);
    P=Q2.HeatinginTime(Tex);
    ResultsPh.push_back(P);
    P=Q3.HeatinginTime(Tex);
    ResultsPh.push_back(P);
    P=Q4.HeatinginTime(Tex);
    ResultsPh.push_back(P);
    P=Q5.HeatinginTime(Tex);
    ResultsPh.push_back(P);
    P=Q6.HeatinginTime(Tex);
    ResultsPh.push_back(P);
    SetOutput2(ResultsPh,(char*)file_out3);

    cout<<"heating"<<endl;


    //------------- CALCUL DE BESOIN EAU CHAUDE -----------------------///

    norm_ecs ECS;
    W=Q1.ECSinTime(ECS,Q1.Qw);
    ResultsPw.push_back(W);
    W=Q2.ECSinTime(ECS,Q2.Qw);
    ResultsPw.push_back(W);
    W=Q3.ECSinTime(ECS,Q3.Qw);
    ResultsPw.push_back(W);
    W=Q4.ECSinTime(ECS,Q4.Qw);
    ResultsPw.push_back(W);
    W=Q5.ECSinTime(ECS,Q5.Qw);
    ResultsPw.push_back(W);
    W=Q6.ECSinTime(ECS,Q6.Qw);
    ResultsPw.push_back(W);
    SetOutput2(ResultsPw,(char*)file_out4);

    cout<<"heating water"<<endl;
    //------------- CALCUL DE BESOIN FROID -----------------------///

    P=Q1.CoolinginTime(Tex);
    ResultsPc.push_back(P);
    P=Q2.CoolinginTime(Tex);
    ResultsPc.push_back(P);
    P=Q3.CoolinginTime(Tex);
    ResultsPc.push_back(P);
    P=Q4.CoolinginTime(Tex);
    ResultsPc.push_back(P);
    P=Q5.CoolinginTime(Tex);
    ResultsPc.push_back(P);
    P=Q6.CoolinginTime(Tex);
    ResultsPc.push_back(P);
    SetOutput2(ResultsPc,(char*)file_out5);

    Q1.Display();
    Q2.Display();
    Q3.Display();
    Q4.Display();

    cout<<"cooling"<<endl;

    //-------Valorisation PAC----//

/*
    FluidCore* fluidin;

    fluidin=new PRfluid("R134a");

    double eta_com=0.9;
    THeatpump PAC(fluidin,eta_com, 5, 5);


    double ThotWater=55+273.15;
    for(int j(0);j<nb_Q;j++){
        COPWater_s.clear();
        for (int i(0);i<Tex.T.size();i++){
            COPWater_s.push_back(PAC.COP_real(ThotWater,Twater[j][i],Tex.T[i]));
        }
        COPWater.push_back(COPWater_s);
    }
    SetOutput2(COPWater,(char*)file_out7);


    cout<<"PAC"<<endl;
    */
    //---------------PRODUCTION ELEC -----------------///
    solar PV;
    vector<double>Ev(8760);
    Ev=PV.ElecinTime(1,Tex);

    SetOutput( Ev,(char*)file_out6);


    cout<<"Well done!"<<endl;
    return 0;
}



void input_analyse(vector<vector<string>> table,int counttype[5]){
    counttype[0]=0;
    counttype[1]=0;
    counttype[2]=0;
    counttype[3]=0;
    counttype[4]=0;

    for(int i=0;i<table.size();i++){
        if(table[i][4]=="env")
            counttype[0]++;
        else if(table[i][4]=="con")
            counttype[1]++;
        else if(table[i][4]=="sou")
            counttype[2]++;
        else if(table[i][4]=="power")
            counttype[3]++;
        else if(table[i][4]=="state")
            counttype[4]++;
        else
            cout<<"Not find"<<endl;
    }

}
