#include <iostream>
#include "hydraulic/HPipe.h"
#include "fluids/PengRobinson.h"
#include "iofile/iofile.h"
#include "hydraulic/EPANETValve.h"
#include "state/Tstate.h"

using namespace std;

//============================================================================================================

//============================================================================================================
int main()
{
    //Input file
    const char sep = '\t';
    const int nbTitleLines = 2; //nombre de lignes avec des titres
    const int nbTitleColumns = 1; //nombre de colonnes avec des titres ie sans chiffre ie
    const char* file_inFLUID = "dataBase/DataFluids.txt";
    const char* file_inPIPE = "hydraulic/PipeData.txt";
    cout << "Test Input Ouput file : " << file_inFLUID << " "<< file_inPIPE << endl << endl;
    TInputFile  InputFileFLUID((char*)file_inFLUID, nbTitleLines, nbTitleColumns, sep);
    TInputFile  InputFilePIPE((char*)file_inPIPE, nbTitleLines, nbTitleColumns, sep);

    //FLUID
    PRfluid fluid("H2O",file_inFLUID);
    double P(0),cv(0),cp(0),u(0),h(0),s(0),sv(0),ho(0),x(0),vl(0),vg(0),Psat(0); //PEU IMPORTE LES PARAMETRES QUE TU VAS DONNER POUR CALCULER
    double ViscDyn(0),ViscCin(0),Ro(0),ThermCond(0),Prandtl(0);
    double v= 0.0006;
    double T= 300;
    fluid.Prop_vT_a(v,T,P,u,h,s,x);
    sv = fluid.sv_vT(v,T);
    if (x == 1)
    {
        fluid.CvCp_vT_g(v,T,cv,cp);
        fluid.TransProp_PT_g(P,T,ViscDyn,ViscCin,Ro,ThermCond,Prandtl);
    }
    else if (x == 0)
    {
        fluid.CvCp_vT_l(v,T,cv,cp);
        fluid.TransProp_PT_l(P,T,ViscDyn,ViscCin,Ro,ThermCond,Prandtl);
    }
    else
    {
        cout << "on est dans le domaine biphasique !" << endl;
        cout << "vl [m3/kg] : " << vl   << endl;
        cout << "vg [m3/kg] : " << vg   << endl;

    }
    cout << Ro << endl;
    cout << ViscDyn << endl;

    //PIPE
    //Déclaration de la conduite
    HPipe pipe("C1", InputFilePIPE);
    PRfluid fluid2("H2O",file_inFLUID);
    TState state;
    state.SetT(280);
    state.SetP(5E6);
    pipe.ThermoPhysicCalculation(state,fluid2);
    cout << pipe.getLength() << endl;
    cout << pipe.getDiameter() << endl;
    cout << pipe.getViscCin() << endl;
    cout << pipe.getRo() << endl;
    cout << pipe.getConduct() << endl;
    cout << pipe.getPrandtl() << endl;
    cout << pipe.getF_Moody() << endl;
    cout << ""<< endl;
    //TEST Valve
    double diameter(0.9); //ft ou 300 mm
    double volumeflowrate(1.4); //cfs ou 0,04 m3/s
    EPANETValve valve (Valve::ValveType::PRV,"V1", diameter,volumeflowrate);
    cout << valve.getFcvHeadLoss(volumeflowrate)<< endl;

    EPANETValve valve2 (Valve::ValveType::PBV,"V2", diameter,volumeflowrate);
    cout << valve2.getPbvHeadLoss(volumeflowrate)<< endl; //ft
    cout << valve2.getFlow()<< endl;
    cout << valve2.getType()<< endl;
    cout << valve2.getDiameter()<< endl;
    return 0;
}
