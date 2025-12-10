#ifndef CONDUITE_LIN_H
#define CONDUITE_LIN_H
#include "Hydraulic.h"
#include "../iofile/iofile.h"
#include "../state/TState.h"
#include "../fluids/PengRobinson.h"
#include "../state/TPhysicState.h"
#include "../state/TThermoState.h"
#include "../Flow/Flow.h"

/*
===============================================================================================================
DESCRIPTION
Dans ce code on applique l'ensemble des fonctions que l'on va utiliser pour faire les calculs sur les pertes de
charges, pertes thermiques, ... pour les conduites de sections circulaires lin�aires

PDC = Pertes de charge
PT = Pertes thermiques
lin = lin�aire
===============================================================================================================
*/


class HPipe: public Hydraulic
{
    public:
        HPipe(){};
        HPipe(const char *iPipeName, TInputFile &InputFile);
        //BasicES
       void 	      init ();
       BasicES&  performance (); //sort tous les résulats
       float        objective(TSequence& seq);


//Pascal
    void ThermoPhysicCalculation (const TState & Point, PRfluid& fluid);
    double MoodyFactor ();
    double Diameter ();

    void getdataBase(string& Type,string & DN); // prends les donnes d'une base de donnees en fonction du type !!! a refaire car ouvre et referme le fichier à chaque pipe !!!
    void getdataBase(TSeqsArray& Records, string& DN);

    void InitConduiteData (const double& iMFlow,const double& iLenght,const double& iDiam,const double& iRuAbs,const double& iDeltaP);
    void BendData (const double& F_Moody);
    double Debit ();
    void PressureDropCalculation (const TState & Point, PRfluid& fluid);
    void DiameterCalculation (const TState & Point, PRfluid& fluid);
    double Heattransfert_mat(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i);
    //Getters
    double getLength();
    double getDiameter();
    double getRugAbs();
    double getDeltaPr();
    double getMassflow();
    double getRe();
    double getPrandtl();
    double getSpeed();
    double getViscDyn();
    double getViscCin();
    double getConduct();
    double getRo();
    double getF_Moody();

    //Setters
    void  setLength(const double& length);
    void  setDiameter(const double& diameter);
    void setRugAbs(const double& rugabs);
    void setMassflow(const double& massflow);
    void setViscDyn(const double& viscdyn);
    void setViscCin(const double& visccin);
    void setConduct(const double& conduct);
    void setRo(const double& ro);


    //public:
    int PipeNumber;
    double Length, Diam, RugAbs, tau, DeltaPr, MassFlow,U_pipe; //VolumeFlow==ancien Qv
    double Re, Prandtl, Speed;
    double ViscDyn, ViscCin, Conduct, Ro, F_Moody, DiamHydro;

// add for network from Matlab
    Flow flow; // peut être une référence?
};

#endif
