#ifndef CONDUITE_LIN2_H
#define CONDUITE_LIN2_H
#include "Hydraulic.h"


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
//const double Pi = 3.14159265359;
//const double g = 9.81;

class HPipe: public Hydraulic
{
    public:
        HPipe(const char *iPipeName, TInputFile &InputFile);
        HPipe();
        HPipe(Flow* iflow);
        HPipe(Flow_physic* iflow);
        ~HPipe();

        void init ();
        BasicES&  performance (); //sort tous les résulats
        float        objective(TSequence& seq);


//Pascal
    void ThermoPhysicCalculation ();
    static double MoodyFactor (double Reynolds, double roughness=0.0);
    double MoodyFactor ();
    double Diameter ();

    void getdataBase(string& Type,string & DN); // prends les donnees d'une base de donnees en fonction du type !!! a refaire car ouvre et referme le fichier à chaque pipe !!!
    void getdataBase(TSeqsArray& Records, string& DN);

    void InitConduiteData(double &iLength, double &iDiam, double &iRuAbs) ;
    void InitConduiteData(double &iMFlow, double &iLength, double &iDiam, double &iRuAbs, double &iDeltaP) ;
    void Initflow(FluidCore*&Fluid, double Tin=300, double Pin=1e5,double MFlow=0);
    void BendData (double& F_Moody);
    double Debit ();
    virtual void PressureDropCalculation ();
    void updateFlow_P();
    void DiameterCalculation ();
    double Heattransfer_mat(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i);
    double Heattransfer_mat(double &Diam,double &Diam_ext,double &Diam_ext_insu, double &Diam_ext_mant,
                            double &U_p,double &U_i, double &U_m);
    double Nusselt_lam(double Reynolds);
    double Nusselt_turb(double Reynolds);
    virtual double Heattransfer_conv();
    double Heattransfer_total(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i);
    double Heattransfer_total(double &Diam,double &Diam_ext,double &Diam_ext_insu, double &Diam_ext_mant,
                              double &U_p,double &U_i, double &U_m);
    double Heattransfer_total(); // expects known U_pipe
    void tauCalculation (double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i);
    void tauCalculation ();      // expects known U_pipe

    void HeatlossCalculation(double Tambient);
    void updateFlow_T();

    void LossCalculation(double Tambient);

    void Display();

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
    double getT();
    double getP();

    //Setters
    void setLength( double& length);
    void setDiameter( double& diameter);
    void setRugAbs( double& rugabs);
    void setMassflow( double massflow);
    void setViscDyn( double& viscdyn);
    void setViscCin( double& visccin);
    void setConduct( double& conduct);
    void setRo( double& ro);


    //public:
    int PipeNumber;
    double Length, Diam, RugAbs, RugRel=0.0, U_pipe=1.e100, U_total; //VolumeFlow==ancien Qv
    // A big value of U_pipe means: pipe material not considered
    double tau, DeltaPr_m, DeltaPr, DeltaT;     // DeltaPr_m is pressure loss per meter
    double Re, Speed, F_Moody;


    TPhysicState *pStateIn;
    Flow* flow;
    bool flow_created=false, state_created=false;
};

#endif
