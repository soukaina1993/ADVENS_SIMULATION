/*******************************************************************************

Expander.h

Performance calculations of volumetric turbine
(MacroexpanderLib.p extension on c++)

Copyright � 2000 LENI-DGM-EPFL.  All rights reserved.

$Id: Expander.cpp,v 1.0 2000/03/14 17:10:00 Malick Kane Exp $

*******************************************************************************/



#ifndef EXPANDER_H
#define EXPANDER_H

//#include<stdio.h>
//#include"Array.h"
//#include<iostream>
#include"turbine.h"
#include"../../fluids/LeeKessler.h"
#include "../../state/TThermoState.h"
#include "../../state/TPhysicState.h"

float rendement  (TSequence& s);
float puissance  (TSequence& s);



//==============================================================================

//void copy (ExpDataRcd * a, ExpDataRcd * b );
/*ExpDataRcd & operator = (ExpDataRcd & a, ExpDataRcd & b ){
   copy(a,b); return a;
}*/


//pascal funtions in MacroExpanderLib.p  library
//==============================================================================

//TODO traduire ces fonctions pascales



//==============================================================================				 double PowerOrder, int Opt);

class TExpander;

/****************************************************************************/
void copy (TExpander & dest, TExpander& src );

class expParams: public ioData{
public:

    expParams (char * ifname, unsigned int nTln=1, unsigned int nTcol=1):
            ioData(ifname, "--", nTln, nTcol)  {}
    expParams (TInputFile & ifile): ioData(ifile) {}

    virtual ~expParams(void){}

    //Get value
    double Vs	(const char * expname)   {return value(expname , "Vs")     ;}
    double VRi	(const char * expname)   {return value(expname , "VRi")    ;}
    double Power	(const char * expname)   {return value(expname , "Power")  ;}
    double EtaTs	(const char * expname)   {return value(expname , "EtaExp") ;}
    double EtaElec	(const char * expname)   {return value(expname , "EtaElec") ;}
    double Dadm	(const char * expname)   {return value(expname , "Dadm")   ;}
    double Length	(const char * expname)   {return value(expname , "Length") ;}
    double Width	(const char * expname)   {return value(expname , "Width")  ;}
    double Height	(const char * expname)   {return value(expname , "Height") ;}
    double Weight	(const char * expname)   {return value(expname , "Weight") ;}
    double NrotNom	(const char * expname)   {return value(expname , "NrotNom");}
    double Cfric1	(const char * expname)   {return value(expname , "Cfric1") ;}
    double Cfric2	(const char * expname)   {return value(expname , "Cfric2") ;}


    //-------------------------------------------------
    void RegisterParams (TExpander & expdata, const char * expname);
    void DefaultParams  (TExpander & expdata);

public:
    //ExpDataRcd * expdata;

};
/****************************************************************************/
class TExpander:public TTurbine {
public:  	
    enum expander {ANYTYPE=0, TRANE=1, MANEUROPE=2, COPELAND=3, ATLAS=4, HITACHI=5};
 
    TExpander (const TExpander  & orig);    
    //TExpander (const expander tp, const TTurbine::SpeedMode sp=TTurbine::FIXED);    	      
    TExpander (LeeKessler* &fluidin1,const TExpander::expander tp	=TExpander::ANYTYPE,const TTurbine::SpeedMode sp	=TTurbine::FIXED	);
    TExpander (LeeKessler::LKFluidRec * pflddata	=NULL,const TExpander::expander tp	=TExpander::ANYTYPE,const TTurbine::SpeedMode sp	=TTurbine::FIXED);
    TExpander(const char *iExpanderName,const TExpander::expander tp	=TExpander::ANYTYPE,const TTurbine::SpeedMode sp	=TTurbine::FIXED);

    
    virtual ~TExpander (void) 		{delete ivbles;}
    TExpander * clone  (void) 		{return new TExpander(*this);}
    void copy 	     (const TExpander & orig);
    
   
    //Set/Get expander 
    void expregister ();
    void expselector ();
    double vri 	    (void)		{return VRi;}
    double vri 	    (double VR)		{return VRi=VR;}
    double nExpander  (void)		{return NTins;}	//Ins exp in parallel
    double nExpander  (double Nins)		{return NTins=Nins;}//Ins. exp in parallel
    double nUsedExp   (void)		{return NTfon;}
    double nUsedExp   (double Nfon)		{return NTfon=Nfon;}
    double nSpeed	    (void)		{return Nvar;}
    double nSpeed	    (double Nv)		{return Nvar=Nv;}

    double Section(double Diam);
    double Diameter_cal(double Vs, double Nr);

    //ExpDataRcd & selectexp (expander AType);

    expander type (const expander AType)	{idxtp	= (unsigned int) AType; return exptp = AType;}
    expander type () 			{return exptp;  }


    //LeeKessler::LKFluidRec * fluid   ()		{return fluidinside;}
    //LeeKessler::LKFluidRec * fluid   (LeeKessler::LKFluidRec & fld)	{SetFluid (fld); return pfluid = &fld; }

    void setFluid(LeeKessler* fluidin1){fluidinside=fluidin1;}

    //Set/Get variables
    TSequence & variable (void);
    TSequence & variable (TSequence& seq);
  
   //overloading methods    
    void 	init  	   (); 			//initialize expander   
    BasicES   & 	performance ();			//expander performance 
    void    	display 	   ();			//expander display       
    float power 		   (TSequence& seq);	//electric power objective
    float efficiency   	   (TSequence& seq);	//global efficiency objective
    float cost	   	   (TSequence& seq);	//inverstment cost objective
    float objective	   (TSequence& seq);	//Select the appropiate OF
    
    //Get objective value
    float power 	   	   () 	{return ElecPower;}
    float efficiency  	   () 	{return EtaG;}
    float cost	   	   () 	{return ElecPower;}
    
    //Set/Get the default optimizer
  
    gao * optimizer (gao & theGA);

    void  optimize  (TTurbine::Criteria ob=TTurbine::POWER, gao::type GA=gao::simple);
    void SetFluid(LeeKessler::LKFluidRec & fld);


    // void InitMacroExpanderLib(); // TODO find this function
    void SimpleExpander ();
    void InitExpanderData ( float Pin,float Tin,float  Xin,float Pout, float Mflow,float  Nrot);
    void MacroExpanderM ();
    void MacroExpanderN ();
    void InitExpanderParameters();

    void VolumetricExpander();
    void ExpanderCapacity( double Mflow);


    double racine_moodyDiag(double borneInf, double borneSup, double ecart);
    double racine_suction(double borneInf, double borneSup, double ecart);
    double suction(double Pc);

    double moodyDiag(double moodyFactor);

    double mechanicalLoss(double ANrot, double NrotNom, double Cfric1, double Cfric2, double Nr);
    double mechanicalLoss(double ANrot);
    double leakage(double deltaP, double ro, double viscDyn);
    void ThermoPhysicCalculation(TPhysicState *pstate, TState *State);

    double Speed, Re, MoodyFactor,DeltaP;
    TPhysicState *pstateIn=new TPhysicState;
    TPhysicState *pstateOut=new TPhysicState;
public:
    TThermoState downstate ;
    TThermoState upstate ;


    expander  	      	exptp;			//expander type

    LeeKessler* fluidinside;
  //  LeeKessler::LKFluidRec * 		pfluid; 			//fluid data pointer
    static gao *	      	ga;			//default ga optimizer
    expParams * 	params= new expParams("expparam.txt")  ;			//expander parameters


    double Cin = 14;
protected:
    bool getflag () {return flag;}
private:

    bool flag;			//flag to distroy pdata memory or not



public:
    double  Length, Width, Height, Weight,
            Cfric1, Cfric2,
            NrotNom,
//external init input variables
    InitMassFlow, InitNrot,
//internal input variables
    Aadm, Diameter,
            Vsuction, VRi, PRi;
    double	NTins=1, NTfon=1, Nvar;
//external input variables
    double	alphaH, dTsur,
//external output variables}

            PR, CR,
//internal output variables}
    InternalHeatLoss,
            LeakageMassFlow, TurbineMassFlow,
//Specific parameters
    Mvol, Cone, mhu,
//Polytropic parameters
    EtaP, ZetaQ, Eta, Kapa, Emdh, Ework, Epoly, Np;



};




//==============================================================================
#endif	// EXPANDER_H