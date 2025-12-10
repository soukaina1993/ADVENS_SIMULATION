//
// Created by robin.lemaire on 21.06.2021.
//

#ifndef FLUIDS_TEXCHANGER_H
#define FLUIDS_TEXCHANGER_H

#include <iostream>
#include <string>
#include "../BasicES.h"
#include "../../iofile/iofile.h"
#include "../../hydraulic/HPipe.h"
#include "../../state/TPhysicState.h"
#include "../../state/TThermoState.h"
#include "../../fluids/FluidCore.h"
#include "../../Flow/Flow_physic.h"



#include"../../gao/gao.h"
#include "../../Flow/Flow.h"

/****************************************************************************/


class exchParams: public ioData{
public:

    exchParams (char * ifname, unsigned int nTln=1, unsigned int nTcol=1):
            ioData(ifname, "--", nTln, nTcol)  {}
    exchParams (TInputFile & ifile): ioData(ifile) {}


    //Get value
    double Getk (const char * exchname)   {return value(exchname, (char*)"k");}
    double Geta	(const char * exchname)   {return value(exchname , (char*)"a");}
    double Getb	(const char * exchname)   {return value(exchname , (char*)"b");}
    double Getc	(const char * exchname)   {return value(exchname , (char*)"c");}
    double Getd	(const char * exchname)   {return value(exchname , (char*)"d");}
    double Gete	(const char * exchname)   {return value(exchname , (char*)"e");}
    double Getf	(const char * exchname)   {return value(exchname , (char*)"f");}
    double Getg (const char * exchname)   {return value(exchname , (char*)"g");}
    double Geth (const char * exchname)   {return value(exchname , (char*)"h");}
};


class TExchanger : public BasicES	
{
public:


    enum Flowtype {PARALLEL	=0, COUNTER	=1, CROSS=2};   // 3 types de flux: co-courant, contre-courant et courants croisés
    enum Criteria  {RENDEXERG	=0, POLL	=1, COST	=2};  //3 aspects à optimiser : rendement exergétique, pollution et coût
    enum Type{DEFAULT=0, COND=1, EVAP=2, EVAPCOND=3};      // 4 types d'échangeurs : sans changement de phases, condeuseur, évaporateur, évapocondenseur
    enum Exchanger {ANYTYPE=0, P1=1, EV=2};


    explicit TExchanger(TExchanger::Exchanger tp	=TExchanger::ANYTYPE, TExchanger:: Flowtype	flw=TExchanger::PARALLEL, TExchanger::Type sp	=TExchanger::DEFAULT);
    explicit TExchanger (FluidCore* &fluidin1, TExchanger::Exchanger tp	=TExchanger::ANYTYPE, TExchanger:: Flowtype	flw=TExchanger::PARALLEL, TExchanger::Type sp	=TExchanger::DEFAULT);
    ~TExchanger() override {delete params;}


    Flowtype flow (const Flowtype FType)	{return flw = FType;}
    Flowtype flow () 			{return flw;  }

    Exchanger ex (const Exchanger FType)	{return exchtp = FType;}
    Exchanger ex () 			{return exchtp;  }

    Criteria crit(const Criteria Crit) {return obj=Crit;}
    Criteria crit ()            {return obj;  }

    Type type(const Type AType) {return typ=AType;}
    Type type ()             {return typ;}

    void copy(TExchanger &src);

    //// Fonctions de calcul
    double deltaTlm();

    double getSurface();

    double getPower();

    double getThi();

    double getTho();

    double getTci();

    double getTco();

    double getCr();

    double efficiency();

    double getPinch();

    void getSpecHeat();

    void getMassFlow();

    double getNTU();

    double getSurfaceByNTU();

    void CpCalculation();

    double efficiencyByNUT();

    void outletTemp();

    void deltaT();

    virtual void design(){}


    //// Fonctions BasicES

    void init() override;
    TSequence & variable ();
    TSequence & variable (TSequence& seq);

    BasicES & performance();
    void  optimize  (TExchanger::Criteria ob=TExchanger::RENDEXERG, gao::type GA=gao::simple);

    void display();

    float objective(TSequence& seq);


    //// Fonctions pour les classes-filles

    virtual void kCalculation() {}

    ////Set & Get
    void hexregister ();

    double pinch 	    ()		{return Pinch;}
    double pinch 	    (double dT)		{return Pinch=dT;}
    double totArea       ()		{return Area;}
    double totArea  (double A)		{return Area=A;}
    double exchCoeff ()     {return k;}
    double exchCoeff (double ik)     {return k=ik;}
    double cost    ()		{return Cost;}
    double cost (double cost)		{return Cost=cost;}
    double heat     ()		{return Q;}
    double heat  (double heat)		{return Q=heat;}
    double capaciteRatio       ()		{return Cr;}
    double capaciteRatio  (double cr)		{return Cr=cr;}
    double nut       ()		{return NUT;}
    double nut  (double n)		{return NUT=n;}
    double eff       ()		{return Efficiency;}
    double eff  (double E)		{return Efficiency=E;}
    double heatCoeff    ()		{return k;}
    double heatCoeff (double kl)		{return k=kl;}
    double dTlog        ()      {return DeltaTlog;}
    double dTlog    (double dTlog)  {return DeltaTlog=dTlog;}
    double LeftDt       ()		{return LeftDeltaT;}
    double LeftDt  (double lDt)		{return LeftDeltaT=lDt;}
    double RightDt       ()		{return RightDeltaT;}
    double RightDt  (double RDt)		{return RightDeltaT=RDt;}

    void RegisterParams(const char * exchname);
    void DefaultParams ();

        ////Donnees

    Flow_physic Hot;
    Flow_physic Cold;

    double Pinch;
    double Area;
    double k;
    double Cost, Q,
            Cr, NUT, Efficiency;
    double DeltaTlog, LeftDeltaT, RightDeltaT;

    //TODO pressure loss;

        ////Parametres

    double a, b, c, d, e, f, g, h;


    Flowtype	flw;				//flow type
    Criteria	obj;				//optimization criteria
    Type        typ;                        //Type
    Exchanger  	exchtp;			//expander type
    gao	     *	bcga;  				//basic ga optimizer
    exchParams  *params	= new exchParams((char*) "exchparam.txt")  ;			//exchanger parameters

    static gao *ga;
    FluidCore* fluidinside;

protected:
    unsigned int  idxtp;				//Turbine type index
    unsigned int  idxsp;				//speed index
    unsigned int  idxnb;				//Turbine number in parallel
};


#endif //FLUIDS_TEXCHANGER_H
