//
// Created by robin.lemaire on 28.07.2021.
//

#include "TTubular.h"

/****************************************************************************/

TTubular::TTubular(TExchanger::Exchanger tp, TExchanger::Flowtype flw, TExchanger::Type typ,  TTubular::fluidIn fld)
{
    flow (flw);
    type (typ);
    fluidI (fld);

    Pinch	 = 0;
    Cost =0;
    Q = 0;
    Area = 0;
    k = 0;
    Cr = 0;
    NUT = 0;
    Efficiency=0;

    *Hot.Massflow=0;
    Hot.pStateIn->T=0;
    Hot.pStateOut->T=0;
    Hot.pStateIn->x=0;
    Hot.pStateOut->x=0;
    Hot.SpecHeat=0;

    Hot.pStateIn->P=0;

    *Cold.Massflow=0;
    Cold.pStateIn->T=0;
    Cold.pStateOut->T=0;
    Cold.pStateIn->x=0;
    Cold.pStateOut->x=0;
    Cold.SpecHeat=0;
    Cold.StateIn->P=0;
}

/****************************************************************************/

void TTubular::kCalculation(){
    convExchCoeffCalculation();
    k = 1 / (1/hI*DiamE/DiamI+1/hE);
}

/****************************************************************************/

void TTubular::convExchCoeffCalculation() {
	
	//// Issu du livre de these fourni par M. Kane
	
	
    double Re;
    double Nu;
    double Gi;  //debit massique surfacique intérieur
    Gi = *Hot.Massflow / (3.14 * (DiamI / 2) * (DiamI / 2));
    double Ge;  //debit massique surfacique extérieur
    Ge = *Cold.Massflow / (3.14 * (DiamH / 2) * (DiamH / 2));



    if (fld==0)//// si le fluide chaud est a l interieur des tubes
    {
        if (typ == 0)     // sans changement de phase
        {
            //hI
            Re = Gi * DiamI / Hot.pStateIn->ViscDyn;
            Nu = 0.023 * pow(Re, 0.8) * pow(Hot.pStateIn->Prandtl, 1.0 / 3.0);
            hI = Nu * Hot.pStateIn->ThermCond / DiamI;

            //hE
            Re = Ge * DiamH / Cold.pStateIn->ViscDyn;
            Nu = 0.36 * pow(Re, 0.55) * pow(Cold.pStateIn->Prandtl, 1.0 / 3.0);
            hE = Nu * Cold.pStateIn->ThermCond / DiamH;
        } else if (typ == 1) //condenseur
        {
            //hI

            double Prh=HotChange.tStateliquid.u*Hot.SpecHeat/HotChange.pSatState.ThermCond;  //Prandt number for liquid, hot
            double Prc=ColdChange.tStateliquid.u*Cold.SpecHeat/ColdChange.pSatState.ThermCond ; //Prandt number for liquid, cold

            double Xtt;
            Xtt = pow((HotChange.tStategaz.u/HotChange.tStateliquid.u), 0.1) *
                  pow(((1 - Hot.pStateIn->x) / Hot.pStateIn->x), 0.9) *
                  pow((HotChange.tStategaz.u/HotChange.tStategaz.v)/(HotChange.tStateliquid.u/HotChange.tStateliquid.v), 0.5);
            Re = Gi * (1 - Hot.pStateIn->GetX()) * DiamI / HotChange.tStateliquid.u;
            Nu = Prh * pow(Re, 0.9) * pow(0.15 * (pow(Xtt, -1) + 2.85 * pow(Xtt, -0.476)), 1.15) /
                 (5 * Prh + 5 * log (1 + 5 * Prh) +
                 2.5 * log (0.0031 * pow(Re, 0.812)));
            hI = Nu * HotChange.pSatState.ThermCond / DiamI;



            //hE
            Re = Ge * (1 - Cold.pStateIn->GetX()) * DiamH / ColdChange.tStateliquid.u;
            Nu = 0.3 * pow(Re, 0.6) * pow(Prc, 0.4) *
                 sqrt((ColdChange.tStateliquid.u/ColdChange.tStateliquid.v)/(ColdChange.tStategaz.u/ColdChange.tStategaz.v) + 1);
            hE = Nu * ColdChange.pSatState.ThermCond/ DiamH;
        } else if (typ == 2) //evap
        {
            //hI
            double Prh=HotChange.tStateliquid.u*Hot.SpecHeat/HotChange.pSatState.ThermCond;  //Prandt number for liquid, hot
            double Prc=ColdChange.tStateliquid.u*Cold.SpecHeat/ColdChange.pSatState.ThermCond ; //Prandt number for liquid, cold


            Re = Gi * DiamI * 0.5 / HotChange.tStateliquid.u;
            cout<<"DiamI= "<<DiamI;
            cout<<"Gi= "<<Gi;
            cout<<"Re= "<<Re;
            Nu = 0.06 * pow((HotChange.tStateliquid.u/HotChange.tStateliquid.v)/(HotChange.tStategaz.u/HotChange.tStategaz.v), 0.28) * pow(Re, 0.85) *
                 pow(Prh, 0.4);
            hI = Nu * HotChange.pSatState.ThermCond / DiamI;

            //hE
            Re = Ge * (1 - Cold.pStateIn->GetX()) * DiamH / ColdChange.tStateliquid.u;
            Nu = 0.36 * pow(Re, 0.35) * pow(Prc, 1 / 3) / pow(1 - 0.5 /*alphaH*/, 0.63);
            hE = Nu * ColdChange.pSatState.ThermCond / DiamH;
        }
    }
    else if(fld==1)////si le fluide froid est à l'intérieur des tubes
    {
        if (typ == 0)     // sans changement de phase
        {
            //hI
            Re = Gi * DiamI / Cold.pStateIn->ViscDyn;
            Nu = 0.023 * pow(Re, 0.8) * pow(Cold.pStateIn->Prandtl, 1.0 / 3.0);
            hI = Nu * Cold.pStateIn->ThermCond / DiamI;

            //hE
            Re = Ge * DiamH / Hot.pStateIn->ViscDyn;
            Nu = 0.36 * pow(Re, 0.55) * pow(Hot.pStateIn->Prandtl, 1.0 / 3.0);
            hE = Nu * Hot.pStateIn->ThermCond / DiamH;
        } else if (typ == 1) //condenseur
        {
            //hI
            double Prh=HotChange.tStateliquid.u*Hot.SpecHeat/HotChange.pSatState.ThermCond;  //Prandt number for liquid, hot
            double Prc=ColdChange.tStateliquid.u*Cold.SpecHeat/ColdChange.pSatState.ThermCond ; //Prandt number for liquid, cold

            double Xtt;
            Xtt = pow((ColdChange.tStategaz.u/ColdChange.tStateliquid.u), 0.1) *
                  pow(((1 - Cold.pStateIn->x) / Cold.pStateIn->x), 0.9) *
                  pow((ColdChange.tStategaz.u/ColdChange.tStategaz.v)/(ColdChange.tStateliquid.u/ColdChange.tStateliquid.v), 0.5);
            Re = Gi * (1 - Cold.pStateIn->GetX()) * DiamI / ColdChange.tStateliquid.u;
            Nu = Prc * pow(Re, 0.9) * pow(0.15 * (pow(Xtt, -1) + 2.85 * pow(Xtt, -0.476)), 1.15) /
                 (5 * Prc + 5 * log (1 + 5 * Prc) +
                  2.5 * log (0.0031 * pow(Re, 0.812)));
            hI = Nu * ColdChange.pSatState.ThermCond / DiamI;

            //hE
            Re = Ge * (1 - Hot.pStateIn->GetX()) * DiamH / HotChange.tStateliquid.u;
            Nu = 0.3 * pow(Re, 0.6) * pow(Prh, 0.4) *
                 sqrt((HotChange.tStateliquid.u/HotChange.tStateliquid.v)/(HotChange.tStategaz.u/HotChange.tStategaz.v) + 1);
            hE = Nu * HotChange.pSatState.ThermCond/ DiamH;
        } else if (typ == 2) //evap
        {
            //hI

            double Prh=HotChange.tStateliquid.u*Hot.SpecHeat/HotChange.pSatState.ThermCond;  //Prandt number for liquid, hot
            double Prc=ColdChange.tStateliquid.u*Cold.SpecHeat/ColdChange.pSatState.ThermCond ; //Prandt number for liquid, cold

            Re = Gi * DiamI * 0.5/ ColdChange.tStateliquid.u;
            Nu = 0.06 * pow((ColdChange.tStateliquid.u/ColdChange.tStateliquid.v)/(ColdChange.tStategaz.u/ColdChange.tStategaz.v), 0.28) * pow(Re, 0.85) *
                 pow(Prc, 0.4);
            hI = Nu * ColdChange.pSatState.ThermCond / DiamI;

            //hE
            Re = Ge * (1 - Hot.pStateIn->GetX()) * DiamH / HotChange.tStateliquid.u;
            Nu = 0.36 * pow(Re, 0.35) * pow(Prh, 1 / 3) / pow(1 - 0.5 /*alphaH*/, 0.63);
            hE = Nu * HotChange.pSatState.ThermCond / DiamH;
        }
    }
}

/****************************************************************************/

void TTubular::init()
{
	TExchanger::init();
	DiamE=a;
	DiamH=b;
	DiamI=c;
}


