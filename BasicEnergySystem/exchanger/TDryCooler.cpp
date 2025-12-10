//
// Created by robin.lemaire on 30.07.2021.
//

#include "TDryCooler.h"


/****************************************************************************/

void TDryCooler::kCalculation(){
    convExchCoeffCalculation();
    k = 1 / (1/ng/hE+Se/hI/Si);
}

/****************************************************************************/

void TDryCooler::convExchCoeffCalculation() {
    Se=M_PI*DiamE*DiamE/4;
    Si=M_PI*DiamI*DiamI/4;

    double Lc=HotChange.tStategaz.h-HotChange.tStateliquid.h;;  // Lc chaleur latente de cond
    hI=0.555*pow((HotChange.tStateliquid.u/HotChange.tStateliquid.v*(HotChange.tStateliquid.u/HotChange.tStateliquid.v-HotChange.tStategaz.u/HotChange.tStategaz.v)
            *9.81*pow(HotChange.pSatState.ThermCond,3)*Lc/(HotChange.tStategaz.u*(HotChange.pSatState.GetT()-Hot.pStateWall->GetT())*DiamI)), 1.0/4.0);

    double alpha=Cold.pStateIn->ThermCond/(Cold.pStateIn->ViscDyn/Cold.pStateIn->ViscCin*Cold.SpecHeat);
    double Pa=2.71*Beta*(Hot.pStateWall->GetT()-Cold.pStateIn->T)*9.81/pow(alpha*Hot.pStateIn->ViscDyn, 1.0/4.0);

    double Gr=9.81*Beta*Hot.pStateIn->T-Hot.pStateOut->T*pow(Lc,3)/Hot.pStateIn->ViscCin;

    double Nu=0.201*pow((Gr*Hot.pStateIn->Prandtl*Pa/DiamI), 1.0/3.0);
    hE=Nu*Hot.pStateIn->ThermCond/DiamE;

    double a=B/A;
    double b=1.28*A/Rex*sqrt(a-0.2);
    double la=Rex*(b-1);
    double Eps=2*hE/Cold.pStateIn->ThermCond*eAil;
    double EpsAil=tanh(Eps*la*(1+0.35*log(b)))/(Eps*la*(1+0.35*log(b)));
    double SAil=2*n*(4*A*B-3.1415*DiamI*DiamI/4); //surface d'ailettes par mètre de longueur
    double SNet=Se*(1-n*eAil); //surface nette des tubes par mètre de longueur
    //double Stot=SNet+SAil; //surface d'échange totale des tubes par mètre de longueur

    ng=1-(1-EpsAil)*SAil/(SNet+SAil);
}

/****************************************************************************/

void TDryCooler::init()
{
    TExchanger::init();
    DiamI=a;
    DiamE=b;
    A=c;
    B=d;
    Rex=e;
    eAil=f;
    n=g;
    Beta=h;

    Cold.SpecHeat=1005; //Chaleur spécifique de l'air à 300 K (27 °C), ne varie pas significativement dans les domaines de T utilisés
}

/****************************************************************************/