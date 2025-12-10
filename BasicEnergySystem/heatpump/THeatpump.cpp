//
// Created by lucile.schulthe on 31.03.2022.
//

#include "THeatpump.h"

// pointer to existing fluid (e.g. if all HP have same fluid)
THeatpump::THeatpump(FluidCore* &fluidin1, double ieta_comp, double ipinchCond, double ipinchEvap,
                     THeatpump::Type typ, THeatpump::Heattype heat) :
            fluidin(fluidin1), eta_comp(ieta_comp), pinchCond(ipinchCond), pinchEvap(ipinchEvap), typ(typ), heat(heat) {
    cond = new FlowPhaseChange(fluidin);
    evap = new FlowPhaseChange(fluidin);
}

// creation of new fluid
THeatpump::THeatpump(string fluidin1, double ieta_comp, double ipinchCond, double ipinchEvap,
                     THeatpump::Type typ, THeatpump::Heattype heat) :
            eta_comp(ieta_comp), pinchCond(ipinchCond), pinchEvap(ipinchEvap), typ(typ), heat(heat) {
    fluid_define(fluidin,fluidin1);    //creates a PRfluid --> call destructor !
    Fluid_created=true;
    cond = new FlowPhaseChange(fluidin);
    evap = new FlowPhaseChange(fluidin);
}

THeatpump::~THeatpump(){
    delete cond;
    delete evap;
    if (Fluid_created) {
        delete fluidin;
    }
}


void THeatpump::init (){}

double THeatpump::COP_max(const double & Thot,const double &Tcold){
    if (Thot <= Tcold) {
        cerr << "COP_max cannot be computed: Thot <= Tcold" << endl;
        exit(1);
    }
    return Thot / (Thot-Tcold);
}

void THeatpump::COP_real(double &Thot, double &Tcold, double &Text){
    getexergy_old(Thot, Tcold, Text);
    COP = eta_ex * COP_max(Thot,Tcold);
}

void THeatpump::COP_real(double &Thot, double &Tcold){
    if (cond != nullptr && evap != nullptr) {
        double Tevap = Tcold - pinchEvap;
        double Tcond = Thot + pinchCond;
        compute_cycle(Tcond, Tevap);
        COP = -cond->getDeltah() / (cond->StateIn->h - evap->StateOut->h);
        if (COP <= 0) {
            cerr << "COP = " << COP << endl;
            exit(1);
        }
    }
    else {
        cerr << "COP_real cannot be computed" << endl;
        exit(1);
    }
}

void THeatpump::getexergy_old(double &Thot, double& Tcold, double& Text) {
    double Tevap = Tcold - pinchEvap;
    double Tcond = Thot + pinchCond;
    if (Thot==0 || Tcold==0 || Tevap==0) {
        cerr << "Temperature is zero. Division by zero" << endl;
        exit(1);
    }
    double dS_cp = DeltaS_Cp(Tcond, Tevap);
    eta_ex = 1 - (Text / Thot) * (1 - (eta_comp * dS_cp * (Thot / Tcold - 1)) / (Tcond / Tevap - 1));
    if (eta_ex < 0 || eta_ex > 1){
        cerr << "something wrong: eta_ex = " << eta_ex << endl;
        exit(1);
    }
}

void THeatpump::getexergy_old(double Text){
    getexergy_old(Hot->StateOut->T, Cold->StateIn->T, Text);
}


/* void THeatpump::getexergy(double Text){
    double Tevap=Cold->StateIn->T-pinchEvap;
    double Tcond=Hot->StateOut->T+pinchCond;
    double dS_cp = DeltaS_Cp(Tcond,Tevap);
    eta_ex =  1;
    cerr << "getexergy not implemented, eta_ex=1" << endl;
    //eta_0r=1-(T_evap/T2_log*(1-eta_ks)-eta_ks*dhv/dhks);
} */


double THeatpump::getCost(){
    return Cost;
}

void THeatpump::compute_cycle(double &Tcond, double &Tevap){
    if (cond != nullptr && evap != nullptr) {
        //cond->Tsat(Tcond);   //TODO not needed?
        //evap->Tsat(Tevap);
        evap->StateOut->x = 1;

        cond->TSatCalculation(Tcond);    //calcul des états thermo de saturation
        evap->TSatCalculation(Tevap);    //calcul des états thermo de saturation

        //TODO sous-refroidissement, surchauffe??
        cond->StateOut->SetState(cond->tStateliquid);  //copy l'etat liquid pour la sortie du condenseur
        evap->StateOut->SetState(evap->tStategaz);     //copy l'etat gazeux pour la sortie de l'evap

        evap->StateIn->h = cond->StateOut->h;     // vanne isenthalpe
        evap->StateIn->x = (evap->StateIn->h - evap->tStateliquid.h) /
                           (evap->tStategaz.h - evap->tStateliquid.h);  // qualité du fluide en 4
        evap->StateIn->s =
                evap->StateIn->x * (evap->tStategaz.s - evap->tStateliquid.s) + evap->tStateliquid.s; // entropie en 4

        evap->StateIn->P = evap->StateOut->P;     // isobar
        cond->StateIn->P = cond->StateOut->P;

        double h2s = fluidin->h_Ps(cond->StateIn->P, evap->StateOut->s); // point isentropique
        cond->StateIn->h = evap->StateOut->h + (h2s - evap->StateOut->h) / eta_comp;
    }
}

double THeatpump::DeltaS_Cp(double &Tcond, double &Tevap){
    if (cond != nullptr && evap != nullptr) {
        compute_cycle(Tcond, Tevap);
        fluidin->Prop_Ph_g(cond->StateIn->P, cond->StateIn->h, cond->StateIn->T, cond->StateIn->v,
                           cond->StateIn->u, cond->StateIn->s);

        double Cv, Cp_cond, Cp_evap;
        fluidin->CvCp_vT(cond->StateIn->v, cond->StateIn->T, Cv, Cp_cond);
        fluidin->CvCp_vT(evap->StateOut->v, evap->StateOut->T, Cv, Cp_evap);

        double cp = (Cp_evap + Cp_cond) / 2;
        return (evap->StateOut->s - evap->StateIn->s) / cp;
    }
    else {
        cerr << "DeltaS_Cp cannot be computed" << endl;
        return 0;
    }
}


/*void THeatpump::getMassFlow(){
    if (!(Hot->StateOut->T==0 || Cold->StateIn->T==0)) {
        if(COP==0)
            COP_real(Hot->StateOut->T,Cold->StateIn->T);

        if (*Hot->Massflow==0)
            *Hot->Massflow = COP/(COP-1) * Cold->getPower();

        if (*Cold->Massflow==0)
            *Cold->Massflow = (COP-1)/COP * Hot->getPower();
    }
    else {
        cerr << "Températures non définies, impossible de calculer le débit"<<endl;
        exit(1);
    }
}*/



BasicES & THeatpump::performance(){}
float THeatpump::objective(TSequence& seq){return 0;}
void THeatpump::display(){
    cout <<"/****************************************************************************/" << endl;
    cout<< endl;
    cout<<"                         Heat pump                           "<<endl;
    if(heat==0)
        cout << "Heat power                     [kW] :       " << Hot->getPower() / 1E3 << endl;
    else if(heat==1)
        cout << "Cool power                     [kW] :       " << Cold->getPower() / 1E3 << endl;
    else{
        cout << "Heat power                     [kW] :       " << Hot->getPower() / 1E3 << endl;
        cout << "Cool power                     [kW] :       " << Cold->getPower() / 1E3 << endl;
    }
    cout<<"Cost                          [CHF] :       "<<Cost<<endl;
    cout<<"COP                      [-] :       "<<COP<<endl;
    cout<<endl;

    cout<<"                         Hot fluid "<<endl;
    cout<<"Fluid working pressure        [bar] :       "<<Hot->StateIn->P/1E5<<endl;
    cout<<"Inlet Temperature               [C] :       "<<Hot->StateIn->T -273.15<<endl;
    cout<<"Outlet Temperature              [C] :       "<<Hot->StateOut->T-273.15<<endl;
    cout<<"Mass flow rate               [kg/s] :       "<<Hot->Massflow<<endl;
    cout<<"Specific heat              [J/kg/K] :       "<<Hot->SpecHeat<<endl;
    cout<<endl;

    cout<<"                         Cold fluid"<<endl;
    cout<<"Fluid working pressure        [bar] :       "<<Cold->StateIn->P/1E5<<endl;
    cout<<"Inlet Temperature               [C] :       "<<Cold->StateIn->T -273.15<<endl;
    cout<<"Outlet Temperature              [C] :       "<<Cold->StateOut->T-273.15<<endl;
    cout<<"Mass flow rate               [kg/s] :       "<<Cold->Massflow<<endl;
    cout<<"Specific heat              [J/kg/K] :       "<<Cold->SpecHeat<<endl;
    cout<<endl;
    cout <<"/****************************************************************************/" << endl;
    cout << endl;
}


