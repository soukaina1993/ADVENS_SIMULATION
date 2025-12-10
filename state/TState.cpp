//
// Created by Yolaine Adihou on 09/03/2021.
//

#include "TState.h"
#include <iostream>
#include "../iofile/iofile.h"


TState::TState(const double iAv, const double iAP, const double iAT, const double iAh,
               const double iAs, const double iAx) : v(iAv), P(iAP), T(iAT), h(iAh), s(iAs), x(iAx) {}

void TState::SetState(const TState& AState) {
    v = AState.v ;
    P = AState.P ;
    T = AState.T ;
    h = AState.h ;
    s = AState.s ;
    x = AState.x ;
}

void TState::SetState(TState* AState) {
    v = AState->v ;
    P = AState->P ;
    T = AState->T ;
    h = AState->h ;
    s = AState->s ;
    x = AState->x ;
}

void TState::GetState(TState& AState) {
    AState.v = v ;
    AState.P = P ;
    AState.T = T ;
    AState.h = h ;
    AState.s = s ;
    AState.x = x ;
}

void TState::GetState(TState *AState) {
    AState->v = v ;
    AState->P = P ;
    AState->T = T ;
    AState->h = h ;
    AState->s = s ;
    AState->x = x ;
}

void TState::GetState(const char* file_in, int index) {
    TInputFile  InputFile((char*)file_in, 1, 1, '\t');
    InputFile.open();
    int nbRecords = InputFile.nRecords();
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    Records = *InputFile.GetRecords();
    TSequence* PtrSeq = Get(Records, index);

    T = PtrSeq->Get(0);
    P = PtrSeq->Get(1);
    v = PtrSeq->Get(2);
    h = PtrSeq->Get(3);
    s = PtrSeq->Get(4);
    x = PtrSeq->Get(5);

}

void TState::Display() {
    std::cout << "P     [kPa] : " << P / 1E3 << std::endl;
    std::cout << "T     [degC] : " << T - 273.15 << std::endl;
    std::cout << "v     [m^3/kg]: " << v << std::endl;
    std::cout << "u     [kJ/kg] : " << u / 1E3 << std::endl;
    std::cout << "h     [kJ/kg] : " << h / 1E3 << std::endl;
    std::cout << "s     [kJ/kg K] : " << s / 1E3 << std::endl;
    std::cout << "x     [-] : " << x << std::endl;
}


void TState::State_PT(FluidCore* fluid){
    double ps=fluid->P_T_s(T);
    if(P>ps)
        fluid->Prop_PT_g(P,T,v,u,h,s);
    else if (P==ps){
        if (x==-1 &&h==0) x=0;
        fluid->Prop_Tx_s(T,x,v,ps,u,h,s);}
    else
        fluid->Prop_PT_l(P,T,v,u,h,s);
}

void TState::State_Ps(FluidCore* fluid){
    double p,stemp;
    v=fluid->v_Ps_a(P,s);
    T=fluid->T_vP(v,P);
    fluid->Prop_vT(v,T,p,u,h,stemp);
}

void TState::State_Ph(FluidCore* fluid) {
    double p, htemp;
    v = fluid->v_Ph_a(P, h);
    T = fluid->T_vP(v, P);
    fluid->Prop_vT(v, T, p, u, htemp, s);
}

void TState::StateCalculation(FluidCore* fluid) {
    if (T!=0){
        if (v!=0){
            fluid->Prop_vT(v,T,P,u,h,s);
        }
        else if(h!=0){
            double htemp;
            v=fluid->v_Th_a(T,h);
            fluid->Prop_vT(v,T,P,u,htemp,s);
        }
        else if(P!=0){
            double ps=fluid->P_T_s(T);
            if(P>ps)
                fluid->Prop_PT_g(P,T,v,u,h,s);
            else if (P==ps){
                if (x==-1 &&h==0) x=0;
                fluid->Prop_Tx_s(T,x,v,ps,u,h,s);}
            else
                fluid->Prop_PT_l(P,T,v,u,h,s);
        }
        else if(s!=0){
            double stemp;
            v=fluid->v_Ts_a(T,s);
            fluid->Prop_vT(v,T,P,u,h,stemp);
        }

        else if(x!=-1){
            fluid->Prop_Tx_s(T,x,v,P,u,h,s);
        }
        else {
            cout << "Pas encore implementer"<<endl;
        }
    }
    else if(P!=0){
        if (v!=0){
            double p;
            T=fluid->T_vP(v,P);
            fluid->Prop_vT(v,T,p,u,h,s);
        }

        else if (h!=0){
            double p,htemp;
            v=fluid->v_Ph_a(P,h);
            T=fluid->T_vP(v,P);
            fluid->Prop_vT(v,T,p,u,htemp,s);
        }
        else if(s!=0){
            double p,stemp;
            v=fluid->v_Ps_a(P,s);
            T=fluid->T_vP(v,P);
            fluid->Prop_vT(v,T,p,u,h,stemp);
        }

        else if(x!=-1){
            double p;
            T=fluid->T_P_s(P);
            fluid->Prop_Tx_s(T,x,v,p,u,h,s);
        }
        else {
            cout << "Pas encore implementer"<<endl;
        }
    }
    else
        cout <<"error, need T or P to calculate"<<endl;

}void TState::calcul_state(FluidCore* &fluid){
    init(fluid);
    fluid->Prop_vT(v,T,P,u,h,s);
}
void TState::init(FluidCore* &fluid){
    if (T!=0) {
        if (v != 0) {}
        else if (h != 0)
            v = fluid->v_Th_a(T, h);
        else if (s != 0)
            v = fluid->v_Ts_a(T, s);
        else if (P != 0) {
            double ps = fluid->P_T_s(T);
            if (P < ps)
                v = fluid->v_PT_g(P, T);
            else if (P > ps)
                v = fluid->v_PT_l(P, T);
            else
                cerr << "Not defined" << endl;
        }
        else
            v = fluid->v_Tx_a(T, x);
    }
    else if(P!=0){
        if (v!=0)
            T=fluid->T_vP(v,P);
        else if (h!=0){
            v=fluid->v_Ph_a(P,h);
            T=fluid->T_vP(v,P);
        }
        else if (s!=0){
            v=fluid->v_Ps_a(P,s);
            T=fluid->T_vP(v,P);
        }
        else if (u!=0)
            cerr << "Pas encore implementer"<<endl;
        else {
            T=fluid->T_P_s(P);
            v=fluid->v_Tx_a(T,x);}
    }
    else
        cerr <<"error, need T or P to calculate"<<endl;
}


void TState::init_all(FluidCore* &fluid){
    if (T!=0){
        if (v!=0){
            fluid->Prop_vT(v,T,P,u,h,s);
        }
        else if(P!=0){
            double ps=fluid->P_T_s(T);
            if(P<ps){
                fluid->Prop_PT_g(P,T,v,u,h,s);}
            else if (P==ps){
                if (x==-1)
                    x=0;       //TODO
                fluid->Prop_Tx_s(T,x,v,ps,u,h,s);}
            else
                fluid->Prop_PT_l(P,T,v,u,h,s);
        }
        else if(h!=0){
            double htemp;
            v=fluid->v_Th_a(T,h);
            fluid->Prop_vT(v,T,P,u,htemp,s);
        }
        else if(s!=0){
            double stemp;
            v=fluid->v_Ts_a(T,s);
            fluid->Prop_vT(v,T,P,u,h,stemp);
        }

        else if(x!=-1){
            fluid->Prop_Tx_s(T,x,v,P,u,h,s);
        }
        else {
            cerr << "Pas encore implementer"<<endl;
        }
    }
    else if(P!=0){
        if (v!=0){
            double p;
            T=fluid->T_vP(v,P);
            fluid->Prop_vT(v,T,p,u,h,s);
        }

        else if (h!=0){
            double p,htemp;
            v=fluid->v_Ph_a(P,h);
            T=fluid->T_vP(v,P);
            fluid->Prop_vT(v,T,p,u,htemp,s);
        }
        else if(s!=0){
            double p,stemp;
            v=fluid->v_Ps_a(P,s);
            T=fluid->T_vP(v,P);
            fluid->Prop_vT(v,T,p,u,h,stemp);
        }

        else if(x!=-1){
            double p;
            T=fluid->T_P_s(P);
            fluid->Prop_Tx_s(T,x,v,p,u,h,s);
        }
        else {
            cerr << "Pas encore implementer"<<endl;
        }
    }
    else
        cerr <<"error, need T or P to calculate"<<endl;
}

void TState::init_all_x(FluidCore* &fluid){
    if (T!=0){
        if (v!=0){
            fluid->Prop_vT_a(v, T,P,u,h,s,x);
        }
        else if(P!=0){
            double ps=fluid->P_T_s(T);
            if(P<ps){
                fluid->Prop_PT_g(P ,T,v,u,h,s);}
            else if (P==ps){
                if (x==-1)
                    cerr << "No unique solution in biphasic domain" << endl;
                else
                    fluid->Prop_Tx_s(T,x,v,ps,u,h,s);
            }
            else
                fluid->Prop_PT_l(P,T,v,u,h,s);
        }
        else if(h!=0){
            double htemp;
            v=fluid->v_Th_a(T,h);
            fluid->Prop_vT_a(v,T,P,u,htemp,s, x);
        }
        else if(s!=0){
            double stemp;
            v=fluid->v_Ts_a(T,s);
            fluid->Prop_vT_a(v,T,P,u,h,stemp, x);
        }
        else if(x!=-1){
            fluid->Prop_Tx_s(T,x,v,P,u,h,s);
        }
        else {
            cerr << "Pas encore implementer"<<endl;
        }
    }
    else if(P!=0){
        if (v!=0){
            double p;
            T=fluid->T_vP(v,P);
            fluid->Prop_vT_a(v,T,p,u,h,s, x);
        }
        else if (h!=0){
            double p,htemp;
            v=fluid->v_Ph_a(P,h);
            T=fluid->T_vP(v,P);
            fluid->Prop_vT_a(v,T,p,u,htemp,s, x);
        }
        else if(s!=0){
            double p,stemp;
            v=fluid->v_Ps_a(P,s);
            T=fluid->T_vP(v,P);
            fluid->Prop_vT_a(v,T,p,u,h,stemp, x);
        }
        else if(x!=-1){
            double p;
            T=fluid->T_P_s(P);
            fluid->Prop_Tx_s(T,x,v,p,u,h,s);
        }
        else {
            cerr << "Pas encore implementer"<<endl;
        }
    }
    else
        cerr <<"error, need T or P to calculate"<<endl;
}

void TState::Display(string& var1){
    if(var1=="T") cout<<T<<' ';
    else if(var1=="P") cout<<P<<' ';
    else if(var1=="x") cout<<x<<' ';
    else if(var1=="v") cout<<v<<' ';
    else if(var1=="h") cout<<h<<' ';
    else if(var1=="s") cout<<s<<' ';
    else if(var1=="u")cout<<u<<' ';
    else if(var1=="res_all"){
        cout<<T<<' '<<P<<' '<<h<<' '<<s<<' '<<v<<' '<<u<<' '<<x<<' ';
    }
    else if(var1=="all"){
        Display();}
    else{
        cerr<<"error: not T, P, x, v, h, s, u"<<endl;
    }}

void TState::assign(string &var1, double &arg1){
    if(var1=="T") T=arg1;
    else if(var1=="P") P=arg1;
    else if(var1=="x") x=arg1;
    else if(var1=="v")
        v=arg1;
    else if(var1=="h")
        h=arg1;
    else if(var1=="s")
        s=arg1;
    else if(var1=="u")
        u=arg1;
    else
        cerr<<"error: not T, P, x, v ,h, s, u"<<endl;
}