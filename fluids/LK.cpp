
// Created by lucile.schulthe on 14.06.2022.
//

#include "LK.h"
#define repeat do
#define until(exp) while(!(exp))

void LK2::LK_Sonic_PT_g(double Psup, double Tsup, double &vLaval,double &PLaval, double &TLaval, double &Csonic){
    int i=0;
    double hsup, Ssup, hl, sl;
    double DPV, DPT, CV, CP, Ci;

    //output: vLaval, hsup, Ssup, DPV, DPT, CV, CP, Csonic
    LK_Prop_PT_g(Psup, Tsup, vLaval, hsup, Ssup, DPV, DPT, CV, CP, Csonic);

    do{
        Ci = Csonic;
        hl= hsup -sqr(Csonic) / 2 / 1e3;
        iiLK_Prop_hs_g(hl, Ssup, vLaval, PLaval, TLaval);
        LK_Prop_PT_g(Psup, Tsup, vLaval, hsup, Ssup, DPV, DPT, CV, CP, Csonic);
        if(i >= 50){
            cout<<"error LK_Sonic_PT_g"<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
    }until(abs(Csonic - Ci) / Csonic < 0.001);


}

void LK2::LK_TransProp_PT_g (double P,double T, double &ViscDyn,double & ViscCin,double & Ro,double & ThermCond, double & Prandtl){
// ViscDyn in [Pa s],ViscCin in [m2/s],ThermCond in [W/(m K)]}
// from "The properties of gases and liquids" of Reid & Prausnitz, Chap. 9.4,9.6,10.5}
    k = 0;
    double v, h, s, DPV, DPT, Cv, Cp, SV, Tr, Ror, Mur, Tstar, Vcu, Mwu,
    Visc0, ViscStar, Visc2Star, OmegaV, Fc,
    G1, G2, y, Psi, Ag, Bg, Z, q;

    LK_Prop_PT_g(P, T, v, h, s, DPV, DPT, Cv, Cp, SV);
    Cv = Cv * 1000;//  ==> J/kg K}
    Cp = Cp * 1000;//  ==> J/kg K}
// *** viscosities ***}
    Ro = 1 / v;
    Ror = vc * Ro;
    Mwu = Mw * 1E-3;// [kg/mol]}
    Vcu = 1E3 * vc * Mw;// [cm3/mole]}
    Tr = T / Tc;
    Tstar = 1.2593 * Tr;
// eq. (9-4.3)}
    OmegaV = 1.16145 * XpwrY(Tstar, -0.14874) + 0.52487 * exp(-0.7732 * Tstar) + 2.16178 * exp(-2.43787 * Tstar);
// eq. (9-4.11)}
    Mur = 131.3 * DipM / sqrt(Vcu * Tc);
// eq. (9-4.10)}
    Fc = 1 - 0.2756 * w + 0.059035 * sqr(Mur) + k;
// eq. (9-6.18)}
    y = vc * Ro / 6;
// eq. (9-6.19)}
    G1 = (1 - y / 2) / XpwrI(1 - y, 3);
// eq. (9-6.20)}
    G2 = (E(1) / y * (1 - exp(-E(4) * y)) + E(2) * G1 * exp(E(5) * y) + E(3) * G1) / (E(1) * E(4) + E(2) + E(3));
// eq. (9-6.21)}
    Visc2Star = E(7) * sqr(y) * G2 * exp(E(8) + E(9) / Tstar + E(10) / sqr(Tstar));
// eq. (9-6.17)}
    ViscStar = sqrt(Tstar) / OmegaV * (Fc * (1 / G2 + E(6) * y)) + Visc2Star; // [µP]}
// eq. (9-6.16) Chung et al. method}
    ViscDyn = 1E-7 * ViscStar * 36.344 * sqrt(Mw * Tc) / XpwrY(Vcu, 2 / 3);
    ViscCin = ViscDyn / Ro;
// *** Thermal conductivity ***}
// eq. (9-4.9)}
    Visc0 = 40.785 * Fc * sqrt(Mw * T) / XpwrY(Vcu, 2 / 3) / OmegaV; // [µP]}
// eq. (10-3.17)}
    Ag = Cv * Mwu / LKconst.RBAR - 1.5;
    Bg = 0.7862 + w * (-0.7109 + 1.3168 * w);
    Z = 2 + 10.5 * sqr(Tr);
    Psi = 1 + Ag * ((0.215 + 0.28288 * Ag - 1.061 * Bg + 0.26665 * Z) / (0.6366 + Bg * Z + 1.061 * Ag * Bg));
    G2 = (B(1) / y * (1 - exp(-B(4) * y)) + B(2) * G1 * exp(B(5) * y) + B(3) * G1) / (B(1) * B(4) + B(2) + B(3));
    q = 3.586E-3 * sqrt(Tc / Mwu) / XpwrY(Vcu, 2 / 3);
    ThermCond = 31.2 * Visc0 * 1E-7 * Psi / Mwu * (1 / G2 + B(6) * y) + q * B(7) * sqr(y) * sqrt(Tr) * G2;
    Prandtl = Cp * ViscDyn / ThermCond;

}

void LK2::LK_Prop_PT_g(double P, double T,  double &V,  double &H,  double &S,  double &DPV,
                                  double &DPT,  double &CV,  double &CP, double &SV) {

    double DPVo, DPVr, DPTo, DPTr, VR, VRo, VRr, PR, TR,Eo, Er;

    PR= P/kPc;
    TR = T/Tc;
    VRro_PTr_g(PR, TR, VR, VRo, VRr, Eo, Er);

    H_VT(VRo, VRr, PR, TR, Eo, Er, H);
    S_VT(VRo, VRr, TR, Eo, Er, S);
    DPT_VTor(VRo, VRr, TR, DPTo, DPTr, DPT);
    DPV_VTor(VRo, VRr, TR, DPVo, DPVr, DPV);
    CVP_VT(VRo, VRr, Eo, Er, DPVo, DPVr, DPTo, DPTr, TR, CV, CP);
    V = VR*vc/Zc;
    SV = sqrt(-1000*CP/CV *DPV)*V;

}





void LK2::iLK_Prop_Ph_g(double P, double h, double &v, double &T, double &s) {
    double DPT, DPV, CV, CP, SV, hi, Ti;
    int i=0;
    Ts_P(P,T,DPT);
    LK_Prop_PT_g(P, T, v, hi, s, DPV, DPT, CV, CP, SV);
    do{
        Ti= T;
        T= Ti- (hi-h)/CP;
        LK_Prop_PT_g(P, T, v, hi, s, DPV, DPT, CV, CP, SV);
        if(i>50){
            cout<<"error iLK_Prop_Ph_g"<<endl;
            exit(EXIT_FAILURE);
        }
    }until((abs(hi - h) < 1E-6) & (abs(Ti - T) < 1E-4));

}

void LK2::iLK_Prop_vP_g (double v,double P,double & T, double &h, double &s) {
    double Ti, vi, dpv, dpt, cvr, cpr, sv;
    int i;
// initialisation}
    i = 0;
    Ts_P(P, T, dpt);
    LK_Prop_PT_g(P, T, vi, h, s, dpv, dpt, cvr, cpr, sv);
// iteration}
    do {
        Ti = T;
        T = Ti + (vi - v) * dpv / dpt;
        LK_Prop_PT_g(P, T, vi, h, s, dpv, dpt, cvr, cpr, sv);
        if(i>50) {cout<<"error iLK_Prop_vP_g "<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
    }until((abs(Ti - T) < 1E-3) && (abs(vi - v) < 1E-7));
// writeln(i : 0);}
}
void LK2::iLK_Prop_vT_g (double v,double T,double & P,double & h,double & s) {
    double Pi, vi, DPT, DPV, CV, CP, SV;
    int i;
// initialisation}
    i = 0;
        if(v / vc * Zc < 2)
            Ps_T(P, T);
        else
            P = LKconst.RBAR / Mw * T / v;
    LK_Prop_PT_g(P, T, vi, h, s, DPV, DPT, CV, CP, SV);
// itŽration}
    do {
        Pi = P;
        P = Pi - (vi - v) * DPV / 4;    // ajoutŽ par Nicolas le 31.1.97 : /4}
        LK_Prop_PT_g(P, T, vi, h, s, DPV, DPT, CV, CP, SV);
        if (i > 50)  {cout<<"error iLK_Prop_vT_g "<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
    }until((abs(P - Pi) < 1E-4) && (abs(v - vi) < 1E-6));
// writeln(i : 0);}
}

void LK2::iLK_Prop_Th_g (double T, double h,double & v, double &P,double & s) {
    double Pi, hi, DPT, DPV, CV, CP, SV;
    int i;
// initialisation}
    i = 0;
    Ps_T(P, T);
    LK_Prop_PT_g(P, T, v, hi, s, DPV, DPT, CV, CP, SV);
// itŽration}
    do {
        Pi = P;
        P = Pi - (hi - h) / (v + T * DPT / DPV);
        LK_Prop_PT_g(P, T, v, hi, s, DPV, DPT, CV, CP, SV);
        if(i>50){cout<<"error iLK_Prop_Th_g "<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
    }until((abs(P - Pi) < 1E-4) && (abs(h - hi) < 1E-4));
// writeln(i : 0);}
}
void LK2::iLK_Prop_Ts_g (double T,double s,double & v, double &P, double &h) {
    double Pi, si, DPV, DPT, Cv, Cp, Sv;
    int i;
// initialisation}
    i = 0;
    Ps_T(P, T);
    LK_Prop_PT_g(P, T, v, h, si, DPV, DPT, Cv, Cp, Sv);
// itŽration}
    do {
        Pi = P;
        P = Pi - (si - s) * DPV / DPT;
        LK_Prop_PT_g(P, T, v, h, si, DPV, DPT, Cv, Cp, Sv);
        if (i > 50){cout<<"error iLK_Prop_Ts_g "<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
    }until((abs(si - s) < 1E-5) && (abs(Pi - P) < 1E-4));
// writeln(i : 0);}
}




void LK2::VRro_PTr_g(double PR, double TR, double &VR, double &VRo, double &VRr, double &Eo, double &Er) {

    iVor_PTr_g(VRo, PR, TR, CstO);
    iVor_PTr_g(VRr, PR, TR, CstR);
    VR = (1.0 - AF) * VRo + AF * VRr;
    Eo = (CstO.BT + 1.0 - (CstO.GM / sqr(VRo) + CstO.BT + 1.0) / exp(CstO.GM / sqr(VRo))) * CstO.C4 / 2.0 / sqr(TR) / TR / CstO.GM; // +IV+ }
    Er = (CstR.BT + 1.0 - (CstR.GM / sqr(VRr) + CstR.BT + 1.0) / exp(CstR.GM / sqr(VRr))) * CstR.C4 / 2.0 / sqr(TR) / TR / CstR.GM;
}


void LK4::Prop_Ps_bg(double P, double s, double &v, double &T, double &h, double &x){
    LK_Prop_Ps_bg(P/1e3, s/1e3, v,T,h,x);
    h= h*1e3;
}

void LK3::LK_Prop_Ps_bg(double P, double s, double &v, double &T, double &h, double &x) {
    double vGaz, vLiq;
    double hGaz, hLiq;
    double sGaz, sLiq;
    double DPTs, DPVs;

    Ts_P(P,T,DPTs);
    ss_T(sGaz,T);
    if(sGaz < s){
        iLK_Prop_Ps_g(P,s,v,T,h);
    }
    else{
        LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
        x= (s- sLiq) / (sGaz - sLiq);
        v = (1-x)*vLiq + x *vGaz;
        h = (1-x) * hLiq + x *hGaz;
    }
}

void LK2::ss_T(double &ss,double T) {
    double Tr, P, Pr, VR, VRo, VRr, Eo, Er;
    Ps_T(P,T);
    Pr = P / kPc;
    Tr = T/ Tc;
    VRro_PTr_g(Pr,Tr,VR,VRo,VRr,Eo,Er);
    S_VT(VRo, VRr, Tr, Eo, Er, ss);
    //TODO check the pr, tr where does it comes from

}

void LK2::iLK_Prop_Ps_g(double P, double s, double &v, double &T, double &h){
    double si, Ti, DPT, DPV, CV, CP, SV;
    int i=0;
    Ts_P(P,T,DPT);
    LK_Prop_PT_g(P,T,v,h,si, DPV, DPT,CV,CP,SV);
    do{
        Ti=T;
        T=Ti*(1-(si-s)/CP);
        LK_Prop_PT_g(P, T, v, h, si, DPV, DPT, CV, CP, SV);
        if(i> 50){cout<<"error iLK_Prop_Ps_g "<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
    }until((abs(si-s) < 1e-7)&& (abs(Ti-T)<1e-3));

}

void LK2::Ps_T(double &P,double T){
    double Tr, LnPsatSimple, LnPsatRef;
    //TODO check wp
    if(T <= Tc){
        Tr= T /Tc;
        LnPsatRef = sqr((sqr(Tr)*Tr));
        LnPsatSimple=0.169347 * LnPsatRef - 1.28862 * ln(Tr) - 6.09648 / Tr + 5.92714;
        LnPsatRef =0.43577 * LnPsatRef - 13.4721 * ln(Tr) - 15.6875 / Tr + 15.2518;
        P= exp(w * LnPsatRef + LnPsatSimple) *kPc; // unite de kPc, kN/m2

    }
    else{
        P=(4.83002 * w + 5.823942) * (T - Tc) * kPc / Tc + kPc;
    }
}






double LK2::PRs_T(double &Tr, LeeKesslerCst &CSTS) {

    double prs=1.0;
    if (Tr < 0.9999999) {
        prs = exp(CSTS.P1 * sqr(sqr(Tr) * Tr) + CSTS.P2 * ln(Tr) + CSTS.P3 / Tr + CSTS.P4);
    }
    return prs;
}

double LK2::Zs(double &Tr, LeeKesslerCst &CSTS) {

    double zs;

    if (Tr < 0.9999999){
        zs = (1 - CSTS.ZCP) * pow(1 - pow((Tr - 0.3) / 0.7, CSTS.ZA), CSTS.ZB) + CSTS.ZCP;
    }
    else zs = CSTS.ZCP;
    return zs;
}

double LK2::root(double Pr, double Tr, double VRInf, double VRSup, double dVR, LeeKesslerCst &CSTS) {
    double dx, fr, fl, swap, xl, rts;
    int i=0;
    do{
        Por_VTr(fl, VRInf, Tr, CSTS);
        fl = fl - Pr;
        if (fl > 0.2){
            VRInf = VRInf *1.05;
        }
    }until(fl < 0.2);

    do{
        Por_VTr(fr, VRSup,Tr, CSTS);
        fr = fr - Pr;
        if (fr > 0.0){
            VRSup = VRSup *1.05;
        }
    } until (fr < 0.0);

    if (abs(fl)< abs(fr)){
        rts = VRInf;
        xl = VRSup;
        swap = fl;
        fl = fr;
        fr = swap;

    }
    else{
        xl = VRInf;
        rts= VRSup;
    }

    do{
        dx = (xl - rts)*fr / (fr -fl);
        xl = rts;
        fl = fr;
        rts= rts + dx;
        Por_VTr(fr, rts, Tr, CSTS);
        if (i > 100){
            cout<<"error root"<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
        fr = fr -Pr;
    }until((abs(dx)< dVR || fr==0));
    return rts;

}

void LK2::iVor_PTr_g2(double& VR, double PR, double TR, LeeKesslerCst &CSTS){
    int i = 0;
    double V1, Q, iP, DPV;

    VR = TR / PR;
    if (TR > 0.9) {
        if (TR < 1.199)
            VR = VR / Zc;
        if (PR < 0.999)
            VR = TR / PR;
    }
    V1 = VR;
    Q = 1 / VR;

    do {
        do {
            Por_VTr(iP, VR, TR, CSTS);
            if (i > 50) {
                cout << "error iVor_PTr_g2" << endl;
                exit(EXIT_FAILURE);
            }
            i++;
            DPVor_VTr(DPV, VR, TR, CSTS);
            VR = VR / 0.99;
        } until(DPV <= 0.0);
        Q = (iP - PR) / sqr(VR) / DPV + Q;
        V1 = VR;
        if (Q < 0)
            VR = V1 - 0.01;
        else
            VR = 1.0 / Q;
    }until( abs(1 - VR / V1) < 1.0e-6);
}

void LK1::Por_VTr(double &PR, double& VR,double &TR, LeeKesslerCst &CSTS) {
    double b,c,d,ep;

    d= sqr(VR);
    if (88.02*d <= CSTS.GM){
        ep= 0.0;
    }
    else {
        ep = exp(-CSTS.GM / sqr(VR));
    }

    b= -(((CSTS.B4/TR)+CSTS.B3)/TR + CSTS.B2) / TR +CSTS.B1;
    c= ((CSTS.C3 / TR)/TR -CSTS.C2)/TR + CSTS.C1;
    PR= (((CSTS.GM / d + CSTS.BT) * CSTS.C4 * ep / sqr(TR) / TR + (CSTS.D2 / TR + CSTS.D1) / d / VR + c) / VR
            + b + VR) * TR / d;

}
void LK1::DPVor_VTr(double &DPV,double VR, double TR, LeeKesslerCst &CSTS) {
    double b,c,d,cst;

    b= -(((CSTS.B4 /TR)+CSTS.B3)/TR + CSTS.B2)/TR + CSTS.B1;
    c= ((CSTS.C3/TR)/TR -CSTS.C2)/ TR + CSTS.C1;
    d = CSTS.D2 / TR + CSTS.D1;
    if(sqr(VR)*88.02 <= CSTS.GM){
        cst =0.0;
    }
    else{
        cst = exp(-CSTS.GM/ sqr(VR)) * CSTS.C4 / sqr(TR);
    }
    DPV= (((((cst * 2.0 / VR * sqr(CSTS.GM) - 6.0 * d * TR) / VR + cst * (2.0 * CSTS.BT - 5.0) * CSTS.GM) / sqr(VR)
              - (3.0 * c * TR + cst * 3.0 * CSTS.BT)) / VR - 2.0 * b * TR) / VR - TR) / sqr(VR);

}



void LK4::Prop_Ph_bg(double P, double h, double &v, double &T, double &s, double &x){
    P= P/1000;
    h = h/1000;
    iLK_Prop_Ph_bg(P,h,v,T,s,x);
    s= s*1000;

}

void LK3::iLK_Prop_Ph_bg(double P, double h, double &v, double &T, double &s, double &x){
    double vGaz, vLiq;
    double hGaz, hLiq;
    double sGaz, sLiq;
    double DPTs, DPVs;

    if (P <= kPc){
        LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
        if(h < hGaz){
            //biphase area
            x = (h- hLiq) / (hGaz - hLiq);
            v = (1 - x)*vLiq + x * vGaz;
            s = (1 - x)*sLiq + x *sGaz;
        }
        else{
            //overheating area
            iLK_Prop_Ph_g(P,h,v,T,s);
            x=1;
        }
    }
    else{
        //overcritical area
        iLK_Prop_Ph_g(P,h,v,T,s);
        x=1;
    }
}


void LK4::Prop_vs_bg(double v,double s,double &P,double &T,double &h,double &x){

    iLK_Prop_vs_bg(v,s/1e3, P,T,h,x);
    P=P*1000;
    h= h*1000;

}

void LK3::iLK_Prop_vs_bg(double v, double s, double &P, double &T, double &h, double &x) {
    double vi, Ti;
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq, DPTs, DPVs;


    int i=0;
    Ts_v(T,v);
    Ps_T(P,T);
    LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
    if (sGaz > s){
        do{


            Ti= T;
            x = (s-sLiq)/(sGaz-sLiq);
            T= Ti * pow(vi/v,0.05);
            Ps_T(P,T);
            LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
            if(i > 50){
                return ;
            }
            i++;
        } until ((abs(Ti - T) < 0.001) && (abs(vi - v) / v < 0.00001));
        x= (s-sLiq)/(sGaz-sLiq);
        h = hLiq + x*(hGaz-hLiq);

    }
    else{
        iiLK_Prop_vs_g(v,s,P,T,h);
    }



}

void LK2:: iiLK_Prop_vs_g(double v, double s, double &P, double &T, double &h){
    double Pi, si, m ,vi ;
    double DPV, DPT, cv, cp, sv;
    int i= 0;
    Ts_v(T,v);
    Ps_T(P,T);
    do{
        Pi=P;
        // find s = f(v,p)
        do{
            LK_Prop_PT_g(P, T, vi, h, si, DPV, DPT, cv, cp, sv);
            if (i> 100){
                return;
            }
            i++;
            T= (vi-v)*DPV /DPT + T;
        }until(abs(vi-v)< 1e-6);
        //end of s= f(v,p);
        m= cv/T / DPT;
        P= Pi- (si- s)/m;

    } until (abs(P-Pi)< 1.0 && abs(si-s)< 1e-5);
    LK_Prop_PT_g(P, T, vi, h, si, DPV, DPT, cv, cp, sv);

}

void LK2::iiLK_Prop_su_g (double s,double u, double &v,double & P,double & T,double & h){
    double Pinit, newP, newu, So, DPV, DPT, CV, CP, SV, m;
    int i;

    newP = kPc / 50;
    So = s;
    i = 0;
    do {
        Pinit = newP;
// ----- itŽration pour trouver u=f(P,s) -----}
        Ts_P(Pinit, T, DPT);
        LK_Prop_PT_g(Pinit, T, v, h, s, DPV, DPT, CV, CP, SV);
        do {
            T = T * (1 - (s - So) / CP);
            LK_Prop_PT_g(Pinit, T, v, h, s, DPV, DPT, CV, CP, SV);
            if(i>100) {cout<<"error iiLK_Prop_su_g "<<endl;
                exit(EXIT_FAILURE);
            }
            i++;
        } until (abs(s - So) < 1E-7);
// fin du calcul}
        newu = h - newP * v;
        m = -newP / DPV * CV / CP;
        newP = Pinit - (newu - u) / m;
    }until(abs(Pinit - newP) <= 1E-4);

            P = newP;
}


void LK2::Ts_v(double &Ts,double v) {
    double Ti, vi;
    double p, vL, hg,hl, sG, sL, DPT,DPV;
    int i=0;
    Ts= Tc*0.4;
    do{
        Ti= Ts;
        Ps_T(p,Ti);
        LK_Prop_P_s(p, vL, vi, Ti, hg, hl, sL, sG, DPT, DPV);
        if (i >= 100){cout<<"error Ts_v"<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
        Ts= Ti + (v-vi) *DPV /DPT;
    }until((abs(Ti - Ts) / Ts < 1E-3) && (abs(vi - v) / v < 1E-3));

}

void LK4::Prop_vh_bg(double v, double h, double &P, double &T, double &s, double &x){
    iLK_Prop_vh_bg(v,h/1e3,P,T,s,x);

    P= P*1E3;
    s= s*1e3;
}

void LK3::iLK_Prop_vh_bg(double v, double h, double &P, double &T, double &s, double &x){
    double vi, Ti;
    double vGaz, vLiq;
    double hGaz, hLiq;
    double sGaz, sLiq;
    double DPTs, DPVs;
    int i=0;
    Ts_v(T,v);
    Ps_T(P,T);
    LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
    if (hGaz > h){
        do{
            Ti= T;
            vi = (s-sLiq) /(sGaz -sLiq);
            T = Ti * pow(vi/v, 0.05);
            Ps_T(P,T);
            LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
            if (i>50){cout<<"error iLK_Prop_vh_bg"<<endl;
                exit(EXIT_FAILURE);
            }
            i++;
        }until((abs(Ti - T) < 0.001) && (abs(vi - v) / v < 0.0001));
        s= sLiq + x * (sGaz -sLiq);
        x = (h-hLiq) /(hGaz -hLiq);

    }
    else{
        iiLK_Prop_vh_g(v, h, P, T, s);
        x=1;
    }
}


void LK3::LK_Prop_vP_bg (double v,double P,double& T,double&  h,double& s,double&  x) {
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq, DPTs, DPVs;

    LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
    if (v > vGaz){
        // vapor zone}
        iLK_Prop_vP_g(v,P, T,h,s);
        x=1.0;}
    else{
    // biphasic zone}
    x = (v - vLiq) / (vGaz - vLiq);
    s = (1 - x) * sLiq + x * sGaz;
    h = (1 - x) * hLiq + x * hGaz;}
}
void LK3::LK_Prop_vT_bg (double v,double T,double&   P,double&  h, double& s,double&  x) {
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq, DPTs, DPVs;
        if (v / vc * Zc < 2) {
            Ps_T(P, T);
            LK_Prop_P_s(P, vLiq, vGaz, T, hLiq, hGaz, sLiq, sGaz, DPTs, DPVs);
            if (v > vGaz) {
                // vapor zone}
                iLK_Prop_vT_g(v, T, P, h, s);
                x = 1.0;
            }else{
                // biphasic zone}
                x = (v - vLiq) / (vGaz - vLiq);
                s = (1 - x) * sLiq + x * sGaz;
                h = (1 - x) * hLiq + x * hGaz;}
        }else{
            iLK_Prop_vT_g (v, T, P, h, s);
            x=1.0;}
}
void LK3::LK_TransProp_PT_l (double P,double T,double&  ViscDyn,double&  ViscCin,double&  Ro,double&  ThermCond,double&  Prandtl) {
    k = 0;
    double v, Cv, Cp, CpSat,
    Tr, Ror, Tstar, Vcu, MWu,
    ViscStar, Visc2Star, OmegaV, Fc,
    G1, G2, y,ro,h,s,DPv,DPT;

    LK_Prop_PT_l(P,T,v,ro,h,s,DPv,DPT);

    cout<<" v 644 = "<<v<<endl;
    cout<<" P 645 = "<<P<<endl;

    LK_Cvp_PT_l(P, T, Cp, Cv, CpSat);
    Cv = Cv * 1000;// ==> J/kg K}
    Cp = Cp * 1000;// ==> J/kg K}
    cout<<Cp<<endl;
//*** viscosities ***}
    Ro = 1 / v;
    Ror = vc * Ro;
    MWu = Mw * 1E-3;//[kg/mol]}
    Vcu = 1E3 * vc * Mw;//[cm3/mole]}
    Tr = T / Tc;
    Tstar = 1.2593 * Tr;
//eq. (9-4.3)}
    OmegaV = 1.16145 * XpwrY(Tstar, -0.14874) + 0.52487 * exp(-0.7732 * Tstar) + 2.16178 * exp(-2.43787 * Tstar);
//eq. (9-4.11)}
    Mur = 131.3 * DipM / sqrt(Vcu * Tc);
//eq. (9-4.10)}
    Fc = 1 - 0.2756 * w + 0.059035 * sqr(Mur) + k;
//eq. (9-6.18)}
    y = vc * Ro / 6;
//eq. (9-6.19)}
    G1 = (1 - y / 2) / XpwrI(1 - y, 3);
//eq. (9-6.20)}
    G2 = (E(1) / y * (1 - exp(-E(4) * y)) + E(2) * G1 * exp(E(5) * y) + E(3) * G1) / (E(1) * E(4) + E(2) + E(3));
//eq. (9-6.21)}
    Visc2Star = E(7) * sqr(y) * G2 * exp(E(8) + E(9) / Tstar + E(10) / sqr(Tstar));
//eq. (9-6.17)}
    ViscStar = sqrt(Tstar) / OmegaV * (Fc * (1 / G2 + E(6) * y)) + Visc2Star; //[µP]}
//eq. (9-6.16) Chung et al. method}
    ViscDyn = 1E-7 * ViscStar * 36.344 * sqrt(Mw * Tc) / XpwrY(Vcu, 2 / 3.0);
    ViscCin = ViscDyn / Ro;
//--- Thermal Conduction ---}
//eq. (10-9.5)}
    ThermCond = (1.11 / sqrt(Mw)) * (3 + 20 * XpwrY(1 - Tr, 2 / 3.0)) / (3 + 20 * XpwrY(1 - Tnb / Tc, 2 / 3.0));
    Prandtl = Cp * ViscDyn / ThermCond;
}


void LK3::LK_TransProp_vT_l (double v,double T,double&  ViscDyn,double&  ViscCin,double&  Ro,double&  ThermCond,double&  Prandtl) {
    k = 0;
    double Cv, Cp, CpSat,
            Tr, Ror, Tstar, Vcu, MWu,
            ViscStar, Visc2Star, OmegaV, Fc,
            G1, G2, y;
    double P= P_vT_l(v,T)/1000;

    cout<<T<<endl;
    LK_Cvp_PT_l(P, T, Cp, Cv, CpSat);
    Cv = Cv * 1000;// ==> J/kg K}
    Cp = Cp * 1000;// ==> J/kg K}
    cout<<" v 693 = "<<v<<endl;
    cout<<" P 694 = "<<P<<endl;
    cout<<Cp<<endl;
//*** viscosities ***}
    Ro = 1 / v;
    Ror = vc * Ro;
    MWu = Mw * 1E-3;//[kg/mol]}
    Vcu = 1E3 * vc * Mw;//[cm3/mole]}
    Tr = T / Tc;
    Tstar = 1.2593 * Tr;
//eq. (9-4.3)}
    OmegaV = 1.16145 * XpwrY(Tstar, -0.14874) + 0.52487 * exp(-0.7732 * Tstar) + 2.16178 * exp(-2.43787 * Tstar);
//eq. (9-4.11)}
    Mur = 131.3 * DipM / sqrt(Vcu * Tc);
//eq. (9-4.10)}
    Fc = 1 - 0.2756 * w + 0.059035 * sqr(Mur) + k;
//eq. (9-6.18)}
    y = vc * Ro / 6;
//eq. (9-6.19)}
    G1 = (1 - y / 2) / XpwrI(1 - y, 3);
//eq. (9-6.20)}
    G2 = (E(1) / y * (1 - exp(-E(4) * y)) + E(2) * G1 * exp(E(5) * y) + E(3) * G1) / (E(1) * E(4) + E(2) + E(3));
//eq. (9-6.21)}
    Visc2Star = E(7) * sqr(y) * G2 * exp(E(8) + E(9) / Tstar + E(10) / sqr(Tstar));
//eq. (9-6.17)}
    ViscStar = sqrt(Tstar) / OmegaV * (Fc * (1 / G2 + E(6) * y)) + Visc2Star; //[µP]}
//eq. (9-6.16) Chung et al. method}
    ViscDyn = 1E-7 * ViscStar * 36.344 * sqrt(Mw * Tc) / XpwrY(Vcu, 2 / 3);
    ViscCin = ViscDyn / Ro;
//--- Thermal Conduction ---}
//eq. (10-9.5)}
    ThermCond = (1.11 / sqrt(Mw)) * (3 + 20 * XpwrY(1 - Tr, 2 / 3)) / (3 + 20 * XpwrY(1 - Tnb / Tc, 2 / 3));
    Prandtl = Cp * ViscDyn / ThermCond;

}


void LK3::ROr_PTr_YW (double&  ROR,double PR,double TR,double ZCP) {
    if (abs(ZCP - 0.2801) < 0.2)
        ESTD(ROR, PR, TR, LKconst.YWCst9, ZCP);
    else if (abs(ZCP - 0.25) < 0.1)
        ESTD(ROR, PR, TR, LKconst.YWCst5, ZCP);
    else
        ESTD(ROR, PR, TR, LKconst.YWCst3, ZCP);
}


void LK3::VRro_PTr_l (double PR, double TR, double&  VR,double&  VRo,double&  VRr) {
    iVor_PTr_l(VRo, PR, TR, CstO);
    iVor_PTr_l(VRr, PR, TR, CstR);
    VR = (1.0 - AF) * VRo + AF * VRr;
}


void LK3::iVor_PTr_l (double & VR, double PR,double TR, LeeKesslerCst &CSTS){
    double V1, DPV, iP;
    int i=0;
    ROr_PTr_YW(VR, PR, TR, CSTS.ZCP); //
//    ROr_PTr_YW(VR, PR, TR, Zc); // TODO ici c'est Zcp et non zc
    VR = 1.0 / VR * Zc;
    do {
        Por_VTr(iP, VR, TR, CSTS);
        DPVor_VTr(DPV, VR, TR, CSTS);
        if(i>100){                  // cornelia.blanke: changed i>50 to i>100
            cout<<"error iVor_PTr_l"<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
        V1 = (PR - iP) / DPV;
        VR = VR + V1;
    }until(abs(V1 / VR) < 1.0e-5);
}

void LK3::ESTD (double & ROR, double PR, double TR, LKYWCstRec YWCst, double ZCP){

//	PROVIDES ESTIMATE OF DENSITY IN COMPRESSED LIQUID REGION USING YEN}
//	AND WOODS EQUATION A.I.CH.E.,JL.,VOL.12,  1966 FOR GIVEN ZcP	}
//}
//	CALLS SUBROUTINE Ps_T														}
//}
//	OUTPUT:	ROR	REDUCED DENSITY												}
//	INPUT:	PR		REDUCED PRESSURE												}
//				TR		REDUCED TEMPERATURE											}
//				DZcON	RECORD SUPPLYING YEN AND WOODS CONSTANTS				}
//				ZcP	CRITICAL COMPRESSIBILITY FACTOR, (Pc*vc/R*Tc)	}
    double i1, i2, i3, i4, // replace respectively A&E&I , B&F&J , G&K , D&H&L }
    PS, DELPR, DTR, DROR27, DELZc;

//P = PR * Pc;}
//T = TR * Tc;}
    Ps_T(PS, TR * Tc);
    DELPR = (PR - (PS / kPc));
    DTR = exp(ln(1.0 - TR) / 3.0); // this is third root of former DTR -IV- }
    if (DELPR >= 0.1e-05) {

        i1 = (((-2.198 * DTR + 3.699) * DTR - 0.646) * DTR - 1.626) * DTR + 0.714;
        i2 = exp(ln(TR) * 2.0967) / (1.0 + 0.8 * exp(ln(-ln(TR)) * 0.441)) * 0.268;
        i3 = exp(ln(1.01 - TR) * 0.75) * exp(-7.848 * (1.01 - TR)) * 4.221 + 0.05;
        i4 = (((-47.38 * DTR + 114.44) * DTR - 103.79) * DTR + 45.22) * DTR - 10.6;
        if (DELPR < 0.2)
            DROR27 = (i1 + i2 * ln(0.2) + i3 * exp(i4 * 0.2)) * DELPR / 0.2;
        else
            DROR27 = i1 + i2 * ln(DELPR) + i3 * exp(i4 * DELPR);
        i1 = (((YWCst.A5 * DTR + YWCst.A4) * DTR + YWCst.A3) * DTR + YWCst.A2) * DTR + YWCst.A1;
        i2 = (((YWCst.B5 * DTR + YWCst.B4) * DTR +YWCst. B3) * DTR + YWCst.B2) * DTR + YWCst.B1;
        i3 = ((YWCst.C4 * TR + YWCst.C3) * TR + YWCst.C2) * TR + YWCst.C1;
        i4 = (((YWCst.D5 * DTR + YWCst.D4) * DTR + YWCst.D3) * DTR + YWCst.D2) * DTR + YWCst.D1;
        if (DELPR < 0.2)
            DELZc = (exp(0.2 * i4) * i3 + ln(0.2) * i2 + i1) * DELPR / 0.2;
        else
            DELZc = exp(i4 * DELPR) * i3 + ln(DELPR) * i2 + i1;
        if (abs(ZCP - 0.27) < 0.1)
            DELZc = 0.0;
    }; // IF DELPRÉ }
    i1 = ((-1522.06 * ZCP + 989.625) * ZCP - 214.578) * ZCP + 17.4425;
    if (ZCP <= 0.26)
        i2 = ((641.0 * ZCP + 501.0) * ZCP- 402.063) * ZCP + 60.2091;
    else
        i2 = ((-384.211 * ZCP + 107.4844) * ZCP + 13.6377) * ZCP - 3.28257;
    DTR = (((0.93 - i2) * sqr(DTR) + i2) * DTR + i1) * DTR + 1.0;
    if (DELPR >= 0.1e-05)
        ROR = DTR + DROR27 + DELZc;
    else
        ROR = DTR;
}

void LK2::iiLK_Prop_vh_g(double v, double h, double &P, double &T, double &s){
    double vi, Pi, Ti,hi, m;
    double DPV, DPT, CV, CP, SV;
    int i=0;
    Ts_v(T,v);
    Ps_T(P,T);
    do{
        Pi=P;
        do{
            Ti = T;
            LK_Prop_PT_g(P, Ti, vi, hi, s, DPV, DPT, CV, CP, SV);
            if(i > 100){
                return;
            }
            T= Ti + (vi -v)*DPV/DPT;
        } until ((abs(Ti - T) < 1E-4) && (abs(vi - v) < 1E-8));
        m= CV/DPT + v;
        P= Pi- (hi -h)/m;
    } until ((abs(P - Pi) < 1E-4) & (abs(hi - h) < 1E-4));
    LK_Prop_PT_g(P, T, vi, hi, s, DPV, DPT, CV, CP, SV);

}


void LK1::DPT_VTor (double VRo,double VRr, double TR, double &DPTo,double &DPTr,double &DPT){
        DPTor_VTr(DPTo, VRo, TR,CstO);
        DPTor_VTr(DPTr, VRr, TR,CstR);
        DPT = (1.0 - AF) * DPTo + AF *DPTr;
    DPT =DPT *kPc/Tc;

}

void LK1::DPTor_VTr (double & DPT, double VR,double TR, LeeKesslerCst &CSTS){
    double EP, CTS, BOT, COT;
    EP = exp(-CSTS.GM / sqr(VR));
    CTS = (EP * 2.0 * CSTS.C4 / TR) / sqr(TR);
    BOT = ((2.0 * CSTS.B4 / TR + CSTS.B3) / sqr(TR)) + CSTS.B1;
    COT = CSTS.C1 - 2.0 * CSTS.C3 / sqr(TR) / TR;
    DPT = ((((CSTS.D1 / 5.0 / VR - (2.0 * CSTS.C4 / sqr(TR) / TR) * CSTS.GM * EP) / sqr(VR) + (COT / 0.2e+01 - CTS * CSTS.BT)) / VR + BOT) / VR + 0.1e+01) / VR;

}

void LK1::DPV_VTor (double VRo, double VRr, double TR, double &DPVo, double &DPVr,double  &DPV){
    DPVor_VTr(DPVo, VRo, TR, CstO);
    DPVor_VTr(DPVr, VRr, TR, CstR);
    DPV = (1.0 - AF) * DPVo + AF * DPVr;
    DPV = DPV * Zc * kPc / vc;
}


void LK2::iiLK_Prop_hs_g (double h,double s, double &v, double &P, double &T){
    double Pi, m, Ti, si, hi,DPV, DPT, CV, CP, SV;
    int i=0;
    P = kPc / 10;
    iLK_Prop_Ph_g(P, h, v, T, si);
    do{
        Pi = P;
        //----- iteration pour trouver s=f(P,h) -----
        do{
            Ti = T;
            LK_Prop_PT_g(Pi, T, v, hi, si, DPV, DPT, CV, CP, SV);
            if (i>=100){cout<<"error iiLK_Prop_hs_g"<<endl;
                exit(EXIT_FAILURE);
            }
            i++;
            T = Ti - (hi - h) / CP;
        }until((abs(hi - h) < 1E-5) && (abs(Ti - T) < 1E-3));
        //------ fin de l'iteration intermediaire -------
        m = -v / T;
        P = Pi - (si - s) / m;
    }until((abs(Pi - P) < 1E-3) && (abs(si - s) < 1E-3));
    LK_Prop_PT_g(P, T, v, h, s, DPV, DPT, CV, CP, SV);
}


void LK2::H_VT(double VRo,double  VRr,double  PR,double  TR,double  Eo,double  Er,double  &H){
    double HO, DHO, DHR, DSH, DH, T;
    T = TR * Tc;
    HO = T*(ALPHA+T*(BETA/2.0+T*(GAMMA/3.0+DELTA/4.0*T)))-Tref*(ALPHA+Tref*(BETA/2.0+Tref*(GAMMA/3.0+DELTA/4.0*Tref)))+ HREFO;
    DHO = -(CstO.B2 + (2.0 * CstO.B3 + 3.0 * CstO.B4 / TR) / TR) / VRo - (CstO.C2 - 3.0 * CstO.C3 / sqr(TR)) / 2.0 / sqr(VRo) + CstO.D2 / 5.0 / sqr(sqr(VRo)) / VRo + (PR * VRo/TR -1.0 + 3.0 * Eo ) * TR;
    DHR = -(CstR.B2 + (2.0 * CstR.B3 + 3.0 * CstR.B4 / TR) / TR) / VRr - (CstR.C2 - 3.0 * CstR.C3 / sqr(TR)) / 2.0 / sqr(VRr) + CstR.D2 / 5.0 / sqr(sqr(VRr)) / VRr + (PR * VRr/TR -1.0 + 3.0 * Er ) * TR;
    DH = (1.0 - AF) * DHO + AF * DHR;
    DSH = DH * LKconst.RBAR * Tc / Mw;
    H = HO + DSH;
}

void LK2::S_VT(double VRo,double  VRr,double TR,double Eo,double Er, double  &S){
    double SO, POR, DSO, DSR, DS, DSS, T;
    T = TR * Tc;
    SO = (sqr(T) * T - sqr(Tref) * Tref) * DELTA / 3.0 + (sqr(T) - sqr(Tref)) * GAMMA / 2.0 + ln(T / Tref) * ALPHA + (T - Tref) * BETA + SREFO;
    POR =LKconst.PO / kPc;
    DSO = (-(CstO.B1 + (CstO.B3 + 2.0 * CstO.B4 / TR) / sqr(TR)) - (CstO.C1 - 2.0 * CstO.C3 / sqr(TR) / TR) / VRo / 2.0 - CstO.D1 / 5.0 / sqr(sqr(VRo))) / VRo + 2.0 * Eo + ln(POR * VRo / TR);
    DSR = (-(CstR.B1 + (CstR.B3 + 2.0 * CstR.B4 / TR) / sqr(TR)) - (CstR.C1 - 2.0 * CstR.C3 / sqr(TR) / TR) / VRr / 2.0 - CstR.D1 / 5.0 / sqr(sqr(VRr))) / VRr + 2.0 * Er + ln(POR * VRr / TR);
    DS = (1.0 - AF) * DSO + AF * DSR;
    DSS = DS * LKconst.RBAR / Mw;
    S = SO + DSS;
}

void  LK2::CVP_VT (double VRo,double  VRr,double  Eo,double  Er,double  DPVO,double  DPVR,double  DPTO,double  DPTR, double  TR, double  &CV, double &CP){
    double T, CPO, CVO,DCPO, DCPR, DCP,DCVO, DCVR, DCV;
    T = TR * Tc;
    CPO = ((DELTA * T + GAMMA) * T + BETA) * T + ALPHA;
    CVO = CPO - LKconst.RBAR / Mw;
    // calcul de cv
    DCVO = 3.0 * CstO.C3 / sqr(VRo) / sqr(TR) / TR - (3.0 * CstO.B4 / TR + CstO.B3) * 2.0 / VRo / sqr(TR) + 6.0 * Eo;
    DCVR = 3.0 * CstR.C3 / sqr(VRr) / sqr(TR) / TR - (3.0 * CstR.B4 / TR + CstR.B3) * 2.0 / VRr / sqr(TR) + 6.0 * Er;
    DCV = (1.0 - AF) * DCVO + AF * DCVR;
    CV = CVO - DCV * LKconst.RBAR / Mw;
    // calcul de cp
    DCPO = sqr(DPTO) / DPVO * TR + DCVO + 1.0;
    DCPR = sqr(DPTR) / DPVR * TR + DCVR + 1.0;
    DCP = (1 - AF) * DCPO + AF * DCPR;
    CP = CPO - DCP * LKconst.RBAR / Mw;
}

void LK2::Ts_P(double P, double &Ts, double &DPTS){
    double Tinit, newP;
    int i=0;
    Ts = Tc / (-3.0 / 7.0 * ln(P / kPc) / ln(10) / (1.0 + w) + 1.0);
    do{
        Tinit = Ts;
        Ps_T(newP, Tinit);
        DPTS = DPTS_PT_sg(newP, Tinit);
        if (i>=100){cout<<"error Ts_P"<<endl;
            cout<<Ts<<endl;
            cout<<P<<endl;
            //exit(EXIT_FAILURE);       //TODO: check problem with convergence!
            return;
        }
        i++;
        Ts = Tinit - (newP - P) / DPTS;
    }until( (abs(Ts - Tinit)< 0.001) && ((abs(P - newP) / P )< 0.0001));

}

double LK2::DPTS_PT_sg (double Ps, double Ts){
    double Tr6,DPTS;
    if (Ts <= Tc){
        Tr6 = sqr(sqr(Ts / Tc)) / sqr(Tc) * Ts;
        DPTS= (((15.6875 * Tc / Ts - 13.4721) / Ts + 2.61462 * Tr6) * w + 6.09648 * Tc / sqr(Ts) - 1.28862 / Ts + 1.016082 * Tr6) * Ps;
    }
    else
        DPTS= (4.83002 * w + 5.823942) * kPc / Tc;
    return DPTS;
}

void LK2::LK_Prop_P_s (double P, double &vLiq, double &vGaz, double &Ts, double &hLiq, double &hGaz, double &sLiq, double &sGaz, double &DPTs, double &DPVs){
    double CVG, CPG, DPT, SV, dhVap, dsVap;
    if(P <= kPc){
        Ts_P(P, Ts, DPTs);
        LK_Prop_PT_g(P, Ts, vGaz, hGaz, sGaz, DPVs, DPT, CVG, CPG, SV);
        vLiq = 1.0 / (ADA1 * Ts + ADA0 - 1.0 / vGaz);
        dhVap = (vGaz - vLiq) * Ts * DPTs;
        hLiq = (hGaz - dhVap);
        dsVap = dhVap / Ts;
        sLiq = (sGaz - dsVap);
    }
    else
        cout<<"LK_Prop_P_s inadapte : P>Pc !!!"<<endl;
}


void LK2::iVor_PTr_g(double &VR, double PR,double TR, LeeKesslerCst &CSTS ){
    //méthode de Newton
    //estimation de VR

    double VRi, VRs, PRi, DPV;
    int i=0;
    PRi = PRs_T(TR, CSTS);
    VRs = Zs(TR, CSTS) * TR / PRi;
    VR = VRs;
    do{
        VRi = VR;
        Por_VTr(PRi, VR, TR, CSTS);
        DPVor_VTr(DPV, VR, TR, CSTS);
        if (i >=100) {cout<<"error iVor_PTr_g"<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
        if (DPV > -1E-50){
            if (PR > 1)
                VR = root(PR, TR, 0.5 * CSTS.ZCP, 0.75 * CSTS.ZCP, 1E-5, CSTS);
            else
                VR = root(PR, TR, VRs, TR / PR, 1E-5, CSTS);
            break;
        }
        else
            VR=VRi - (PRi - PR)/DPV;
    } until ((abs(VR - VRi) / VR < 1e-10) && (abs(PR - PRi) / PR < 1e-10));
}


void LK3::InitLKFluid(){
    double PREF, VREF;
    double DPTRF, VREFO, VREFR, RHGRF, ER, EO, DHRFO, DHRFR, DH, DSRFO, DSRFR, DS, RHFRF, VFREF, TCM, PCM, VCM, HCM, SCM, DPVM, DPTM, CV, CP, SV, VFC, HFC, SFC, HFG, SFG;
    double POR, Tnbr,LnTnbr, PL;

    // TODO bizarre ici les deux PREF


    PREF = Tnb;
    PREF = Tc;
    Tnbr = Tnb / Tc;
    LnTnbr=ln(Tnbr);
    ER = sqr(sqr(Tnbr) * Tnbr);

    double wp = (ln(LKconst.PO / kPc) + 6.09648 / Tnbr + LnTnbr * 1.28862 - 0.169347 * ER - 5.92714) /(-15.6875 / Tnbr - 13.4721 *LnTnbr + 0.43577 * ER + 15.2518);

    Tref = LKconst.TDTM;
    Ps_T(PREF, Tref);
    Ts_P(PREF, Tref, DPTRF);

    AF = wp / LKconst.WR;
    Zc = kPc * vc * Mw / LKconst.RBAR / Tc;
    PREF = PREF / kPc;
    Tref = Tref / Tc;
    POR = LKconst.PO / kPc;

    VRro_PTr_g(PREF, Tref, VREF, VREFO, VREFR, EO, ER);



    VREF = VREF * vc / Zc;
    RHGRF = 1.0 / VREF;
    EO = CstO.GM / sqr(VREFO); // i
    EO = (1.0 - (CstO.BT + 1.0 + EO) / exp(EO) + CstO.BT) / sqr(Tref) / Tref * CstO.C4 / CstO.GM / 2.0;
    ER = CstR.GM / sqr(VREFR);
    ER = (1.0 - (CstR.BT + 1.0 + ER) / exp(ER) + CstR.BT) / sqr(Tref) / Tref * CstR.C4 / CstR.GM / 2.0; // -IV-
    DHRFO = (((3.0 * CstO.B4 / Tref + 2.0 * CstO.B3) / Tref + CstO.B2) / VREFO - (3.0 * CstO.C3 / sqr(Tref) - CstO.C2) / 2.0 / sqr(VREFO) - CstO.D2 / 5.0 / sqr(sqr(VREFO)) / VREFO - 3.0 * Tref * EO - (PREF * VREFO) + Tref); // -IV-
    DHRFR = (((3.0 * CstR.B4 / Tref + 2.0 * CstR.B3) / Tref + CstR.B2) / VREFR - (3.0 * CstR.C3 / sqr(Tref) - CstR.C2) / 2.0 / sqr(VREFR) - CstR.D2 / 5.0 / sqr(sqr(VREFR)) / VREFR - 3.0 * Tref * ER - (PREF * VREFR) + Tref);
    DH = (1.0 - AF) * DHRFO + AF * DHRFR;
    DSRFO = ((2.0 * CstO.B4 / Tref + CstO.B3) / sqr(Tref) + CstO.B1) / VREFO - (CstO.C3 / sqr(Tref) / Tref - 0.5 * CstO.C1) / sqr(VREFO) + 0.2 * CstO.D1 / sqr(sqr(VREFO)) / VREFO - 2.0 * EO - ln(POR * VREFO / Tref); //   -IV-
    DSRFR = ((2.0 * CstR.B4 / Tref + CstR.B3) / sqr(Tref) + CstR.B1) / VREFR - (CstR.C3 / sqr(Tref) / Tref - 0.5 * CstR.C1) / sqr(VREFR) + 0.2 * CstR.D1 / sqr(sqr(VREFR)) / VREFR - 2.0 * ER - ln(POR * VREFR / Tref);
    DS = (1.0 - AF) * DSRFO + AF * DSRFR;

    PREF = PREF * kPc;
    Tref = Tref * Tc;

    Ps_T(PL, TL);
    ADA1 = (2 / vc - (roL + PL / TL / LKconst.RBAR * Mw)) / (Tc - TL);
    ADA0 = 2 / vc - ADA1 * Tc;

    RHFRF = ADA1 * Tref + ADA0 - RHGRF;
    VFREF = 1.0 / RHFRF;
    href = Tref * DPTRF * (VREF - VFREF) + 200; // 200 kJ/kg
    sref= ((href - 200) / Tref) + 1;//1 kJ/kg K
    HREFO = href + DH * LKconst.RBAR * Tc / Mw;
    SREFO = sref + DS * LKconst.RBAR / Mw;
    TCM = Tc - 2;

            Ps_T(PCM, TCM);
            LK_Prop_PT_g(PCM, TCM, VCM, HCM, SCM, DPVM, DPTM, CV, CP, SV);
            LK_Prop_P_sl(PCM, VFC, HFC, SFC, HFG, SFG, VCM, HCM, SCM);
            hc = (HCM + HFC) / 2.0;
            sc = (SCM + SFC) / 2.0;
    }


void LK3::LK_Prop_PT_l (double P,double T,double& V, double& RO,double&  H,double&  S,double&  DPV,double&  DPT) {
    double PS, PRs, VG, HG, SG, CV, CP, SV, VF, HF, SF, HFG, SFG,VRos, EPO, EPFO, DHOL, DSOL,
    VRrs, EPR, EPFR, DHRL, DSRL,
    VRo, VRr,VRs, VR, PR, TR, DPVo, DPVr, DPTo, DPTr;

    Ps_T(PS, T);
    PRs = PS / kPc;
    PR = P / kPc;
    TR = T / Tc;
    if ((TR > 0.9999) )
        cout<<"LK3:: Out of domain:   LK_Prop_PT_l IS WRONG SUBROUTINE USE LK_Prop_PT_g or Twoph"<<endl;
    else{

//calcul de v}
        VRro_PTr_l(PR, TR, VR, VRo, VRr);
        V = VR * vc / Zc;

        LK_Prop_PT_g(PS, T, VG, HG, SG, DPV, DPT, CV, CP, SV);
    //TODO check
        LK_Prop_P_sl(PS, VF, HF, SF, HFG, SFG, VG, HG, SG);
        CV = 0.0;
        CP = 0.0;
        SV = 0.0;

    //calcul de v sat.}
        VRro_PTr_l(PRs, TR, VRs, VRos, VRrs);
        V = V * VF / (VRs * vc / Zc);
    //calcul de h}
        RO = 1.0 / V;
        EPO = exp(-CstO.GM / sqr(VRo));
        EPFO = exp(-CstO.GM / sqr(VRos));
        EPR = exp(-CstR.GM / sqr(VRr));
        EPFR = exp(-CstR.GM / sqr(VRrs));
        DSOL = 1 / sqr(sqr(VRo)) / VRo - 1 / sqr(sqr(VRos)) / VRos; // i }
        DSRL = 1 / sqr(sqr(VRr)) / VRr - 1 / sqr(sqr(VRrs)) / VRrs; // i }
        DHOL = ((CstO.BT + 1.0) * (EPO - EPFO) - CstO.GM * (EPO / sqr(VRo) - EPFO / sqr(VRos))) * (-1.5) * CstO.C4 / CstO.GM / sqr(TR) - (3.0 *CstO.B4 / sqr(TR) + 2.0 * CstO.B3 / TR + CstO.B2) * (1.0 / VRo - 1.0 / VRos) + (1.0 / sqr(VRo) - 1.0 / sqr(VRos)) * (3.0 * CstO.C3 / sqr(TR) - CstO.C2) / 2.0 + (DSOL) * CstO.D2 / 5.0 + TR * (PR * VRo / TR - PRs * VRos / TR);
        DHRL = ((CstR.BT + 1.0) * (EPR - EPFR) - CstR.GM * (EPR / sqr(VRr) - EPFR / sqr(VRrs))) * (-1.5) * CstR.C4 / CstR.GM / sqr(TR) - (3.0 * CstR.B4 / sqr(TR) + 2.0 * CstR.B3 / TR + CstR.B2) * (1.0 / VRr - 1.0 / VRrs) + (1.0 / sqr(VRr) - 1.0 / sqr(VRrs)) * (3.0 * CstR.C3 / sqr(TR) - CstR.C2) / 2.0 + (DSRL) * CstR.D2 / 5.0 + TR * (PR * VRr / TR - PRs * VRrs / TR);// -IV- : CHECK THEM ?IV? }
        H = HF + ((1.0 - AF) * DHOL + AF * DHRL) * LKconst.RBAR * Tc / Mw;
    //calcul de s}
        DSOL = (1.0 / VRos - 1.0 / VRo) * ((CstO.B4 * 2.0 / TR + CstO.B3) / sqr(TR) + CstO.B1) + (1.0 / sqr(VRo) - 1.0 / sqr(VRos)) * (2.0 / sqr(TR) / TR * CstO.C3 - CstO.C1) - (DSOL) * CstO.D1 - (CstO.BT + 1.0) * (EPO - EPFO) / sqr(TR) / TR * CstO.C4 / CstO.GM - (EPO / sqr(VRo) - EPFO / sqr(VRos)) / sqr(TR) / TR * CstO.C4 + ln(VRo / VRos);
        DSRL = (1.0 / VRrs - 1.0 / VRr) * ((CstR.B4 * 2.0 / TR + CstR.B3) / sqr(TR) + CstR.B1) + (1.0 / sqr(VRr) - 1.0 / sqr(VRrs)) * (2.0 / sqr(TR) / TR * CstR.C3 - CstR.C1) - (DSRL) * CstR.D1 - (CstR.BT + 1.0) * (EPR - EPFR) / sqr(TR) / TR * CstR.C4 / CstR.GM - (EPR / sqr(VRr) - EPFR / sqr(VRrs)) / sqr(TR) / TR * CstR.C4 + ln(VRr / VRrs); // -IV- : CHECK THEM ?IV? }
        S = SF + ((1.0 - AF) * DSOL + AF * DSRL) * LKconst.RBAR / Mw;
    //calcul de DPV}
        DPV_VTor(VRo, VRr, TR, DPVo, DPVr, DPV);
    //calcul de DPT}
        DPT_VTor(VRo, VRr, TR, DPTo, DPTr, DPT);
    }
}
double LK3::LK_v_PT_l (double P,double T) {
    double PS, VRo, VRr, VR;
    cout<<"TLK3"<<T<<endl;
    Ps_T(PS,T);

    PS = PS / kPc;
    P = P / kPc;
    T = T / Tc;
    if ((P - PS) < 0){
        cout<<"Out of domain:   LK_v_PT_l IS WRONG SUBROUTINE USE LK_Prop_PT_g or Twoph"<<endl;
        return 0;}
    else if (T > 0.9999){
        cout<<"Out of domain:   LK_v_PT_l IS WRONG SUBROUTINE USE LK_Prop_PT_g or Twoph"<<endl;
        return 0;}
    else{
        VRro_PTr_l(P, T, VR, VRo, VRr);
        return VR * vc / Zc;
    }
}
void LK3::iLK_T_Ph_l (double P,double H,double& T) {
    int K;
    double HT, DPT, V, RO, S, DPV, CPL, CVL, CPF;
        K = 1;
        HT = H;
    if (P <= kPc * 0.9999)
                Ts_P(P, T, DPT);
    else{
    T = Tc * 0.97;
    T = T - 0.1;}

    do { //1}
        LK_Prop_PT_l(P, T, V, RO, H, S, DPV, DPT);
        if (abs(HT - H) < 0.1) break;
        LK_Cvp_PT_l(P, T, CPL, CVL, CPF);
        T = T + (HT - H) / CPL;
        K = K + 1;
    }until(K > 20);
//writeln('iLK_T_Ph_l IterFailed TO CONVER');}
}
void LK3::iLK_Th_Ps_l (double& T,double& H, double P,double S) {
    int I;
    double DPT, SN, V, RO, DPV, CPL, CVL, CPF;
        I = 1;
    if (P <= kPc * 0.9999)
                Ts_P(P, T, DPT);
    else{
    T = Tc * 0.97;
    SN = S;
    T = T - 0.1;}
    do {
        LK_Prop_PT_l(P, T, V, RO, H, S, DPV, DPT);
        if (abs(S - SN) < 0.2e-03){cout<<"error iLK_Th_Ps_l "<<endl;
            exit(EXIT_FAILURE);
        }
        T = T * (1.0 + (SN - S) / CPL);
        I = I + 1;
    }until (I > 20);
//writeln('iLK_Th_Ps_l FAILED TO CONVERGE');}
}
void LK3::LK_Cvp_PT_l (double P,double T, double& CPL, double& CVL, double& CPF) {
    double PS, TB, TU, PB, PU, V, H, S, DPE, DPF, CV, CP, SV, VFU, HFU, SF, HFG, SFG, VFB, HFB, //}
    CPSF, VF, DPG, RO, DPVS, DPTS, CVF, VRo, VRr, VR, VOS, VRS, EPO, EPFO, EPR, EPFR, //}
    DCVO, DCVR, DCV, DPVO, DPVR, DPV, DPTO, DPTR, DPT, DCPO, DCPR, DCP,
    TR, PR; // reduced to critical values ( TR = T/Tc ) }

    if (T > Tc * 0.999)
        cout<<" ATTEMPT HAS BEEN MADE TO CALCULATE LIQUID,SPEC HEAT IN THE GASEOUS REGION"<<endl;
    else{
        /*
        Ps_T(PS, T);
        if(((PS - P) / PS > 1.0E-6))
            cout<<" ATTEMPT HAS BEEN MADE TO CALCULATE,SATD LIQUID SPEC HEAT WHEN T>TC"<<endl;*/
        if(kPc<P){
            cout<<" ATTEMPT HAS BEEN MADE TO CALCULATE,SATD LIQUID SPEC HEAT WHEN T>TC"<<endl;
        }
        else{
            TB = T - 0.1;
            TU = T + 0.1;
            Ps_T(PU, TU);
            Ps_T(PB, TB);
            LK_Prop_PT_g(PU, TU, V, H, S, DPE, DPF, CV, CP, SV);
            LK_Prop_P_sl(PU, VFU, HFU, SF, HFG, SFG, V, H, S);
            LK_Prop_PT_g(PB, TB, V, H, S, DPE, DPF, CV, CP, SV);
            LK_Prop_P_sl(PB, VFB, HFB, SF, HFG, SFG, V, H, S);
            CPSF = (HFU - HFB) * 5.0;
            VF = (VFU + VFB) / 2.0;
            Ps_T(PS, T);
            Ts_P(PS, T, DPG);
            PS = PS * 1.0000011; // (1.0 + 1.1E-6) }
            LK_Prop_PT_l(PS, T, V, RO, H, S, DPVS, DPTS);
            CPF = CPSF - (T * DPTS / DPVS + VF) * DPG;
            CVF = sqr(DPTS) * T / DPVS + CPF;
            TR = T / Tc;
            PR = P / kPc;
            PS = PS / kPc;

            VRro_PTr_l(PR, TR, VR, VRo, VRr);
            VRro_PTr_l(PS, TR, VR, VOS, VRS);
        //calcul de cv}
            EPO = exp(-CstO.GM / sqr(VRo));
            EPFO = exp(-CstO.GM / sqr(VOS));
            EPR = exp(-CstR.GM / sqr(VRr));
            EPFR = exp(-CstR.GM / sqr(VRS));
            DCVO = ((((1.0 + CstO.BT) * (EPFO - EPO) / CstO.GM + EPFO / sqr(VOS) - EPO / sqr(VRo)) * CstO.C4 + (1.0 / sqr(VRo) - 1.0 / sqr(VOS)) * CstO.C3) / TR * 3.0 - (6.0 * CstO.B4 / TR + 2.0 * CstO.B3) * (1.0 / VRo - 1.0 / VOS)) / sqr(TR); // +IV+ }
            DCVR = ((((1.0 + CstR.BT) * (EPFR - EPR) / CstR.GM + EPFR / sqr(VRS) - EPR / sqr(VRr)) * CstR.C4 + (1.0 / sqr(VRr) - 1.0 / sqr(VRS)) * CstR.C3) / TR * 3.0 - (6.0 * CstR.B4 / TR + 2.0 * CstR.B3) * (1.0 / VRr - 1.0 / VRS)) / sqr(TR);
            DCV = (1.0 - AF) * DCVO + AF * DCVR;
            CVL = CVF - DCV * LKconst.RBAR / Mw;
        //calcul de cp}
            DPT_VTor(VRo, VRr, TR, DPTO, DPTR, DPT);
            DPV_VTor(VRo, VRr, TR, DPVO, DPVR, DPV);
            DCPO = sqr(DPTO) / DPVO * TR + DCVO + 1;
            DCPR = sqr(DPTR) / DPVR * TR + DCVR + 1;
            DCP = (1.0 - AF) * DCPO + AF * DCPR;
            CPL = CPF - DCP * LKconst.RBAR / Mw;
}}}



void LK3::LK_Prop_P_sl(double P, double &vf, double& hf, double& sf, double& hfg, double& sfg, double vg, double hg,
                       double sg) {
    double T, DPT;
            Ts_P(P, T, DPT);
        vf = 1.0 / (ADA1 * T + ADA0 - 1.0 / vg);
    sfg = (vg - vf) * DPT;
    hfg = sfg * T;
    sf = sg - sfg;
    hf = hg - hfg;
}


void LK1::Init_LKConst(){

    LKconst.TDTM = 273.15; //;{  T AT WHICH SAT. LIQUID ENTHALPY=200KJ/KG AND SAT. LIQUID ENTROPY=1 KJ/KG K}
    LKconst.RBAR = 8.314; //{[J/(mol K)]}
    LKconst.PO = 101.3; // {reference pressure for entropy = 1 atm. [kN/m2]	-RZ-}
    LKconst.WR = 0.3978;
    //{simple fluid}


    CstO.BT = 0.65392;
    CstO.GM = 0.60167e-01;
    CstO.B1 = 0.1181193;
    CstO.B2 = 0.265728;
    CstO.B3 = 0.15479;
    CstO.B4 = 0.30323e-01;
    CstO.C1 = 0.236744e-01;
    CstO.C2 = 0.186984e-01;
    CstO.C3 = 0.0;
    CstO.C4 = 0.42724e-01;
    CstO.D1 = 0.155488e-04;
    CstO.D2 = 0.623689e-04;
    CstO.ZCP = 0.2901e-00;
    CstO.ZA = 2.40005;
    CstO.ZB = 0.529594;
    CstO.P1 = 0.169347;
    CstO.P2 = -1.28862;
    CstO.P3 = -6.09648;
    CstO.P4 = 5.92714;
    //{reference fluid}
    CstR.BT = 1.226e-00;
    CstR.GM = 0.3754e-01;
    CstR.B1 = 0.2026579e-00;
    CstR.B2 = 0.331511e-00;
    CstR.B3 = 0.27655e-01;
    CstR.B4 = 0.203488e-00;
    CstR.C1 = 0.313385e-01;
    CstR.C2 = 0.503618e-01;
    CstR.C3 = 0.16901e-01;
    CstR.C4 = 0.41577e-01;
    CstR.D1 = 0.48736e-04;
    CstR.D2 = 0.740336e-05;
    CstR.ZCP = 0.2551e-00;
    CstR.ZA = 3.49445;
    CstR.ZB = 0.561744;
    CstR.P1 = 0.3426963;
    CstR.P2 = -6.64782138;
    CstR.P3 = -12.3369675;
    CstR.P4 = 11.994306;
//YWCst3
    LKconst.YWCst3.A1 = 0.89e-01;
    LKconst.YWCst3.A2 = -0.4344e-00;
    LKconst.YWCst3.A3 = 0.7915e-00;
    LKconst.YWCst3.A4 = -0.7654e-00;
    LKconst.YWCst3.A5 = 0.3367e-00;
    LKconst.YWCst3.B1 = 0.674e-01;
    LKconst.YWCst3.B2 = -0.6109e-01;
    LKconst.YWCst3.B3 = 0.6261e-01;
    LKconst.YWCst3.B4 = -0.2378e-00;
    LKconst.YWCst3.B5 = 0.1665e-00;
    LKconst.YWCst3.C1 = -0.1393e-01;
    LKconst.YWCst3.C2 = -0.3459e-02;
    LKconst.YWCst3.C3 = -0.1611e-00;
    LKconst.YWCst3.C4 = 0.0e-00;
    LKconst.YWCst3.D1 = -0.655e+01;
    LKconst.YWCst3.D2 = 0.78027e+01;
    LKconst.YWCst3.D3 = 0.15344e+02;
    LKconst.YWCst3.D4 = -0.3704e+02;
    LKconst.YWCst3.D5 = 0.20169e+02;

//YWCst5
    LKconst.YWCst5.A1 = 0.933e-01;
    LKconst.YWCst5.A2 = -0.3445e-00;
    LKconst.YWCst5.A3 = 0.4042e-00;
    LKconst.YWCst5.A4 = -0.2083e-00;
    LKconst.YWCst5.A5 = 0.5473e-01;
    LKconst.YWCst5.B1 = 0.22e-01;
    LKconst.YWCst5.B2 = -0.3363e-02;
    LKconst.YWCst5.B3 = -0.796e-01;
    LKconst.YWCst5.B4 = 0.8546e-01;
    LKconst.YWCst5.B5 = -0.217e-01;
    LKconst.YWCst5.C1 = 0.1937e-01;
    LKconst.YWCst5.C2 = -0.3055e-01;
    LKconst.YWCst5.C3 = 0.631e-01;
    LKconst.YWCst5.C4 = 0.0e-00;
    LKconst.YWCst5.D1 = -0.16e+02;
    LKconst.YWCst5.D2 = 0.30699e+02;
    LKconst.YWCst5.D3 = 0.19645e+02;
    LKconst.YWCst5.D4 = -0.81305e+02;
    LKconst.YWCst5.D5 = 0.47031e+02;

//YWCst9
    LKconst.YWCst9.A1 = -0.817e-01;
    LKconst.YWCst9.A2 = 0.3274e-00;
    LKconst.YWCst9.A3 = -0.5014e-00;
    LKconst.YWCst9.A4 = 0.387e-00;
    LKconst.YWCst9.A5 = -0.1342e-00;
    LKconst.YWCst9.B1 = -0.23e-01;
    LKconst.YWCst9.B2 = -0.124e-01;
    LKconst.YWCst9.B3 = 0.1625e-00;
    LKconst.YWCst9.B4 = -0.2135e-00;
    LKconst.YWCst9.B5 = 0.8643e-01;
    LKconst.YWCst9.C1 = 0.5626e-01;
    LKconst.YWCst9.C2 = -0.3518e-00;
    LKconst.YWCst9.C3 = 0.6194e-00;
    LKconst.YWCst9.C4 = -0.3809e-00;
    LKconst.YWCst9.D1 = -0.21e+02;
    LKconst.YWCst9.D2 = 0.55174e+02;
    LKconst.YWCst9.D3 = -0.33637e+02;
    LKconst.YWCst9.D4 = -0.28109e+02;
    LKconst.YWCst9.D5 = 0.26277e+02;

}


void LK4::Readfluid(const char* iFluidName){
    const char sep = '\t';
    const int nbTitleLines = 2; //nombre de lignes avec des titres
    const int nbTitleColumns = 1; //nombre de colonnes avec des titres ie sans chiffre ie
    const char* file_in = "dataBase/DataFluids.txt";
    TInputFile  InputFile((char*)file_in, nbTitleLines, nbTitleColumns, sep);
    InputFile.open();

    TOutputFile OutputFile((char*)"iofdata.out");
    OutputFile.open();
    int nbRecords = InputFile.nRecords();
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    Records = *InputFile.GetRecords();

    FluidName = (char*)iFluidName;
    TSequence* PtrSeq = Get(Records, FluidName);

    vc = PtrSeq->Get(0);
    Pc = PtrSeq->Get(1);
    kPc=Pc/1e3;
    Tc = PtrSeq->Get(2);
    Tnb = PtrSeq->Get(3);
    Mw = PtrSeq->Get(4);
    w = PtrSeq->Get(5);  //dit w dans PG
    DipM = PtrSeq->Get(6);

    // cpo coefficients
    ALPHA = PtrSeq->Get(7)/1000;
    BETA = PtrSeq->Get(8)/1000;
    GAMMA = PtrSeq->Get(9)/1000;
    DELTA = PtrSeq->Get(10)/1000;
    roL = PtrSeq->Get(11); //pas de 5eme coefficent cp dans le fichier d'entrée
    TL= PtrSeq->Get(12);
    Roc = PtrSeq->Get(13);
    Tfp = PtrSeq->Get(14);
    // reference state default initialisation
    // liquid line
    Tref = 273.15;	// [K]
    href = 200.0E3;	// [J/kg]
    sref = 1.0E3;	// [J/(kg K)]

}

void LK4::CriticProp (double & Vcrit,double & Pcrit,double & Tcrit,double & MM,double &Zcrit){
    Vcrit = vc;
    Pcrit = kPc * 1000;//  ==> N/m2 }
    Tcrit = Tc;
    Zcrit = Zc;
    MM = Mw;
}

void LK4::DomainProp (double &  Vmin, double & Vmax,double &  Pmin, double & Pmax, double & Tmin, double & Tmax){
    double Vcrit, Pcrit, Tcrit, MM, Zcrit;
    CriticProp(Vcrit, Pcrit, Tcrit, MM, Zcrit);
    Pmin = 0.0001 * Pcrit;
    Pmax = 10 * Pcrit;
    Tmin = 0.3 * Tcrit;
    Tmax = 4 * Tcrit;
    Vmax = v_PT_g(Pmin, Tmax);
    Vmin = v_PT_l(Pmax, Tmin);
}



int LK4::Zone_PT(double P,double T){
    double Vcrit, Pcrit, Tcrit, MM, Zcrit;
    CriticProp(Vcrit, Pcrit, Tcrit, MM, Zcrit);
    if (T > Tcrit)
    return 3;
    else {
        if(P <= Psat_T(T)) {
            if (P == Psat_T(T))
                return 1;
            else
                return 2;}

        else
            return 0;
    }
}


// %%%%%%%%%%%% Zone Gazeuse %%%%%%%%%%%%}

//------- input v,P -------}


void LK4::Prop_vP_g (const double & v,const double &  P,double& T,double& h,double& s) {
    iLK_Prop_vP_g(v, P / 1E3, T, h, s);
    h = h * 1E3; //{ ==> N/m2 }
    s = s * 1E3; //{ ==> J/kg K }
}
double LK4::T_vP_g (const double & v,const double &  P) {
    double T,h, s;
    Prop_vP_g(v, P, T, h, s);
    return T;
}
double LK4::h_vP_g (const double & v,const double &  P) {
    double  T,h, s;
    Prop_vP_g(v, P, T, h, s);
    return h;
}
double LK4::s_vP_g (const double & v,const double &  P) {
    double T, h,s;
    Prop_vP_g(v, P, T, h, s);
    return s;
}
double LK4::u_vP_g (const double & v,const double &  P) {
    double T, h, s;
    Prop_vP_g(v, P, T, h, s);
    return h - v * P;
}

//------- input v,T -------}

void LK4::Prop_vT_g (const double & v,const double &  T,double& P,double&h, double&s) {
    iLK_Prop_vT_g(v, T, P, h, s);
    P= P * 1E3;// ==> N/m2 }
    h= h * 1E3;// ==> J/kg }
    s= s * 1E3;// ==> J/kgK }
}
double LK4::P_vT_g (const double & v,const double &  T) {
    double VR, PR, TR;
        VR = v / vc * Zc;
    TR = T / Tc;
    Pr_vTr_g(VR, PR, TR);
    return PR * kPc * 1E3;
}
double LK4::h_vT_g (const double & v, const double & T) {
    return h_PT_g(P_vT_g(v, T), T);
}
double LK4::s_vT_g (const double & v,const double &  T) {
    return s_PT_g(P_vT_g(v, T), T);
}
double LK4::u_vT_g (const double & v,const double &  T) {
    double P;
    P= P_vT_g(v, T);
    return h_PT_g(P, T) - v * P;
}

//------- input v,h -------}

void LK4::Prop_vh_g (const double & v,const double &  h,double& P,double& T,double& s) {
    iiLK_Prop_vh_g(v, h / 1E3, P, T, s);
    P= P * 1E3;
    s= s * 1E3;
}
double LK4::P_vh_g (const double & v,const double &  h) {
    double P,T, s;
    Prop_vh_g(v, h, P, T, s);
    return P;
}
double LK4::T_vh_g (const double & v,const double &  h) {
    double P,T, s;
    Prop_vh_g(v, h, P, T, s);
    return T;
}
double LK4::s_vh_g (const double & v,const double &  h) {
    double P, T,s;
    Prop_vh_g(v, h, P, T, s);
    return s;
}
double LK4::u_vh_g (const double & v,const double &  h) {
    double P, T, s;
    Prop_vh_g(v, h, P, T, s);
    return h - v * P;
}

//------- input v,s -------}

void LK4::Prop_vs_g (const double & v,const double &  s,double& P,double& T,double& h) {
    iiLK_Prop_vs_g(v, s / 1E3, P, T, h);
    P= P * 1E3;
    h= h * 1E3;
}
double LK4::P_vs_g (const double & v, const double & s) {
    double P,T, h;
    Prop_vs_g(v, s, P, T, h);
    return P;
}
double LK4::T_vs_g (const double & v,const double &  s) {
    double P,T, h;
    Prop_vs_g(v, s, P, T, h);
    return T;
}
double LK4::h_vs_g (const double & v,const double &  s) {
    double P, T,h;
    Prop_vs_g(v, s, P, T, h);
    return h;
}
double LK4::u_vs_g (const double & v,const double &  s) {
    double P, T, h;
    Prop_vs_g(v, s, P, T, h);
    return h - v * P;
}


void LK4::GTH_PT_g (double P,double T,double& Cv,double& Cp,double& Sv,double&Av, double&Bp,double& Gt){
    //grandeurs thermiques diverses}
    double v, VR, VRo, VRr, Eo, Er, PR, TR, DPTo, DPTr, DPT, DPVo, DPVr, DPV;
    P= P / 1000; // ==> kN/m2 }
    PR= P / kPc;
    TR= T / Tc;
    VRro_PTr_g(PR, TR, VR, VRo, VRr, Eo, Er);
    v= VR * vc / Zc;
 //TODO convergence
    DPV_VTor(VRo, VRr, TR, DPVo, DPVr, DPV);
    DPT_VTor(VRo, VRr, TR, DPTo, DPTr, DPT);
    CVP_VT(VRo, VRr, Eo, Er, DPVo, DPVr, DPTo, DPTr, TR, Cv, Cp);
    Cv= Cv * 1E3;//=> J/kg K}
    Cp= Cp * 1E3;
    Sv= sqrt(-1000 * Cp / Cv * DPV) * v;
    Av= P / T / DPT;
    Gt= -P / v / DPV;
    Bp= Av / Gt;
}

void LK4::Sonic_PT_g (double Psup, double Tsup,double& vLaval,double&PLaval, double& TLaval,double& Csonic){
    //calcul la vitesse du son Csonic et la pression de Laval}

    LK_Sonic_PT_g(Psup / 1E3, Tsup, vLaval, PLaval, TLaval, Csonic);
    PLaval= PLaval * 1E3; // ==> N/m2 }
}
double LK4::Cv_PT_g (const double & P, const double & T){
    double Cv, Cp, Sv, Av, Bp, Gt;
    GTH_PT_g(P, T, Cv, Cp, Sv, Av, Bp, Gt);
    return Cv;
}
double LK4::Cp_PT_g (const double & P, const double & T){
    // Calcule la chaleur massique isobare en J/kg K en fonction de la pression en Pa  et de la *) temperature en K . *)
    double Cv, Cp, Sv, Av, Bp, Gt;

    GTH_PT_g(P, T, Cv, Cp, Sv, Av, Bp, Gt);
    return Cp;
}

void LK4::TransProp_PT_g (const double & P,const double &  T,double& ViscDyn, double& ViscCin, double& Ro, double&ThermCond,double& Prandtl) {
    //ViscDyn in [Pa s],ViscCin in [m2/s],ThermCond in [W/(m K)]}
    //from "The properties of gases and liquids" of Reid & Prausnitz, Chap. 9.4,9.6,10.5}
    LK_TransProp_PT_g(P / 1E3, T, ViscDyn, ViscCin, Ro, ThermCond, Prandtl);
}
double LK4::SurfacTension(const double & P,const double &  T) {
    // Surface tension correlation [N/m] from "The properties of gases and liquids" of Reid & Prausnitz, Chap. 12.3}
    double Pt, Q,Tbr;
    Pt= P;
    Tbr= Tnb/Tc;
    Q= 0.1196*(1.0+Tbr*ln(kPc/101.325)/(1.0-Tbr))-0.279;
    return 1.0E-3*pow(kPc*1.0E-2,2.0/3.0)*pow(Tc,1.0/3.0)*Q*pow(1.0-T/Tc,11.0/9.0);
}
//------- input P,h -------}

void LK4::Prop_Ph_g (const double & P,const double &  h,double& v,double& T,double& s) {
    iLK_Prop_Ph_g(P/ 1000, h/ 1000, v, T, s);
    s= s * 1000; // ==> J/kg K }
}
double LK4::v_Ph_g (const double & P,const double &  h) {
    double v,T, s;
    Prop_Ph_g(P, h, v, T, s);
    return v;
}
double LK4::T_Ph_g (const double & P,const double &  h) {
    double v,T, s;
    Prop_Ph_g(P, h, v, T, s);
    return T;
}
double LK4::s_Ph_g (const double & P,const double &  h) {
    double v, T,s;
    Prop_Ph_g(P, h, v, T, s);
    return s;
}
double LK4::u_Ph_g (const double & P,const double &  h) {
    double v, T, s;
    Prop_Ph_g(P, h, v, T, s);
    return  (h - P * v);
}

//------- input P,s -------}

void LK4::Prop_Ps_g (const double & P,const double &  s,double& v,double& T,double& h) {
    iLK_Prop_Ps_g(P/1.0E3, s/1.0E3, v, T, h);
    h= h *  1.0E3;// ==> J/kg }
}
double LK4::v_Ps_g (const double & P,const double &  s) {
    double v,T, h;
    Prop_Ps_g(P, s, v, T, h);
    return v;
}
double LK4::T_Ps_g (const double & P,const double &  s) {
    double v,T, h;
    Prop_Ps_g(P, s, v, T, h);
    return T;
}
double LK4::h_Ps_g (const double & P,const double &  s) {
    double v, T,h;
    Prop_Ps_g(P, s, v, T, h);
    return h;
}
double LK4::u_Ps_g (const double & P, const double & s) {
    double v, T, h;
    Prop_Ps_g(P, s, v, T, h);
    return h - v * P;
}

//------- input T,h -------}

void LK4::Prop_Th_g (const double & T,const double &  h,double& v,double& P,double& s) {
    iLK_Prop_Th_g(T, h / 1E3, v, P, s);
    P= P * 1E3;
    s= s * 1E3;
}
double LK4::v_Th_g (const double & T,const double &  h) {
    double v,P, s;
    Prop_Th_g(T, h, v, P, s);
    return v;
}
double LK4::P_Th_g (const double & T,const double &  h) {
    double v,P, s;
    Prop_Th_g(T, h, v, P, s);
    return P;
}
double LK4::s_Th_g (const double & T,const double &  h) {
    double v, P,s;
    Prop_Ts_g(T, h, v, P, s);
    return s;
}
double LK4::u_Th_g (const double & T,const double &  h) {
    double v, P, s;
    Prop_Th_g(T, h, v, P, s);
    return h - v * P;
}

//------- input T,s -------}

void LK4::Prop_Ts_g (const double & T, const double & s,double& v,double& P,double& h) {
    iLK_Prop_Ts_g(T, s / 1E3, v, P, h);
    P= P * 1E3;
    h= h * 1E3;
}
double LK4::v_Ts_g (const double & T, const double & s) {
    double v,P, h;
    Prop_Ts_g(T, s, v, P, h);
    return v;
}
double LK4::P_Ts_g (const double & T,const double &  s) {
    double v,P, h;
    Prop_Ts_g(T, s, v, P, h);
    return P;
}
double LK4::h_Ts_g (const double & T,const double &  s) {
    double v, P,h;
    Prop_Ts_g(T, s, v, P, h);
    return h;
}
double LK4::u_Ts_g (const double & T,const double &  s) {
    double v, P, h;
    Prop_Ts_g(T, s, v, P, h);
    return h - v * P;
}

//------- input T,u -------}

//------- input h,s -------}
void LK4::Prop_hs_g (const double & h,const double &  s,double& v,double& P,double& T) {
    iiLK_Prop_hs_g(h / 1E3, s / 1E3, v, P, T);
    P= P * 1E3;
}
double LK4::P_hs_g (const double & h,const double &  s) {
    double v, P, T;
    Prop_hs_g (h, s, v, P, T);
    return P;
}
//------- input s,u -------}
void LK4::Prop_su_g (const double & s,const double &  u,double& v, double&P,double& T,double& h) {
    iiLK_Prop_su_g(s / 1E3, u / 1E3, v, P, T, h);
    P= P * 1E3;// ==> N/m2 }
    h= h * 1E3;// ==> J/kg }
}
double LK4::v_su_g (const double & s,const double &  u) {
    double v,P, T, h;
    Prop_su_g(s, u, v, P, T, h);
    return v;
}
double LK4::P_su_g (const double & s,const double &  u) {
    double v,P, T, h;
    Prop_su_g(s, u, v, P, T, h);
    return P;
}
double LK4::T_su_g (const double & s,const double &  u) {
    double v, P,T, h;
    Prop_su_g(s, u, v, P, T, h);
    return T;
}
double LK4::h_su_g (const double & s, const double & u) {
    double v, P, T,h;
    Prop_su_g(s, u, v, P, T, h);
    return h;
}












void LK4::Prop_PT_g (const double & Pi, const double & T,double&  v, double& h, double& s){
    double dpv, dpt, sv,cp,cv;
    double P = Pi / 1000; // ==> kN/m2 }
    LK_Prop_PT_g(P, T, v, h, s, dpv, dpt, cv, cp, sv);
    h = h * 1000; // ==> J/kg }
    s = s * 1000; // ==> J/kg }
}

void LK4::PropE_PT_g (const double & Pi, const double & T,double& v, double& h, double& s,double&  Cv,double&  Cp){
    double dpv, dpt, sv,cp,cv;
    double P = Pi / 1000; // ==> kN/m2 }
    LK_Prop_PT_g(P, T, v, h, s, dpv, dpt, cv, cp, sv);
    h = h * 1000; // ==> J/kg }
    s = s * 1000; // ==> J/kg }
    cv = cv * 1000; //==> J/kg K}
    cp = cp * 1000; // ==> J/kg K}
}


double LK4::v_PT_g (const double & P,const double & T){
    double VR, PR, TR;
    PR = P / kPc / 1000; //  ==> kN/m2 }
    TR = T / Tc;
    Vr_PTr_g(VR, PR, TR);
    return VR * vc / Zc;
}


double LK4::h_PT_g (const double & Pi,const double & Ti){
    double VR, VRo, VRr, Eo, Er, H;
    double P = Pi / kPc / 1000; // { ==> kN/m2 }
    double T = Ti / Tc;
    VRro_PTr_g(P, T, VR, VRo, VRr, Eo, Er);
    // todo convergence
    H_VT(VRo, VRr, P, T, Eo, Er, H);
    return H * 1E3;
}


double LK4::s_PT_g (const double & Pi,const double & Ti){
    double VR, VRo, VRr, Eo, Er, s;
    double P = Pi / kPc ; // { ==> kN/m2 }
    double T = Ti / Tc;
    VRro_PTr_g(P, T, VR, VRo, VRr, Eo, Er);
    // todo convergence
    S_VT(VRo, VRr, T, Eo, Er, s);
    return s * 1E3;
}


double LK4::u_PT_g (const double & P,const double & T){
    return h_PT_g(P, T) - v_PT_g(P, T) * P;
}

double LK4::Z_PT_g (const double & P,const double & T){
    return P * v_PT_g(P, T) / T;
}


double LK4::Zr_PT_g (const double & P,const double & T){
    return P * v_PT_g(P, T) / (8314 / Mw) / T;
}


void LK4::Prop_P_s (const double &Psat, double& vl,double& vg,double &T, double &hl, double &hg, double &sl, double &sg){
    double DPTs, DPVs;
    double P_temp = Psat / 1000; // ==> kN/m2 }
    LK_Prop_P_s(P_temp, vl, vg, T, hl, hg, sl, sg, DPTs, DPVs);
    hl = hl * 1E3;// ==> J/kg }
    hg = hg * 1E3; // ==> J/kg }
    sl = sl * 1E3;// ==> J/kgK }
    sg = sg * 1E3; // ==> J/kgK }
}

double LK4::v_T_sl (double T){
    return LK_vs_T_sl(T);
}

double LK4::Psat_T (double T){
    double P;
    Ps_T(P, T);
    return P  * 1000; // { ==> N/m2 }
}


double LK4::hsat_T (double T){
    double hs;
    hs_T(hs, T);
    return hs  * 1000; // { ==> N/m2 }
}


double LK4::Tsat_h (double h){
    double Ts;
    Ts_h(Ts, h/1E3);
    return Ts;
}


double LK4::Tsat_v (double v){
    double Ts;
    Ts_v(Ts, v);
    return Ts;
}

double LK4::Tsat_P (double P){
    double Ts,DPT;
    Ts_P(P / 1E3, Ts, DPT);
    return Ts;
}

double LK4::dPsat_PT_sg (double Ps,double Ts){
    return DPTS_PT_sg(Ps / 1E3, Ts) * 1E3;
}


void LK4::Prop_Px_b (double P,double x, double & v,double & Tsat,double & h,double & s) {
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq;
    Prop_P_s(P, vGaz, vLiq, Tsat, hGaz, hLiq, sGaz, sLiq);
    v= (1 - x) * vLiq + x * vGaz;
    h= (1 - x) * hLiq + x * hGaz;
    s= (1 - x) * sLiq + x * sGaz;
}
void LK4::Prop_PTx_bg (double P, double &T,double & v,double & h,double & s,double & x) {
    //input:(P,T) ou (P,x)}

    if (T > Tsat_P(P)){
            Prop_PT_g(P, T, v, h, s);
            x= 1.0;}
    else	Prop_Px_b(P, x, v, T, h, s);
}

void LK4::Prop_vP_bg (double v,double P, double &T,double & h,double & s,double &x) {
    LK_Prop_vP_bg (v, P/ 1E3, T, h, s, x);
    h= h * 1E3;
    s= s * 1E3;}

void LK4::Prop_vT_bg (double v,double T, double& P,double & h, double &s,double &x) {
    LK_Prop_vT_bg (v, T, P, h, s, x);
    P= P * 1E3;
    h= h * 1E3;
    s= s * 1E3;
}


double LK4::x_Ph (double P,double h) {
    double vGaz, vLiq, Tsat, hGaz, hLiq, sGaz, sLiq;
    Prop_P_s(P, vLiq, vGaz, Tsat, hLiq, hGaz, sLiq, sGaz);
    if (h <= hLiq)
        return 0;
    else if ((hLiq < h) && (h < hGaz))
        return (h - hLiq) / (hGaz - hLiq);
    else
        return 1;
}

//%%%%%%%%%%%% Zone Liquide %%%%%%%%%%%%}

void LK4::GTH_PT_l (double P,double T,double &  Cv,double &  Cp,double &   Sv,double &   Av,double &   Bp, double &  Gt){
//grandeurs thermiques diverses}
    double v, h, s, CpSat,
    DPV, DPT;
    PropE_PT_l (P, T, v, h, s, Cv, Cp, CpSat, DPV, DPT);
    Sv= sqrt(-1000 * Cp / Cv * DPV) * v;
    Av= P / T / DPT;
    Gt= -P / v / DPV;
    Bp= Av / Gt;
}


void LK4::Prop_PT_l (const double & P,const double & T, double &v, double & h, double & s) {
    double DPV, DPT, Cv, Cp, CpSat;
    PropE_PT_l (P, T,v, h, s, Cv, Cp, CpSat, DPV, DPT);
}

void LK4::PropE_PT_l (double P,double T, double &v, double & h, double & s,double &  Cv,double &  Cp,double &  CpSat,double &  DPV,double &  DPT){
    double RO, CPL, CVL, CPF;
    P= P / 1000; // ==> kN/m2 }
    LK_Prop_PT_l(P, T, v, RO, h, s, DPV, DPT);
    h= h * 1000;// ==> J/kg }
    s= s * 1000;// ==> J/kg }
    LK_Cvp_PT_l(P, T, CPL, CVL, CPF);
    Cv= CVL * 1000;// ==> J/kg K}
    Cp= CPL * 1000;// ==> J/kg K}
    CpSat= CPF * 1000;// ==> J/kg K}
}

double LK4::v_PT_l (const double & P,const double & T) {

    cout<<"LK4"<<endl;
    double Ptemp=P/ 1e3;
    return LK_v_PT_l(Ptemp , T);
}
double LK4::h_PT_l (const double & Pi,const double & T) {
    double  RO, DPV, DPT, v, s, h;
    double P= Pi / 1000; // ==> kN/m2 }
    LK_Prop_PT_l(P, T, v, RO, h, s, DPV, DPT);
    return h * 1000;// ==> J/kg }
}
double LK4::T_Ph_l (const double & Pi,const double & hi) {
    double T;
    double P= Pi / 1000; // ==> kN/m2 }
    double h= hi / 1000; // ==> kJ/kg }
    iLK_T_Ph_l(P, h, T);
    return T;
}
double LK4::v_Ph_l (const double & P,const double & h) {
    return v_PT_l(P, T_Ph_l(P, h));
}
void LK4::Th_Ps_l (const double & Pi,const double & si,  double & T, double &h) {
    double P= Pi / 1000; // ==> kN/m2 }
    double s= si / 1000; // ==> kJ/kg }
    iLK_Th_Ps_l(T, h, P, s);
    h= h * 1000;// ==> J/kg }
}
void LK4::TransProp_PT_l (const double & P,const double & T,  double & ViscDyn,  double &ViscCin,  double &Ro,  double &ThermCond, double & Prandtl) {
    //ViscDyn in [Pa s],ViscCin in [m2/s],ThermCond in [W/(m K)]}
    //from "The properties of gases and liquids" of Reid & Prausnitz, Chap. 9.4,9.6,10.5}
    cout<<"LK4:TransProp_PT_l"<<endl;
    LK_TransProp_PT_l(P / 1E3, T, ViscDyn, ViscCin, Ro, ThermCond, Prandtl);
}








void LK2::Vr_PTr_g(double &VR, double &PR,double & TR){
    double VRo, VRr;
    iVor_PTr_g(VRo, PR, TR, CstO);
    iVor_PTr_g(VRr, PR, TR, CstR);
    VR = (1.0 - AF) * VRo + AF * VRr;
}

void LK2::Pr_vTr_g(double &VR, double &PR,double & TR){
    double P, Pi, vi, T,DPVo, DPVr, dPv, v, VRo, VRr, Eo, Er;
    int i;
    // initialisation}
    i = 0;
    T = TR * Tc;
    v = VR * vc / Zc;
    if (VR < 2)
        Ps_T(P, T);
    else
        P = TR / VR * kPc;
    do{
    Pi = P;
    VRro_PTr_g(P / kPc, TR, VR, VRo, VRr, Eo, Er);
    DPV_VTor(VRo, VRr, TR, DPVo, DPVr, dPv);
    if (i>100)    {cout<<"error Pr_vTr_g "<<endl;
        exit(EXIT_FAILURE);
    }
    i++;
    vi = VR * vc / Zc;
    P = Pi - (vi - v) * dPv;
    }until((abs(P - Pi) < 1E-4) && (abs(v - vi) < 1E-6));
    PR = P / kPc;
}


void LK3::Pr_vTr_l(double &VR, double &PR,double & TR){
    double P, Pi, vi, T,DPVo, DPVr, dPv, v, VRo, VRr, Eo, Er;
    int i;
    // initialisation}
    i = 0;
    T = TR * Tc;
    v = VR * vc / Zc;
    if (VR < 2)
        Ps_T(P, T);
    else
        P = TR / VR * kPc;
    do{
        Pi = P;
        VRro_PTr_l(P / kPc, TR, VR, VRo, VRr);
        DPV_VTor(VRo, VRr, TR, DPVo, DPVr, dPv);
        if (i>100)   {cout<<"error Pr_vTr_l "<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
        vi = VR * vc / Zc;
        P = Pi - (vi - v) * dPv;
    }until((abs(P - Pi) < 1E-4) && (abs(v - vi) < 1E-6));
    PR = P / kPc;
}


void LK2::hs_T(double &hs, double T){
    double TR, P, PR, VR, VRo, VRr, Eo, Er;
    Ps_T(P, T);
    PR = P / kPc;
    TR = T / Tc;
    VRro_PTr_g(PR, TR, VR, VRo, VRr, Eo, Er);
    H_VT(VRo, VRr, PR, TR, Eo, Er, hs);
}

void  LK2::Ts_h(double &Ts,double h){
    double m, v, P, newT, Tinit, newh, S, DPV, DPT, DPTS, CV, CP, SV;
    int i;
    newT = Tc * 0.6;
    i = 0;
    do {
        Tinit = newT;
//-- -- -calcul de h, et de la pente-- -- -
        Ps_T(P, newT);
        LK_Prop_PT_g(P, newT, v, newh, S, DPV, DPT, CV, CP, SV);
        DPTS = DPTS_PT_sg(P, newT);
        if (i>50){cout<<"error Ts_h "<<endl;
            exit(EXIT_FAILURE);
        }
        i++;
        m = (v + newT * DPT / DPV) * DPTS + CP;
//-- -- -fin de ce calcul-- -- --
        newT = Tinit - (newh - h) / m;
    }until(abs(newT - Tinit) < 0.1);
    Ts = newT;

}




double LK3::LK_vs_T_sl(double T){
    double VR, VG, PR, TR, PS;
    Ps_T(PS, T);
    PR = PS / kPc;
    TR = T / Tc;
    Vr_PTr_g(VR, PR, TR);
    VG = VR * vc / Zc;
    return 1 / (ADA1 * T + ADA0 - 1 / VG);
}