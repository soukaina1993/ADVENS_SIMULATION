//
// Created by Yolaine Adihou on 17/04/2021.
//

#include "MWater.h"
using namespace std;
//--------------------------------------------------------
#pragma mark		=== TOOL FUNCTIONS ===
//--------------------------------------------------------

void MWater::init(){
    // k = 0.0;
    // vc = v_PT_g(Pc,Tc);
    r = R*1e3/Mw;
    Zc = Pc*vc/r/Tc;

// set FluidAprox class
    Set(Pc,Tc,Zc,w);
    double Ps = P_T_s(Tref);
    h_cst = href-h_PT_l(Ps,Tref);
    s_cst = sref-s_PT_l(Ps,Tref);

    double P;
    Prop_vT(vc,Tc,P,uc,hc,sc);
}

#pragma mark		=== BASE FUNCTIONS ===

void MWater::CvCp_T_o(const double &T, double &cvo, double &cpo){
    CvCp_vT_g(T,cvo,cpo);
}
double MWater::h_T_o(const double &T){
    return T*(c0+T*(0.5*c1+T*(c2/3.0+0.25*c3*T)))-c4/T + h_cst;
}
double MWater::s_T_o(const double &T){
    return T*(c1+T*(0.5*c2+T*c3/3.0))+c0*ln(T)-0.5*c4/sqr(T) + s_cst;
}

void MWater::IdealState(const double &T,
                    double &cpo,
                    double &cvo,
                    double &ho,
                    double &so){
    CvCp_T_o(T,cvo,cpo);
    ho = h_T_o(T);
    so = s_T_o(T);
}

// EOS partial diff.
void MWater::dPvdPT_vT(	const double &v,const double &T,double &dPdv,double &dPdT){
    double Tx, f1, f2, P;
    Tx = TX(T);
    f1 = F1(Tx);
    f2 = F2(Tx);
    P = P_f1(T, f1);
    dPdT = DPDT(P, T, f1, f2);
    dPdv = DPDV(v,T);
    //P1R=P1R0*ROT
}

// thermal factors
void	MWater::AvBpGt_vT(
        const double &v,
        const double &T,
        double &Av,
        double &Bp,
        double &Gt){
    double dPdT, dPdv, P;
    dPvdPT_vT(v,T,dPdv,dPdT);
    P = P_vT(v,T);
    Av = P / T / dPdT;
    Gt = -P / v / dPdv;
    Bp = Av / Gt;
}

MWater::MWater(){
    vc=3.1545741e-3;    //vc=0.00310559;
    Pc=22.089e6;        //Pc=22064000.0;
    Tc=647.286;         //Tc=647.096;
}

void	MWater::CvCp_vT(const double &v,
                   const double &T,
                   double &cv,
                   double &cp){
    double vl, vg, Psat, ul, ug, hl, hg, sl, sg;
    Prop_T_s(T, vl, vg, Psat, ul, ug, hl, hg, sl, sg);

    if (v>vg)       //cornelia.blanke: changed from (v>vl)
        CvCp_vT_g(v,T,cv,cp);
    else if (v<vl)
        CvCp_vT_l(v,T,cv,cp);
    else{
        double cv_g, cp_g,cv_l,cp_l,x;
        x = x_fraction(v,vl,vg);
        CvCp_vT_g(vg,T,cv_g,cp_g);
        CvCp_vT_l(vl,T,cv_l,cp_l);
        cv=cv_g*x+cv_l*(1-x);
        cp=cp_g*x+cp_l*(1-x);
        //fluid.TransProp_PT_g(P,T,ViscDyn,ViscCin,Ro,ThermCond,Prandtl);
    }
}


void MWater::CvCp_vT_g(const double &v,const double &T,double &cv,double &cp){
    CvCp_vT_g(T, cv, cp);
}

void MWater::CvCp_vT_l(const double &v,const double &T,double &cv,double &cp){
    CvCp_vT_l(T, cv, cp);
}

void MWater::CvCp_vT_g(const double &T,double &cv,double &cp){
    cp = Cp_T_g(T);
    cv = cp - r;
}

void MWater::CvCp_vT_l(const double &T,double &cv,double &cp){
    cp = Cp_T_l(T);
    cv = cp;
}

double	MWater::P_vT(const double &v, const double &T){
    double P,u,h,s;
    Prop_vT(v,T,P, u, h,s);
    return P;
}
double	MWater::T_vP(const double &v, const double &P){
    return T_vP_a(v,P);
}

void MWater::v_PT(	const double &P, const double &T, double &vl, double &vg){
    vg = v_PT_g(P,T);
    vl = v_PT_l(P,T);
}

double	MWater::u_vT(	const double &v,
                   const double &T){
    double P,u,h,s;
    Prop_vT(v,T,P, u, h,s);
    return u;
}
double	MWater::s_vT(	const double &v,
                   const double &T){
    double P,u,h,s;
    Prop_vT(v,T,P, u, h,s);
    return s;
}
double	MWater::h_vT(	const double &v,
                   const double &T){
    double P,u,h,s;
    Prop_vT(v,T,P, u, h,s);
    return h;
}

double	MWater::h_vPT(	const double &v,
                    const double &P,
                    const double &T){
    double u,h,s;
    double Psat;
    Prop_vT(v,T,Psat, u, h,s);
    return h;
}
void	MWater::Prop_vT(const double &v,
                const double &T,
                double &P,
                double &u,
                double &h,
                double &s){
    double x;
    Prop_vT_a(v, T, P, u, h, s, x);
}

//    void W_Prop_vT_bg ( double& v, double& T,double& P,double& h,double& s,double& x);
//------------------------------------
// SATURATION LINE FUNCTIONS


void MWater::Prop_P_s(const double &P, double& vl,double&vg,double &T,double& ul,double &ug, double &hl, double &hg, double &sl, double &sg){
    double dpdT, Psat;
    T_P_s(P, T, dpdT);
    //Psat=P;
    Prop_T_s(T,vl,vg,Psat,ul,ug,hl,hg,sl,sg);
}

void MWater::Prop_T_s(const double &T, double&vl,double &vg,double &P,double &ul,double &ug, double &hl, double &hg, double &sl, double &sg) {
    double dpdT, dpv, P1;
    P_T_s(T, P, dpdT);
    v_PT(P,T,vl,vg);
    PWater_vT(vg, T, P1, ug, hg, sg, dpv);
    PWater_vT(vl, T, P1, ul, hl, sl, dpv);
}

void MWater::Prop_PT_l(const double &P,const  double &T,double &v, double &u, double &h, double &s) {
    v = v_PT_l(P,T);
    double Psat;
    Prop_vT(v,T,Psat, u, h,s);
}

double MWater::v_PT_l (const double& P,const double& T)// Calculation of a start value for v in liquid state
{
    return v_T_l(T);
}

double MWater::v_T_l (const double& T)// Calculation of a start value for v in liquid state
{
    double  Df;
    DWater(T, Df);
    return 1/Df;
}

void MWater::Prop_PT_g(const double &P,const  double &T,double &v, double &u, double &h, double &s) {
    v = v_PT_g(P,T);
    double P1;
    Prop_vT(v,T,P1, u, h,s);
}

void MWater::Prop_PT_g(const double &P,const double &T, double &h, double &s) {
    double u, v;
    Prop_PT_g(P, T, v, u, h, s);
}

//void MWater::Prop_vT_g(const double &v, const double &T, double &P, double &u, double &h, double &s) {
//    double dpv;
//    PWater_vT(v,T,P, u, h,s, dpv);
//}

void MWater::Prop_PT_g(double &P, double &T, double &v,double& u, double &h, double &s, double &DPV, double &DPT, double &CV,
                        double &CP, double &SV) {
    double P1;
    Prop_PT_g (P, T,v,u,h,s);
    CP = Cp_T_g(T);
    CV = CP - r;
    //DPVWater (v, P, DPV);
    PWater_vT(v, T, P, u, h, s, DPV);
    SWater (T, P1, DPT);
    SV = sqrt(- CP / CV * DPV) * v;
}

void MWater::TransProp_PT_g(const double &P,const double &T, double &ViscDyn, double &ViscCin, double &Rho, double &ThermCond,
                             double &Prandtl) {
    double v = v_PT_g(P, T);
    TransProp_vT_g(v, T, ViscDyn, ViscCin, Rho, ThermCond, Prandtl);
}

void MWater::FlowProp_vT_g(const double &v,const double &T,double &cp,double &ViscDyn,double &Rho,double &ThermCond,double &Prandtl){
    Rho = 1/v;
    cp = Cp_T_g(T);
    ViscDyn = ViscDyn_T_g(T);
    ThermCond = ThermCond_T_g(T);
    Prandtl = ViscDyn * cp / ThermCond;
}
void MWater::TransProp_vT_g(const double &v,const double &T, double &ViscDyn, double &ViscCin, double &Rho, double &ThermCond,
                           double &Prandtl) {
    double CP;
    FlowProp_vT_g(v, T, CP, ViscDyn, Rho, ThermCond, Prandtl);
    ViscCin = ViscDyn / Rho;
}

double MWater::TX(const double &T) {
    return 0.01 * (T - 338.15);
}
double MWater::F1(const double& Tx) {
    return -7.4192420 + Tx * (2.9721000e-1 + Tx * (-1.1552860e-1 + Tx * (8.6856350e-3 + Tx *
           (1.0940980e-3 + Tx * (-4.3999300e-3 + Tx * (2.5206580e-3 - Tx * 5.2186840e-4))))));
}
double MWater::F2(const double& Tx) {
return 2.9721000e-1 + Tx * (-1.1552860e-1 * 2 + Tx * (8.6856350e-3 * 3 + Tx * (1.0940980e-3 * 4 +
        Tx * (-4.3999300e-3 * 5 + Tx * (2.5206580e-3 * 6 - Tx * 5.2186840e-4 * 7)))));
}
double MWater::P_f1(const double& T, const double& f1){
    return Pc * exp((Tc / T - 1) * f1);
}
double MWater::DPDT(const double &P, const double& T, const double& f1, const double& f2){
    return P * ((-Tc) / T / T * f1 + (Tc / T - 1) * 0.01 * f2);
}

double MWater::Rho_PT_g(const double &P, const double &T) {
    //return 222 / T * P / 101325;        //ideal gas law not generally applicable!
    return 1.0/v_PT_g(P, T);
}
double MWater::Cp_T_g(const double &T) {
    return 2503.2 - 2.5183 * T + 0.0040416 * sqr(T) - 2.1642e-6 * XpwrY(T, 3) + 0.4092e-9 * XpwrY(T, 4);
}
double MWater::ViscDyn_T_g(const double &T) {
    return (-524.17 + 38.433 * T - 3.674e-3 * sqr(T)) * 1e-9;
}
double MWater::ThermCond_T_g(const double &T) {
    return exp(-4.6762 + 2.8947e-3 * T - 0.68587e-6 * sqr(T));
}

void MWater::TransProp_PT_l(const double &P,const double &T, double &ViscDyn, double &ViscCin, double &Rho, double &ThermCond,
                             double &Prandtl) {
    double RHO;
    Rho = Rho_T_l(T);
    TransProp_vT_l(1/Rho, T, ViscDyn, ViscCin, RHO, ThermCond, Prandtl);
}
void MWater::FlowProp_vT_l(const double &v,const double &T,double &cp,double &ViscDyn,double &Rho,double &ThermCond,double &Prandtl){
    Rho = 1.0/v;
    cp = Cp_T_l(T);
    ViscDyn = ViscDyn_T_l(T);
    ThermCond = ThermCond_T_l(T);
    Prandtl = ViscDyn * cp / ThermCond;
}
void MWater::TransProp_vT_l(const double &v,const double &T, double &ViscDyn, double &ViscCin, double &Rho, double &ThermCond,double &Prandtl) {
    double C;
    FlowProp_vT_l(v, T, C, ViscDyn, Rho, ThermCond, Prandtl);
    ViscCin = ViscDyn / Rho;
}

double MWater::Rho_T_l(const double &T) {
    double rho;
    DWater(T, rho);
    return rho;
    //return -43.36 + 10.5 * T - 37.5e-3 * sqr(T) + 57.34e-6 * XpwrY(T, 3) - 34.65e-9 * XpwrY(T, 4);
}
double MWater::Cp_T_l(const double &T) {
    return 31711 - 286.8 * T + 1.105 * sqr(T) - 1.874e-3 * XpwrY(T, 3) + 1.188e-6 * XpwrY(T, 4);
}
double MWater::ViscDyn_T_l(const double &T) {
    return (7003 - 39.7 * T + 0.0768 * sqr(T) - 49.83e-6 * XpwrY(T, 3)) * 1e-6;
}
double MWater::ThermCond_T_l(const double &T) {
    return -0.5278 + 6.246e-3 * T - 8.888e-6 * sqr(T) + 2.169e-9 * XpwrY(T, 3);
}



void MWater::Prop_vs_g(const double &v, const double &s, double &P, double &T, double &u,double &h) {
    double Pi, si, m, vi;
    double dpv, dpt, cv, cp, sv;
    int count = 0;
    T=T_v_s(v);
    SWater(T, P, dpt);
    do {
        Pi = P;
        //iteration pour trouver s=f(v,P)
        do {
            Prop_PT_g(P, T, vi,u, h, si, dpv, dpt, cv, cp, sv);
            //if IterFailed(i, 100, W_Prop_vs_gFailed)
            //then exit(W_Prop_vs_g);
            T += (vi - v) * dpv / dpt;
            count++;
        } while (abs(vi - v) >= 1E-6 && count <=100);
                //fin de l'it�ration intermediaire
        m = cv / T / dpt;
        P = Pi - (si - s) / m;
    } while( (abs(P - Pi) >= 1E0) && (abs(si - s) >= 1E-5));
    Prop_PT_g(P, T, vi,u, h, si, dpv, dpt, cv, cp, sv);
}

double MWater::T_v_s(const double &v) {
    double Ti, vi;
    double P, h, s, u, dpv, dpt, dT;
    int count = 0;

    //initialisation
    Ti = Tc * 0.4;
    //iteration}
    do {
        SWater(Ti, P, dpt);
        vi = v_PT_g(P,Ti);
        PWater_vT(vi, Ti, P, u, h, s, dpv);
        //if IterFailed(i, 100, W_Ts_vFailed)
        //then exit(W_Ts_v);
        dT = (v - vi) * dpv / dpt;
        Ti += dT;
        count++;
    } while ((abs(dT) / Ti >= 1E-3) && (abs(vi - v) / v >= 1E-3) && count<=maxiter);
    if (count > maxiter)
        cerr << "T_v_s not converged" << endl;
    return Ti;
}


double MWater::v_PT_g (const double& P1,const double& T1) // Calculation of a start value for v in gas state
{   double T,P;
    T = T1 / 100;
    P = P1 * 10.1972e-2;
    return 4706 * T / P - 0.668 / XpwrY(T, 2.7) - (437 + 6.31e-3 * P - 1.81e-29 * XpwrY(P, 5)) / XpwrY(T, 8.4)
            - 1.6e-10 * XpwrY(P, 5) / XpwrY(T, 30.5) - 3.26e-43 * XpwrY(P, 25) / XpwrY(T, 147);
}



void MWater::PWater_vT(const double&V, const double&T, double&P, double&U, double&H, double& S, double& DPV)
{
    // TODO strange behaviour close to saturation line
    double D, Q, QR, QT, X, Y, Z, TA, T1, B, A1, A2, A3, A4, A5, A6, A7, A, R1, R2, R3, R4, R5, R6, R7;

    D = 1 / V;
    TA = 1000 / T;            // =Ta/T
    X = D - 1000;            // X=D-D(AJ) where D(AJ) =1000 for J >1
    Y = D - 634;            // Y=D-D(A1) where D(A1)=634
    Z = exp(-4.8e-3 * D);    // Z=EXP(-E*D)
    T1 = TA - 1000/Tc;        // T1=TA-TA(AJ) where TA(AJ)=1.544912=Ta/Tc=1000/647.286 for J=1
    B = TA - 2.5;            // B=TA-TA(AJ) where TA(AJ)=2.5 for  J >1

    A1 = 2.9492937e-2 + (-1.3213917e-4 + (2.7464632e-7 + (-3.6093828e-10 + (3.4218431e-13 + (-2.4450042e-16 +
                                                                                             (1.5518535e-19 +
                                                                                              5.9728487e-24 * Y) * Y) *
                                                                                            Y) * Y) * Y) * Y) * Y +
         Z * (-4.1030848e-1 - 4.1605860e-4 * D);
    A2 = -5.1985860e-3 + (7.7779182e-6 + (-3.3301902e-8 + (-1.6254622e-11 + (-1.7731074e-13 + (1.2748742e-16 +
                                                                                               (1.3746153e-19 +
                                                                                                1.5597836e-22 * X) *
                                                                                               X) * X) * X) * X) * X) *
                         X + Z * (3.3731180e-1 - 2.0988866e-4 * D);
    A3 = 6.8335354e-3 + (-2.6149751e-5 + (6.5326396e-8 - 2.6181978e-11 * X) * X) * X +
         Z * (-1.3746618e-1 - 7.3396848e-4 * D);
    A4 = -1.5641040e-4 + (-7.2546108e-7 + (-9.2734289e-9 + 4.3125840e-12 * X) * X) * X +
         Z * (6.7874983e-3 + 1.0401717e-5 * D);
    A5 = -6.3972405e-3 + (2.6409282e-5 + (-4.7740374e-8 + 5.6323130e-11 * X) * X) * X +
         Z * (1.3687317e-1 + 6.4581880e-4 * D);
    A6 = -3.9661401e-3 + (1.5453061e-5 + (-2.9142470e-8 + 2.9568796e-11 * X) * X) * X +
         Z * (7.9847970e-2 + 3.9917570e-4 * D);
    A7 = -6.9048554e-4 + (2.7407416e-6 + (-5.1028070e-9 + 3.9636085e-12 * X) * X) * X +
         Z * (1.3041253e-2 + 7.1531353e-5 * D);
    A = A2 + B * (A3 + B * (A4 + B * (A5 + B * (A6 + B * A7))));
    Q = A1 + T1 * A;
    QT = -TA / T * (A + T1 * (A3 + B * (2 * A4 + B * (3 * A5 + B * (4 * A6 + B * 5 * A7)))));

    R1 = -1.3213917e-4 + (2.7464632e-7 * 2 + (-3.6093828e-10 * 3 + (3.4218431e-13 * 4 + (-2.4450042e-16 * 5 +
                                                                                         (1.5518535e-19 * 6 +
                                                                                          5.9728487e-24 * 7 * Y) * Y) *
                                                                                        Y) * Y) * Y) * Y +
         Z * (-4.1605860e-4 - 4.8e-3 * (-4.1030848e-1 - 4.1605860e-4 * D));
    R2 = 7.7779182e-6 + (-3.3301902e-8 * 2 + (-1.6254622e-11 * 3 + (-1.7731074e-13 * 4 + (1.2748742e-16 * 5 +
                                                                                          (1.3746153e-19 * 6 +
                                                                                           1.5597836e-22 * 7 * X) * X) *
                                                                                         X) * X) * X) * X +
         Z * (-2.0988866e-4 - 4.8e-3 * (3.3731180e-1 - 2.0988866e-4 * D));
    R3 = -2.6149751e-5 + (6.5326396e-8 * 2 - 2.6181978e-11 * 3 * X) * X +
         Z * (-7.3396848e-4 - 4.8e-3 * (-1.3746618e-1 - 7.3396848e-4 * D));
    R4 = -7.2546108e-7 + (-9.2734289e-9 * 2 + 4.3125840e-12 * 3 * X) * X +
         Z * (1.0401717e-5 - 4.8e-3 * (6.7874983e-3 + 1.0401717e-5 * D));
    R5 = 2.6409282e-5 + (-4.7740374e-8 * 2 + 5.6323130e-11 * 3 * X) * X +
         Z * (6.4581880e-4 - 4.8e-3 * (1.3687317e-1 + 6.4581880e-4 * D));
    R6 = 1.5453061e-5 + (-2.9142470e-8 * 2 + 2.9568796e-11 * 3 * X) * X +
         Z * (3.9917570e-4 - 4.8e-3 * (7.9847970e-2 + 3.9917570e-4 * D));
    R7 = 2.7407416e-6 + (-5.1028070e-9 * 2 + 3.9636085e-12 * 3 * X) * X +
         Z * (7.1531353e-5 - 4.8e-3 * (1.3041253e-2 + 7.1531353e-5 * D));
    QR = R1 + T1 * (R2 + B * (R3 + B * (R4 + B * (R5 + B * (R6 + B * R7)))));

    P = D * 461.51 * T * (1 + D * (Q + D * QR));
//    if (P < 0)
//        cerr << "Negative Pressure in PWater_vT" << endl;
    U = 4.6e4 * ln(T) +
        (1.011246e3 + (8.389300e-1 / 2 + (-2.199890e-4 / 3 + (2.466190e-7 / 4 - 9.704700e-11 / 5 * T) * T) * T) * T) *
        T - 564412.7036970831 - 461.51 * T * T * D * QT + 2.3750207e6;
    H = U + P / D;
    S = -4.6e4 / T + 1.011246e3 * ln(T) +
        (8.389300e-1 + (-2.199890e-4 / 2 + (2.466190e-7 / 3 - 9.704700e-11 / 4 * T) * T) * T) * T - 5727.2441340349 -
        461.51 * ln(D) - 461.51 * D * (Q + T * QT) + 6.6965776e3;
    //   DPV = -461.51 * T * D * D * (1 + 2 * Q * D + 2 * QR * D * D); // faux ( a refaire)
    DPV = R * T *(1 + 2 * Q * vc/V + 2 * QR * vc/V * vc/V); // faux ( a refaire)
}

//-----------------------------------------------------------------------------------------------------------------------------
void MWater::SWater (const double& T,double& P,double& dPdT)		// Saturation pressure P and DPDT from T
{
    double Tx, f1, f2;

    if (T > Tc)
        P = Pc;
    else {
        Tx = TX(T);
        f1 = F1(Tx);
        f2 = F2(Tx);
        P = P_f1(T, f1);
        dPdT = DPDT(P, T, f1, f2);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
void MWater::DWater (const double& T,double& Df) // Density of Saturated Liquid from T
{    double a = exp(ln(1 - T / Tc) / 3.0);
    Df = 317.0 * (1 + (3.6711257 + (-2.8512396e1 + (2.2265240e2 + (-8.8243852e2 + (2.0002765e3 + (-2.6122557e3 +
                                                                                                  (1.8297674e3 -
                                                                                                   5.3350520e2 * a) *
                                                                                                  a) * a) * a) * a) *
                                                   a) * a) * a);
}

//-----------------------------------------------------------------------------------------------------------------------------
void MWater::SatWater(double& T,double& P,double& dpdT,int option) {
    switch (option) {
        case 1: //Saturation pressure p and dpdT at T
            SWater(T, P, dpdT);
            break;
        case 2:  //Saturation temperature T and dpdT at p
            //if (P > 20.089e6)         error in Pascal code ?!
            if (P > Pc)
                T = Tc;
            else {
                double px, dT;
                int count = 1;
                const double err = P_precision * P;
                T = Tc / (1 - ln(P / Pc) / 7.4192420); // Start value

                if (T > Tc)
                    T = Tc - 0.001;
                SWater(T, px, dpdT);

                while ((abs(P - px) >= err) && (count <= maxiter)) {
                    dT = (P - px) / dpdT;
                    if (abs(dT) > 0.1 * T)
                        dT *= 0.1 * T / abs(dT);
                    T += dT;
                    if (T > Tc)
                        T = Tc - 0.001;
                    SWater(T, px, dpdT);
                    count++;
                }
                if (count > maxiter)
                    cerr<< "SatWater not convergent : p, T = "<< P << ", " << T<<endl;
            }
    }
}

double MWater::T_P_s (const double &P) {
    double T, dpdT;
    SatWater(T, const_cast<double &>(P), dpdT, 2);
    return T;
}

double MWater::P_T_s (const double &T){
    if (T > Tc)
        return Pc;
    else {
        double Tx = TX(T);
        double f1 = F1(Tx);
        return P_f1(T, f1);
    }
}

void MWater::P_T_s (const double& T,double& P,double& DPDT){
    SWater (T, P,DPDT);
}

void MWater::T_P_s (const double& P,double& T,double& DPDT){
    SatWater(T, const_cast<double &>(P), DPDT, 2);
}
//-----------------------------------------------------------------------------------------------------------------------------
/*void MWater::SaturateWater (double& T,double& p,double& vg,double& ug,double& hg,double& sg,double& vf,double& uf,double& hf,double& sf, int option) {

    double dpdT, Df;

    SatWater(T, p, dpdT, option);
    vg = v_PT_g(T, p);
    PropertyWater(T, p, vg, ug, hg, sg, option); // Route 1
    DWater(T, Df);
    vf = 1 / Df;
    uf = ug + (vg - vf) * (p - T * dpdT);
    hf = hg - T * (vg - vf) * dpdT;
    sf = sg - (vg - vf) * dpdT;
    //LiquidWater(T, p, vf, uf, hf, sf, 1);// Dangerous
};*/



double MWater::T_vP_a (const double &v, const double& P){
    //TODO what about density anomaly?
    double px,u,h,s;
    double vl, vg, ul, ug, hl, hg, sl, sg;
    int count=1, last_change=1;
    double T, Ttemp, perr, perr_new, dT;
    //T = T_P_s(P);      // saturation temperature at P
    Prop_P_s(P, vl, vg, T, ul, ug, hl, hg, sl, sg);

    if ((vl<=v) && (v<=vg))   //biphasic
        return T;

    Prop_vT(v,T,px, u, h,s);
    perr = p_error(P,px);
    dT = perr;
    while ((perr >= erp) && (count <= maxiter)) {
        if (px > P) {
            Ttemp = T - dT;        // cornelia.blanke: changed formula
        }
        else {
            Ttemp = T + dT;
        }
        Prop_vT(v, Ttemp, px, u, h, s);
        if (((perr_new=p_error(P, px)) > perr) && (dT > 1e-30)) {    // error growing, reduce dT
            dT /= 2.0;
            last_change = count;
            //cout << "Step size reduced" << endl;
        }
        else {
            T = Ttemp;                                      // accept Ttemp, update dT
            dT *=  perr_new/perr;
            perr = perr_new;
            //cout << count << ": " << px << ", " << T << endl;
            count++;
            if (count - last_change > 10) {    // use bigger dT
                dT *= 2.0;
                last_change = count;
                //cout << "Step size enlarged" << endl;
            }
        }
    }
    if (count >maxiter)
        cerr << "T_vP_a not convergent : p, T = "<< px <<", " <<T<<endl;
    return T;
}

double MWater::p_error(const double& p, const double& px){
    return abs(p - px)/p;
}
bool MWater::p_not_converged(const double& p, const double& px){
    return abs(p - px) >= erp * p;
}
bool MWater::s_not_converged(const double& s, const double& sx){
    return abs(s - sx) >= ers * 461.51;
}
bool MWater::h_not_converged(const double& h, const double& hx, const double& T){
    return abs(h - hx) >= erh * 461.51 * T;
}

double MWater::alter_v(const double& v, const double& vmin, const double& vmax){
    double dv;
    if (v <= vc)
        dv = -0.05 * v;
    else
        dv = 0.2 * v;
   if (vmin > 0)
       dv = 0.2 * v;
   if (vmax < 1e30)
       dv = -0.05 * v;
    return dv;
}

void MWater::vT_Ps_a(const double& P,const double& s, double& v, double& T){
    double v1,T1,p1,v2,T2,p2,px,dv,dsdT,u,h,s1,s2,sx,dT,dsdv,dpdT,dpdv,det;
    int count=1;
    // initialization
    T = T_P_s (P);      // saturation temperature at P
    v = v_PT_g(P,T);

    // iterate
    Prop_vT(v,T,px, u, h,sx);
    while ((s_not_converged(s,sx) || p_not_converged(P,px)) && (count <= maxiter)) {
        if (px < 0)
            dv=alter_v(v);
        else {
            dT = 0.001 * T;
            T1 = T + dT;
            v1 = v;
            Prop_vT(v1,T1,p1, u, h,s1);
            dv = 0.00001 * v;
            if (v <= vc)
                dv = -dv;
            v2 = v + dv;
            T2 = T;
            Prop_vT(v2,T2,p2, u, h,s2);

            dsdT = (s1 - sx) / dT;
            dsdv = (s2 - sx) / dv;
            dpdT = (p1 - px) / dT;
            dpdv = (p2 - px) / dv;
            det = dsdT * dpdv - dpdT * dsdv;
            dT = ((s - sx) * dpdv - (P - px) * dsdv) / det;
            dv = (dsdT * (P - px) - dpdT * (s - sx)) / det;
        }
        vregulate(v,dv);
        Tregulate(T,dT);
        T += dT;
        v += dv;
        count++;
        Prop_vT(v,T,px, u, h,sx);
    }
    if (count > maxiter)
        cerr << "vT_Ps_a no convergence" << endl;
}

double MWater::T_Ps_a(const double& P,const double& s){
    double v,T;
    vT_Ps_a(P, s, v, T);
    return T;
}

double MWater::v_Ps_a(const double& P,const double& s){
    double v,T;
    vT_Ps_a(P, s, v, T);
    return v;
}

void MWater::vT_Ph_a(const double& P,const double& h, double& v, double& T){
    double v1,T1,p1,v2,T2,p2,px,dv,u,h1,h2,hx,s,dT,dpdT,dpdv,det,dhdT,dhdv;
    int count=1;

    // initialization
    T = T_P_s (P);      // saturation temperature at P
    v = v_PT_g(P,T);

    // iterate
    Prop_vT(v,T,px, u, hx,s);
    while ((h_not_converged(h, hx, T) || p_not_converged(P, px)) && (count <= maxiter)){
        if (px < 0)
            dv=alter_v(v);
        else {
            dT = 0.001 * T;
            T1 = T + dT;
            v1 = v;
            Prop_vT(v1,T1,p1, u, h1,s);
            dv = 0.00001 * v;
            if (v <= vc)
                dv = -dv;
            v2 = v + dv;
            T2 = T;
            Prop_vT(v2,T2,p2, u, h2,s);
            dhdT = (h1 - hx) / dT;
            dhdv = (h2 - hx) / dv;
            dpdT = (p1 - px) / dT;
            dpdv = (p2 - px) / dv;
            det = dhdT * dpdv - dpdT * dhdv;
            dT = ((h - hx) * dpdv - (P - px) * dhdv) / det;
            dv = (dhdT * (P - px) - dpdT * (h - hx)) / det;
        }
        vregulate(v,dv);
        Tregulate(T,dT);
        T += dT;
        v += dv;
        count++;
        Prop_vT(v,T,px, u, hx,s);
    }
    if (count > maxiter)
        cerr << "vT_Ph_a no convergence" << endl;
}

double MWater::v_Ph_a(const double& P,const double& h){
    double v,T;
    vT_Ph_a(P, h, v, T);
    return v;
}

double MWater::T_Ph_a(const double& P,const double& h){
    double v,T;
    vT_Ph_a(P, h, v, T);
    return T;
}

double MWater::v_PT_a(const double& P,const double& T){
//TODO why is it possible to get a unique v from P,T ?
    double v2,p2,v,dv,u,h,dpdv,s;
    double dvm, dTm, dva, dTa, vt;
    double px = 0;
    int count=1;
    double dT = 0.0;    //never changed ?!
    int kbr=0;
    double dvbf=1.0;
    double vmin = 0;
    double vmax = 1E30;
    double pmin = 0;
    double pmax = 1E30;
    const double dvs1=2 * vc, dvs2=0.7 * vc;
    bool flag=false;
    // initialization
    v = v_PT_g(P,T);
    while (p_not_converged(P,px) && (count <= maxiter)) {
        Prop_vT(v,T,px, u, h, s);
        if (px >= 0) {
            dv = 0.00001 * v;
            if (v <= vc)
                dv = -dv;
            v2 = v + dv;
            Prop_vT(v2,T,p2, u, h, s);
            dpdv = (p2 - px) / dv;
            if (dpdv > 0)
                flag = true;
            else {
                if ((px > P) && (v > vmin)) {
                    vmin = v;
                    pmin = px;
                }
                if ((px < P) && (v < vmax)) {
                    vmax = v;
                    pmax = px;
                }
                if (vmin >= vmax)
                    cerr << "ERROR!";
                if ((vmin > 0) && (vmax < 1E30))
                    kbr = 1;
                if (dpdv == 0) {
                    dvbf = 0.5;
                    flag = true;
                } else {
                    dvbf = 1;
                    dv = (P - px) / dpdv;
                    dT = 0;
                }
            }
        }
        if ((px < 0) || flag) {     // flag==true if dpdv>=0
            if (kbr == 0)
                dv=alter_v(v, vmin, vmax);
            else {
                dpdv = (pmax - pmin) / (vmax - vmin);
                v = vmax;
                px = pmax;
                dv = dvbf * (P - px) / dpdv;
                dT = 0;
                dvbf = 0.5 * dvbf;
            }
        }
        dvm = 0.2 * v;
        if ((v<dvs1) || (v<dvs2))       // ?? sense ??
            dvm = 0.5 * dvm;
        dTm = 0.1 * T;
        if (kbr != 0) {
            vt = v + dv;
            if (!((vt >= vmin) && (vt <= vmax)))
                dv = vmin + (P - pmin) * (vmax - vmin) / (pmax - pmin) - v;
        }
        dva = abs(dv);
        dTa = abs(dT);
        if (dva > dvm)
            dv = dv * dvm / dva;
        if (dTa > dTm)
            dT = dT * dTm / dTa;
        //T += dT;          always dT = 0
        v += dv;
        count++;
    }
    return v;
}

void MWater::vregulate(const double& v,double& dv){
    double dva = abs(dv);
    double dvm = 0.2 * v;
    const double dvs1=2 * vc, dvs2=0.7 * vc;
    if ((v<dvs1) || (v<dvs2))       // ?? sense ??
        dvm = 0.5 * dvm;
    if (dva >dvm)
        dv = dv * dvm / dva;
}

void MWater::Tregulate(const double& T,double& dT){
    double dTm = 0.1 * T;
    double dTa = abs(dT);
    if (dTa > dTm)
        dT = dT * dTm / dTa;

}

double MWater::x_fraction(const double &var, const double &liq, const double &gaz){
    if (gaz < var)
        return -1.0;
    else if (liq > var)
        return -1.0;
    else
        return (var - liq) / (gaz - liq);
}

void MWater::Prop_vT_a(const double &v, const double &T, double &P, double &u, double &h, double &s, double &x) {
    double Tsat, vg, vl, hg, hl, sg, sl, ul, ug;
    double dpv;
    P = P_T_s(T);
    Prop_P_s(P, vl, vg, Tsat, ul, ug, hl, hg, sl, sg);
    if (v > vg){
        PWater_vT(v, T, P, u, h, s, dpv);
        x = -1.0;    //x = 1.0;
    }
    else if (v < vl){
        PWater_vT(v, T, P, u, h, s, dpv);
        x = -1.0;    //x = 0.0;
    }
    else {
        x = (v - vl) / (vg - vl);
        s = (1 - x) * sl + x * sg;
        h = (1 - x) * hl + x * hg;
        u = (1 - x) * ul + x * ug;
    }
}

double MWater::x_fraction(const double & T,const double& P, const double& var, char choice) {
    double DPDT, Psat, vg, Df, vf, x,uf,ug,hg,hf,sg,sf;
    SWater(T, Psat, DPDT);
    if (P < Psat)
        x = -1.0;   //x = 1.0;
    else if (P > Psat)
        x = -1.0;   //x = 0.0;
    else {
        DWater(T, Df);
        vf = 1.0 / Df;

        switch (choice) {
            case 'v': //Saturation pressure p and dpdT at T
                vg = v_PT_g( Psat,T);
                x = x_fraction(var, vf, vg);
                break;
            case 'u':
                Prop_PT_g(Psat,T,vg,ug,hg,sg);
                uf = ug + (vg - vf) * (P - T * DPDT);
                x = x_fraction(var, uf, ug);
                break;
            case 'h':
                Prop_PT_g(Psat,T,vg, ug,hg,sg);
                hf = hg - T * (vg - vf) * DPDT;
                x = x_fraction(var, hf, hg);
                break;
            case 's':
                Prop_PT_g(Psat,T,vg, ug,hg,sg);
                sf = sg - (vg - vf) * DPDT;
                x = x_fraction(var, sf, sg);
                break;
            default:
                cerr << "Invalid choice in x_fraction: " << choice << endl;
                x = -1.0;
        }
    }
    return x;
}

/*
MWater::MWater(const char* iFluidName) {

    const char sep = '\t';
    const int nbTitleLines = 2; //nombre de lignes avec des titres
    const int nbTitleColumns = 1; //nombre de colonnes avec des titres ie sans chiffre ie
    const char* file_in = "dataBase/DataFluids.txt";
    // Debut du traitement
    //cout << "Test Input Ouput file : " << file_in << endl << endl;

    TInputFile  InputFile((char*)file_in, nbTitleLines, nbTitleColumns, sep);
    InputFile.open();



    // Tests
 //   TOutputFile OutputFile((char*)"iofdata.out");
 //   OutputFile.open();
    int nbRecords = InputFile.nRecords();
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    Records = *InputFile.GetRecords();

    FluidName = (char*)iFluidName;
    TSequence* PtrSeq = Get(Records, FluidName);

    vc = PtrSeq->Get(0);
    Pc = PtrSeq->Get(1);
    Tc = PtrSeq->Get(2);
    Tnb = PtrSeq->Get(3);
    Mw = PtrSeq->Get(4);
    w = PtrSeq->Get(5);
    DipM = PtrSeq->Get(6);

    // cpo coefficients
    c0 = PtrSeq->Get(7);
    c1 = PtrSeq->Get(8);
    c2 = PtrSeq->Get(9);
    c3 = PtrSeq->Get(10);
    c4 = PtrSeq->Get(11); //pas de 5eme coefficient cp dans le fichier d'entrée

    // reference state default initialisation
    // liquid line
    Tref = 273.15;	// [K]
    href = 200.0E3;	// [J/kg]
    sref = 1.0E3;	// [J/(kg K)]
    InputFile.close();
    init();
} */

MWater::MWater(const char* iFluidName,const char* file_in ) {
    const char sep = '\t';
    const int nbTitleLines = 2; //nombre de lignes avec des titres
    const int nbTitleColumns = 1; //nombre de colonnes avec des titres ie sans chiffre ie

    TInputFile InputFile((char *) file_in, nbTitleLines, nbTitleColumns, sep);
    InputFile.open();
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    Records = *InputFile.GetRecords();

    FluidName = (char *) iFluidName;
    TSequence *PtrSeq = Get(Records, FluidName);

    vc = PtrSeq->Get(0);
    Pc = PtrSeq->Get(1);
    Tc = PtrSeq->Get(2);
    Tnb = PtrSeq->Get(3);
    Mw = PtrSeq->Get(4);
    w = PtrSeq->Get(5);
    DipM = PtrSeq->Get(6);

    // cpo coefficients
    c0 = PtrSeq->Get(7);
    c1 = PtrSeq->Get(8);
    c2 = PtrSeq->Get(9);
    c3 = PtrSeq->Get(10);
    c4 = PtrSeq->Get(11); //pas de 5eme coefficient cp dans le fichier d'entrée

    // reference state default initialisation
    // liquid line
    Tref = 273.15;    // [K]
    href = 200.0E3;    // [J/kg]
    sref = 1.0E3;    // [J/(kg K)]
    InputFile.close();
    MWater::init();
}


double MWater::DPDV(const double&V,const double&T)
{
    double P,u,h,s,DPV;
    PWater_vT(V, T, P, u, h, s, DPV);
/*    double D, Q, QR, QT, X, Y, Z, TA, T1, B, A1, A2, A3, A4, A5, A6, A7, A, R1, R2, R3, R4, R5, R6, R7,DPV;

    D = 1 / V;
    TA = 1000 / T;            // =Ta/T
    X = D - 1000;            // X=D-D(AJ) where D(AJ) =1000 for J >1
    Y = D - 634;            // Y=D-D(A1) where D(A1)=634
    Z = exp(-4.8e-3 * D);    // Z=EXP(-E*D)
    T1 = TA - 1000/Tc;        // T1=TA-TA(AJ) where TA(AJ)=1.544912=Ta/Tc=1000/647.286 for J=1
    B = TA - 2.5;            // B=TA-TA(AJ) where TA(AJ)=2.5 for  J >1

    A1 = 2.9492937e-2 + (-1.3213917e-4 + (2.7464632e-7 + (-3.6093828e-10 + (3.4218431e-13 + (-2.4450042e-16 +
                                                                                             (1.5518535e-19 +
                                                                                              5.9728487e-24 * Y) * Y) *
                                                                                            Y) * Y) * Y) * Y) * Y +
         Z * (-4.1030848e-1 - 4.1605860e-4 * D);
    A2 = -5.1985860e-3 + (7.7779182e-6 + (-3.3301902e-8 + (-1.6254622e-11 + (-1.7731074e-13 + (1.2748742e-16 +
                                                                                               (1.3746153e-19 +
                                                                                                1.5597836e-22 * X) *
                                                                                               X) * X) * X) * X) * X) *
                         X + Z * (3.3731180e-1 - 2.0988866e-4 * D);
    A3 = 6.8335354e-3 + (-2.6149751e-5 + (6.5326396e-8 - 2.6181978e-11 * X) * X) * X +
         Z * (-1.3746618e-1 - 7.3396848e-4 * D);
    A4 = -1.5641040e-4 + (-7.2546108e-7 + (-9.2734289e-9 + 4.3125840e-12 * X) * X) * X +
         Z * (6.7874983e-3 + 1.0401717e-5 * D);
    A5 = -6.3972405e-3 + (2.6409282e-5 + (-4.7740374e-8 + 5.6323130e-11 * X) * X) * X +
         Z * (1.3687317e-1 + 6.4581880e-4 * D);
    A6 = -3.9661401e-3 + (1.5453061e-5 + (-2.9142470e-8 + 2.9568796e-11 * X) * X) * X +
         Z * (7.9847970e-2 + 3.9917570e-4 * D);
    A7 = -6.9048554e-4 + (2.7407416e-6 + (-5.1028070e-9 + 3.9636085e-12 * X) * X) * X +
         Z * (1.3041253e-2 + 7.1531353e-5 * D);
    A = A2 + B * (A3 + B * (A4 + B * (A5 + B * (A6 + B * A7))));
    Q = A1 + T1 * A;
    //QT = -TA / T * (A + T1 * (A3 + B * (2 * A4 + B * (3 * A5 + B * (4 * A6 + B * 5 * A7)))));

    R1 = -1.3213917e-4 + (2.7464632e-7 * 2 + (-3.6093828e-10 * 3 + (3.4218431e-13 * 4 + (-2.4450042e-16 * 5 +
                                                                                         (1.5518535e-19 * 6 +
                                                                                          5.9728487e-24 * 7 * Y) * Y) *
                                                                                        Y) * Y) * Y) * Y +
         Z * (-4.1605860e-4 - 4.8e-3 * (-4.1030848e-1 - 4.1605860e-4 * D));
    R2 = 7.7779182e-6 + (-3.3301902e-8 * 2 + (-1.6254622e-11 * 3 + (-1.7731074e-13 * 4 + (1.2748742e-16 * 5 +
                                                                                          (1.3746153e-19 * 6 +
                                                                                           1.5597836e-22 * 7 * X) * X) *
                                                                                         X) * X) * X) * X +
         Z * (-2.0988866e-4 - 4.8e-3 * (3.3731180e-1 - 2.0988866e-4 * D));
    R3 = -2.6149751e-5 + (6.5326396e-8 * 2 - 2.6181978e-11 * 3 * X) * X +
         Z * (-7.3396848e-4 - 4.8e-3 * (-1.3746618e-1 - 7.3396848e-4 * D));
    R4 = -7.2546108e-7 + (-9.2734289e-9 * 2 + 4.3125840e-12 * 3 * X) * X +
         Z * (1.0401717e-5 - 4.8e-3 * (6.7874983e-3 + 1.0401717e-5 * D));
    R5 = 2.6409282e-5 + (-4.7740374e-8 * 2 + 5.6323130e-11 * 3 * X) * X +
         Z * (6.4581880e-4 - 4.8e-3 * (1.3687317e-1 + 6.4581880e-4 * D));
    R6 = 1.5453061e-5 + (-2.9142470e-8 * 2 + 2.9568796e-11 * 3 * X) * X +
         Z * (3.9917570e-4 - 4.8e-3 * (7.9847970e-2 + 3.9917570e-4 * D));
    R7 = 2.7407416e-6 + (-5.1028070e-9 * 2 + 3.9636085e-12 * 3 * X) * X +
         Z * (7.1531353e-5 - 4.8e-3 * (1.3041253e-2 + 7.1531353e-5 * D));
    QR = R1 + T1 * (R2 + B * (R3 + B * (R4 + B * (R5 + B * (R6 + B * R7)))));


   //   DPV = -461.51 * T * D * D * (1 + 2 * Q * D + 2 * QR * D * D); // faux ( a refaire)
    DPV = R * T *(1 + 2 * Q * vc/V + 2 * QR * vc/V * vc/V); // faux ( a refaire)        */
    return DPV;
}




/*
 *
 *


void Water::Sonic_PT_g(double &Psup, double &Tsup, double &vLaval, double &PLaval, double &TLaval, double &Csonic) {
    double hsup, Ssup,usup, hl, sl,ul;
    double DPV, DPT, CV, CP, Ci;
    //LK_Prop_PT_g(Psup, Tsup, vLaval, hsup, Ssup, DPV, DPT, CV, CP, Csonic);
    Prop_PT_g (Psup, Tsup, vLaval,usup, hsup, Ssup, DPV, DPT, CV, CP, Csonic);
    do{
    Ci = Csonic;
    hl = hsup - sqr(Csonic) / 2;
    Prop_hs_g(hl, Ssup, vLaval, ul,PLaval, TLaval);
    Prop_PT_g(PLaval, TLaval, vLaval,ul, hl, sl, DPV, DPT, CV, CP, Csonic);
    //if IterFailed (i,50, W_Sonic_PT_gFailed) then	exit(W_Sonic_PT_g);
    } while (abs(Csonic - Ci) / Csonic < 0.001);

}

 void Water::Prop_hs_g(double &h, double &s, double &v,double u, double &P, double &T) {
    double const Pc= 22.089e6;
    double Pi, m, Ti, si, hi;
    double cp;
    P = Pc / 10;
    PhTGazW (P, h, T);
    do {
        Pi = P;
        // iteration pour trouver s=f(P,h)
        do {
            Ti = T;
            Prop_PT_g(Pi, T,v,u, hi, si);
            //if IterFailed(i, 100, W_Prop_hs_gFailed)
            //then exit(W_Prop_hs_g);
            cp = 2503.2 - 2.5183 * T + 0.0040416 * sqr(T) - 2.1642e-6 * XpwrY(T, 3) + 0.4092e-9 * XpwrY(T, 4);
            T = Ti - (hi - h) / cp;
        } while((abs(hi - h) < 1E-5) && (abs(Ti - T) < 1E-3));
        //fin de l'iteration interm�diaire
        m = -v / T;
        P = Pi - (si - s) / m;
    } while ((abs(Pi - P) < 1E-3) && (abs(si - s) < 1E-3));
    Prop_PT_g (P, T,v,u,h,s);
}
*/

/*
void Water::W_Prop_vs_bg(double &v, double &s, double &P, double &T, double &h, double &x) {
    double vi, Ti;
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq, DPTs;

    int i =0;
    T=T_v_s(v);
    SWater(T, P, DPTs);
    Prop_P_s (P, T, hLiq, hGaz, sLiq, sGaz);
    vLiq = v_PT_l (T, P);
    vGaz = v_PT_g (T, P);

    if (sGaz > s)
    {do {
            Ti = T;
            x = (s - sLiq) / (sGaz - sLiq);
            vi = vLiq + x * (vGaz - vLiq);
            T = Ti * XpwrY(vi / v, 0.05);
            SWater(T, P, DPTs);
            Prop_P_s (P, T, hLiq, hGaz, sLiq, sGaz);
            vLiq = v_PT_l (T, P);
            vGaz = v_PT_g (T, P);
            //if (IterFailed (i,100, W_Prop_vs_bgFailed))break;
        } while ((abs(Ti - T) < 0.001) && (abs(vi - v) / v < 0.00001));
        x = (s - sLiq) / (sGaz - sLiq);
        h = hLiq + x * (hGaz - hLiq);}
    else {
        Prop_vs_g(v, s, P, T, h);
        x = 1;}

}

void Water::W_Prop_vP_bg(double &v, double &P, double &T, double &h, double &s, double &x) {
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq;
    Prop_P_s (P, T, hLiq, hGaz, sLiq, sGaz);
    vLiq = v_PT_l (T, P);
    vGaz = v_PT_g (T, P);
    if (v > vGaz)
    {PvThsGazW (P, v, T, h, s);
        x =1.0;}

    else //biohasic zone
    {x = (v - vLiq) / (vGaz - vLiq);
        s = (1 - x) * sLiq + x * sGaz;
        h = (1 - x) * hLiq + x * hGaz;}

}
 void Water::W_GTH_PT_l(double &P, double &T, double &Cv, double &Cp, double &Sv, double &Av, double &Bp, double &Gt) {
    double const R = 8314.3;
    double const W = 18.016;
    double Px, Tx, v, h, s, u;
    double DPT, DPV;

    Px= P;
    Tx= T;
    Prop_PT_g (Px, Tx, h, s, v);
    PWater (Tx, Px, v, u, h, s, DPV);
    SWater (Tx, Px, DPT);
    Cp = 31711 - 286.8 * T + 1.105 * sqr(T) - 1.874e-3 * XpwrY(T, 3) + 1.188e-6 * XpwrY(T, 4);
    //CV := CP - R / W;
    Cv = Cp;
    Sv = sqrt(-1000 * Cp / Cv * DPV) * v;
    Av = P / T / DPT;
    Gt = -P / v / DPV;
    Bp = Av / Gt;
}

 void Water::W_GTH_PT_g(double &P, double &T, double &Cv, double &Cp, double &Sv, double &Av, double &Bp, double &Gt) {
    double const  R = 8314.3;
    double const W = 18.016;
    double Px, Tx, v, h, s, u;
    double DPT, DPV;

    Px= P;
    Tx= T;
    Prop_PT_g (Px, Tx, h, s, v);
    PWater (Tx, Px, v, u, h, s, DPV);
    SWater (Tx, Px, DPT);
    Cp = 2503.2 - 2.5183 * T + 0.0040416 * sqr(T) - 2.1642e-6 * XpwrY(T, 3) + 0.4092e-9 * XpwrY(T, 4);
    Cv = Cp - R / W;
    Sv = sqrt(-1000 * Cp / Cv * DPV) * v;
    Av = P / T / DPT;
    Gt = -P / v / DPV;
    Bp = Av / Gt;

}

 void Water::W_Prop_vh_bg(double &v, double &h, double &P, double &T, double &s, double &x) {
    double vi, Ti;
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq, DPTs;
    int i(0);
    T=T_v_s(v);
    SWater(T, P, DPTs);
    Prop_P_s (P, T, hLiq, hGaz, sLiq, sGaz);
    vLiq = v_PT_l (T, P);
    vGaz = v_PT_g (T, P);

    //LK_Prop_P_s(P, vGaz, vLiq, T, hGaz, hLiq, sGaz, sLiq, DPTs, DPVs);
    if (hGaz > h)
    { Ti = T;
    x = (h - hLiq) / (hGaz - hLiq);
    vi = vLiq + x * (vGaz - vLiq);
    T = Ti * XpwrY(vi / v, 0.05);
    SWater(T, P, DPTs);
    Prop_P_s (P, T, hLiq, hGaz, sLiq, sGaz);
    vLiq = v_PT_l (T, P);
    vGaz = v_PT_g (T, P);}

   // if (IterFailed (i,50, W_Prop_vh_bgFailed))
    //{exit(W_Prop_vh_bg);}

    do {s = sLiq + x * (sGaz - sLiq);
    x = (h - hLiq) / (hGaz - hLiq);} while ((abs(Ti - T) < 0.001) && (abs(vi - v) / v < 0.0001));

    //else
    //{W_Prop_vh_g(v, h, P, T, s);
    //x := 1;}

}

 void Water::W_Prop_Px_b(double &P, double &x, double &v, double &Tsat, double &h, double &s) {
    double vGaz, vLiq, hGaz, hLiq, sGaz, sLiq;
    Prop_P_s (P, Tsat, hLiq, hGaz, sLiq, sGaz);
    vLiq = v_PT_l (Tsat, P);
    vGaz = v_PT_g (Tsat, P);
    v = (1 - x) * vLiq + x * vGaz;
    h = (1 - x) * hLiq + x * hGaz;
    s = (1 - x) * sLiq + x * sGaz;
}

void Water::W_Prop_Ps_bg(double &P, double &s, double &v, double &T, double &h, double &x) {
    double hGaz, hLiq, sGaz, sLiq;
    Prop_P_s (P, T, hLiq, hGaz, sLiq, sGaz);
    if (sGaz < s)
    {PsThGazW (P, s, T, h);
    v = v_PT_g (T, P);
    x = 1.0;}
    else
    {x = (s - sLiq) / (sGaz - sLiq);
    v = (1 - x) * v_PT_l (T, P) + x * v_PT_g (T, P);
    h = (1 - x) * hLiq + x * hGaz;}

}

void Water::W_Prop_Ph_bg(double &P, double &h, double &v, double &T, double &s, double &x) {
    double hGaz, hLiq, sGaz, sLiq,u;
    Prop_P_s (P, T, hLiq, hGaz, sLiq, sGaz);
    if (hGaz < h)
    {PhTGazW (P, h, T);
    Prop_PT_g (P, T,v,u,h,s);
    x = 1.0;}

    else
    {x = (h - hLiq) / (hGaz - hLiq);
    v = (1 - x) * v_PT_l (T, P) + x * v_PT_g (T, P);
    s = (1 - x) * sLiq + x * sGaz;};
}

void Water::W_Prop_PTx_bg(double &P, double &T, double &v, double &h, double &s, double &x) {

    double Tsat, hLiq, hGaz, sLiq, sGaz;
    //input:(P,T) ou (P,x)
    Prop_P_s (P, Tsat, hLiq, hGaz, sLiq, sGaz);
    if (T > Tsat)
    {Prop_PT_g (P, T, h, s);
    v = v_PT_g (T, P);
    x = 1.0;}

    else
    {W_Prop_Px_b(P, x, v, T, h, s);}

}

 void Water::Prop_PT_l(double &P, double &T, double &h, double &s) {
    double psat, pes, vg, ug, hg, sg, dpdT, ves, ues, hes, ses;
    double Vp, Up, Hp, Sp, Vps, Ups, Hps, Sps, pccv, fact1, v, u;

    v = v_PT_l(T, P);
    pccv = 22.089e6;
    PropertyWater(T, P, v, u, h, s, 2);
    Vp = v;
    Up = u;
    Hp = h;
    Sp = s;
    Vps = v;
    SWater(T, psat, dpdT);
    PropertyWater(T, psat, Vps, Ups, Hps, Sps, 2);
    Prop_T_s ( T, pes,hes,hg,ses,sg);
    if (abs(pes - psat) > (1e-3 * psat))
            std::cout<<"error pes psat "<<pes<<" " <<psat<<std::endl;
    fact1 = (ln(pccv) - ln(P)) / (ln(pccv) - ln(psat));
    h = Hp + (hes - Hps) * fact1;
    s = Sp + (ses - Sps) * fact1;

}

 void Water::Prop_vh_g(double &v, double &h, double &P, double &T, double &s) {
    double vi, Pi, Ti, hi, m,ui;
    double dpv, dpt, cv, cp, sv;

    int i = 0;
    T=T_v_s(v);
    SWater(T, P, dpt);
    do {
        Pi = P;
        //iteration pour trouver hi=f(v,Pi)
        do {
            Ti = T;
            Prop_PT_g(P, Ti, vi,ui, hi, s, dpv, dpt, cv, cp, sv);
            //if IterFailed(i, 100, W_Prop_vh_gFailed)
            //then exit(W_Prop_vh_g);
            T = Ti + (vi - v) * dpv / dpt;
        } while ((abs(Ti - T) < 1E-4) && (abs(vi - v) < 1E-8));
        //fin de l'it�ration intermediaire
        m = cv / dpt + v;
        P = Pi - (hi - h) / m;
    } while( (abs(P - Pi) < 1E-4) & (abs(hi - h) < 1E-4));
        Prop_PT_g(P, T, vi,ui, hi, s, dpv, dpt, cv, cp, sv);

}

 void Water::vTPhsLiqW(double &v, double &T, double &P, double &h, double &s) {
    double u;
    Prop_vT_a(v,T,P, u, h,s);
}
 */


//Debut du constructeur-------------------------------------------------
/*
void Glycol::Prop_vT_a(const double&V,const double&T,double&P, double&U, double&H,double& S)
{


}

double Glycol::v_T_l(const double &T) {
    double D,v,a,b,c;
    a=-0.0012;

    b= 0.7619 * N * N - 1.1891 * N - 0.0835;
    c= 177.18 * N + 1000.4;

    D=a*T*T + b*T +c;
    return 1/D;
}

double Glycol::ThermCond(const double &T){
    double a,b,c;
    a=0.000000039563 *N*N*N-0.000000070562*N*N+0.000000042973*N -0.0000000088081;
    b=-0.000026472*N*N*N+0.000046179 *N*N-0.000028958*N+0.0000068746;
            c=0.0044667*N*N*N-0.0074192*N*N+0.0043243*N-0.00065917;
    return a*T*T + b*T +c;
}

double Glycol::Cp(const double &T){
    double a,b;
    a=-0.002*N*N+0.0066*N+0.0003;
            b=-0.7078*N*N-1.1999*N+4.1886;
    return a*T+b;
}

double Glycol::ThermCond(const double &T) {
    double a, b, c;

}

void Prop_PT_l (const double&P,const double& T,double& v, double &u, double &h, double &s);
void TransProp_PT_l ( double& P, double& T,double& ViscDyn,double& ViscCin,double& Rho,double& ThermCond,double& Prandtl);*/