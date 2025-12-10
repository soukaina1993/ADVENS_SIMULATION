//
// Created by lucile.schulthe on 14.06.2022.
//

#ifndef FLUIDS_LK_H
#define FLUIDS_LK_H


#include "../math/maths.h"
#include <string>
#include "../iofile/iofile.h"
#include "LeeK1.h"

struct LeeKesslerCst {
    double BT, GM, B1, B2, B3, B4, C1, C2, C3, C4, D1, D2, ZCP, ZA, ZB, P1, P2, P3, P4;
};

// Consts for Yen and Woods equation in subroutine ESTD
struct LKYWCstRec {
    double A1, A2, A3, A4, A5, B1, B2, B3, B4, B5;
    double C1, C2, C3, C4, D1, D2, D3, D4, D5;
};

struct LKConstRec {
    double PO, WR, RBAR, TDTM;
    LKYWCstRec YWCst3, YWCst5, YWCst9;

};


class LK1:public FluidCore {
public:
    void Init_LKConst();

    void DPT_VTor(double VRo,double  VRr,double TR,double  &DPTo,double &DPTr,double  &DPT);
    void DPV_VTor (double VRo, double VRr, double TR, double &DPVo, double &DPVr,double  &DPV);
    void DPTor_VTr (double & DPT, double VR,double TR, LeeKesslerCst &CSTS);

    void DPVor_VTr(double &DPV,double VR, double TR, LeeKesslerCst &CSTS);
    void Por_VTr(double &PR, double& VR,double &TR, LeeKesslerCst &CSTS);
 //      double Por_VTr(double TR, double VR, LeeKessler::LeeKesslerCst &CSTS); ??

    double TL, roL, Tfp, Roc;
    double ADA0, ADA1, HREFO, SREFO, ALPHA, BETA, GAMMA, DELTA,AF;

//    double EO, ER;
    double kPc;
    LKConstRec LKconst;
protected:
    LeeKesslerCst CstO, CstR;
};

class LK2: public LK1 {
public:
    void LK_Prop_PT_g(double P, double T,  double &V,  double &H,  double &S,  double &DPV, double &DPT,  double &CV,  double &CP, double &SV);

    void iLK_Prop_Ps_g(double P, double s, double &v, double &T, double &h);
    void iLK_Prop_Ph_g(double P, double h, double &v, double &T, double &s);
    void iLK_Prop_vP_g (double v,double P,double & T, double &h, double &s);
    void iLK_Prop_vT_g (double v,double T,double & P,double & h,double & s);
    void iLK_Prop_Th_g (double T, double h,double & v, double &P,double & s);
    void iLK_Prop_Ts_g (double T,double s,double & v, double &P, double &h);

    void iiLK_Prop_hs_g (double h,double s, double &v, double &P, double &T);
    void iiLK_Prop_vs_g(double v, double s, double &P, double &T, double &h);
    void iiLK_Prop_su_g (double s,double u, double &v,double & P,double & T,double & h);
    void iiLK_Prop_vh_g(double v, double h, double &P, double &T, double &s);

    void Vr_PTr_g(double &VR, double &PR,double & TR);
    void Pr_vTr_g(double &VR, double &PR,double & TR);

    void VRro_PTr_g (double PR,double TR, double &VR,double & VRo,double & VRr,double & Eo,double & Er);

    void H_VT(double VRo,double  VRr,double  PR,double  TR,double  Eo,double  Er,double  &H);
    void S_VT(double VRo,double  VRr,double TR,double Eo,double Er, double  &S);
    void CVP_VT (double VRo,double  VRr,double  Eo,double  Er,double  DPVO,double  DPVR,double  DPTO,double  DPTR, double  TR, double  &CV, double &CP);

    void LK_Prop_P_s (double P, double &vLiq, double &vGaz, double &Ts, double &hLiq, double &hGaz, double &sLiq, double &sGaz, double &DPTs, double &DPVs);
    void Ts_P(double P, double &Ts, double &DPTS);
    void Ts_v(double &Ts,double v);
    void ss_T(double &ss,double T);
    void Ps_T(double &P,double T);
    void hs_T(double &hs, double T);
    void Ts_h(double &Ts,double h);

    double DPTS_PT_sg (double Ps, double Ts);

    void LK_Sonic_PT_g(double pSup, double Tsup, double &Laval,double &PLaval, double &TLaval, double &Csonic);
    void LK_TransProp_PT_g (double P,double T, double &ViscDyn,double & ViscCin,double & Ro,double & ThermCond, double & Prandtl);



private:
    double PRs_T(double &Tr, LeeKesslerCst &CSTS);
    double Zs(double &Tr, LeeKesslerCst &CSTS);
    double root(double Pr, double Tr, double VRInf, double VRSup, double dVR, LeeKesslerCst &CSTS);
    void iVor_PTr_g2(double& VR, double PR, double TR, LeeKesslerCst &CSTS);
    void iVor_PTr_g(double &VR, double PR,double TR, LeeKesslerCst &CSTS );




};

class LK3: public LK2 {
public:
    void InitLKFluid();

    void LK_Prop_PT_l (double P,double T,double& V, double& RO,double&  H,double&  S,double&  DPV,double&  DPT);
    double LK_v_PT_l (double P,double T);
    void iLK_T_Ph_l (double P,double H,double& T);
    void iLK_Th_Ps_l (double& T,double& H, double P,double S);
    void LK_Cvp_PT_l (double P,double T, double& CPL, double& CVL, double& CPF);

    void LK_Prop_P_sl(double P, double &vf, double& hf, double& sf, double& hfg, double& sfg, double vg, double hg,
                      double sg);
    double LK_vs_T_sl(double T);

    void LK_Prop_Ps_bg(double P, double s, double &v, double &T, double &h, double &x);
    void iLK_Prop_Ph_bg(double P, double h, double &v, double &T, double &s, double &x);
    void iLK_Prop_vs_bg (double v,double s,double &P,double &T,double &h,double &x);
    void iLK_Prop_vh_bg(double v, double s, double &P, double &T, double &h, double &x);

    void LK_Prop_vP_bg (double v,double P,double& T,double&  h,double& s,double&  x);
    void LK_Prop_vT_bg (double v,double T,double&   P,double&  h, double& s,double&  x);
    void LK_TransProp_PT_l (double P,double T,double&  ViscDyn,double&  ViscCin,double&  Ro,double&  ThermCond,double&  Prandtl);
    void LK_TransProp_vT_l(double P,double T,double&  ViscDyn,double&  ViscCin,double&  Ro,double&  ThermCond,double&  Prandtl);

    void ROr_PTr_YW (double&  ROR,double PR,double TR,double ZCP);
    void VRro_PTr_l (double PR, double TR, double&  VR,double&  VRo,double&  VRr);
    void iVor_PTr_l (double & VR, double PR,double TR, LeeKesslerCst &CSTS);
    void ESTD (double & ROR, double PR, double TR, LKYWCstRec YWCst, double ZCP);

    void Pr_vTr_l(double &VR, double &PR,double & TR);
};


class LK4: public LK3 {
public:
//%%%%%%%%%%%% Initialisation %%%%%%%%%%%%
    void InitFluid(){InitLKFluid();};
    void Readfluid(const char* iFluidName);
    void InitLeeKesler(){Init_LKConst();};


//%%%%%%%%%%%% General %%%%%%%%%%%%

    void CriticProp (double & Vcrit,double & Pcrit,double & Tcrit,double & MM,double &Zcrit);
    void DomainProp (double &  Vmin, double & Vmax,double &  Pmin, double & Pmax, double & Tmin, double & Tmax);
    int Zone_PT(double P,double T); //return la zone

    // %%%%%%%%%%%% Zone Gazeuse %%%%%%%%%%%%}

    //------- input v,P -------}

    void Prop_vP_g (const double & v,const double &  P,double& T,double& h,double& s);
    double T_vP_g (const double & v,const double &  P);
    double h_vP_g (const double & v,const double &  P);
    double s_vP_g (const double & v,const double &  P);
    double u_vP_g (const double & v,const double &  P);

  //------- input v,T -------}

    void Prop_vT_g (const double & v,const double &  T,double& P,double&h, double&s);
    double P_vT_g (const double &v,const double &  T);
    double h_vT_g (const double & v, const double & T);
    double s_vT_g (const double & v,const double &  T);
    double u_vT_g (const double & v,const double &  T);

    //------- input v,h -------}

    void Prop_vh_g (const double & v,const double &  h,double& P,double& T,double& s);
    double P_vh_g (const double & v,const double &  h);
    double T_vh_g (const double & v,const double &  h);
    double s_vh_g (const double & v,const double &  h);
    double u_vh_g (const double & v,const double &  h);

    //------- input v,s -------}

    void Prop_vs_g (const double & v,const double &  s,double& P,double& T,double& h);
    double P_vs_g (const double & v, const double & s);
    double T_vs_g (const double & v,const double &  s);
    double h_vs_g (const double & v,const double &  s);
    double u_vs_g (const double & v,const double &  s);


    //------- input P,T -------

    void Prop_PT_g ( const double & P, const double & T,double&  v, double& h, double& s);
    void PropE_PT_g (const double & P, const double & T,double& v, double& h, double& s,double&  Cv,double&  Cp);
    double v_PT_g (const double & P,const double & T);
    double h_PT_g (const double & P,const double & T);
    double s_PT_g (const double & P,const double & T);
    double u_PT_g (const double & P,const double & T);

    void GTH_PT_g (double P,double T,double& Cv,double& Cp,double& Sv,double&Av, double&Bp,double& Gt);
    void Sonic_PT_g (double Psup, double Tsup,double& vLaval,double&PLaval, double& TLaval,double& Csonic);
    double Cv_PT_g (const double & P, const double & T);
    double Cp_PT_g (const double & P, const double & T);
    double Z_PT_g (const double & P,const double & T);
    double Zr_PT_g (const double & P,const double & T);
    void TransProp_PT_g (const double & P,const double &  T,double& ViscDyn, double& ViscCin, double& Ro, double&ThermCond,double& Prandtl);
    double SurfacTension(const double & P,const double &  T);
    //------- input P,h -------}

    void Prop_Ph_g (const double & P,const double &  h,double& v,double& T,double& s);
    double v_Ph_g (const double & P,const double &  h);
    double T_Ph_g (const double & P,const double &  h);
    double s_Ph_g (const double & P,const double &  h);
    double u_Ph_g (const double & P,const double &  h);

    //------- input P,s -------}

    void Prop_Ps_g (const double & P,const double &  s,double& v,double& T,double& h);
    double v_Ps_g (const double & P,const double &  s);
    double T_Ps_g (const double & P,const double &  s);
    double h_Ps_g (const double & P,const double &  s);
    double u_Ps_g (const double & P, const double & s);

    //------- input T,h -------}

    void Prop_Th_g (const double & T,const double &  h,double& v,double& P,double& s);
    double v_Th_g (const double & T,const double &  h);
    double P_Th_g (const double & T,const double &  h);
    double s_Th_g (const double & T,const double &  h);
    double u_Th_g (const double & T,const double &  h);

    //------- input T,s -------}

    void Prop_Ts_g (const double & T, const double & s,double& v,double& P,double& h);
    double v_Ts_g (const double & T, const double & s);
    double P_Ts_g (const double & T,const double &  s);
    double h_Ts_g (const double & T,const double &  s);
    double u_Ts_g (const double & T,const double &  s);

    //------- input T,u -------}

    //------- input h,s -------}
    void Prop_hs_g (const double & h,const double &  s,double& v,double& P,double& T);
    double P_hs_g (const double & h,const double &  s);
    //------- input s,u -------}
    void Prop_su_g (const double & s,const double &  u,double& v, double&P,double& T,double& h);
    double v_su_g (const double & s,const double &  u);
    double P_su_g (const double & s,const double &  u);
    double T_su_g (const double & s,const double &  u);
    double h_su_g (const double & s, const double & u);

    //%%%%%%%%%%%% Zone Saturee %%%%%%%%%%%%}

    void Prop_P_s (const double &Psat, double& vl,double& vg,double &T, double &hl, double &hg, double &sl, double &sg);
    double v_T_sl (double T);
    double Psat_T (double T);

    double hsat_T (double T);
    double Tsat_h (double h);
    double Tsat_v (double v);

    double Tsat_P (double P);
    double dPsat_PT_sg (double Ps,double Ts);

    //%%%%%%%%%%%% Zone Biphase %%%%%%%%%%%%}
    void Prop_Ps_bg(double P, double s, double &v, double &T, double &h, double &x);
    void Prop_Ph_bg(double P, double h, double &v, double &T, double &s, double &x);
    void Prop_vs_bg (double v,double s,double &P,double &T,double &h,double &x);
    void Prop_vh_bg(double v, double h, double &P, double &T, double &s, double &x);
    void Prop_Px_b (double P,double x, double & v,double & Tsat,double & h,double & s);
    void Prop_PTx_bg (double P, double &T,double & v,double & h,double & s,double & x);

    void Prop_vP_bg (double v,double P, double &T,double & h,double & s,double &x);
    void Prop_vT_bg (double v,double T, double& P,double & h, double &s,double &x);


    double x_Ph (double P,double h);

    //%%%%%%%%%%%% Zone Liquide %%%%%%%%%%%%}
    void GTH_PT_l (double P,double T,double &  Cv,double &  Cp,double &   Sv,double &   Av,double &   Bp, double &  Gt);
    void Prop_PT_l (const double & P,const double & T, double &v, double & h, double & s);
    void PropE_PT_l (double P,double T, double &v, double & h, double & s,double &  Cv,double &  Cp,double &  CpSat,double &  DPV,double &  DPT);

    double v_PT_l (const double & P,const double & T);
    double h_PT_l (const double & P,const double & T);
    double T_Ph_l (const double & P,const double & h);
    double v_Ph_l (const double & P,const double & h);
    void Th_Ps_l (const double & P,const double & s,  double & T, double &h);
    void TransProp_PT_l (const double & P,const double & T,  double & ViscDyn,  double &ViscCin,  double &Ro,  double &ThermCond, double & Prandtl);
};

#endif //FLUIDS_LK_H
