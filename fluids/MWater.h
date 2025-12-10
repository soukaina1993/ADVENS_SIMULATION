//
// Created by Yolaine Adihou on 17/04/2021.
//
#pragma once
#ifndef FLUIDSMAIN_CPP_WATER_H
#define FLUIDSMAIN_CPP_WATER_H
#include "../math/maths.h"
#include "FluidCore.h"
#include <string>
#include "../iofile/iofile.h"

using namespace std;

// const double TOK = 273.15;
// const double sqrt2=sqrt(2.0);
const double erh = 0.00001;
const double erp = 0.00001;
const double ers = 0.00001;


class MWater:public FluidCore {
public:
    explicit MWater(const char *iFluidName, const char* file_in="dataBase/DataFluids.txt");
    //MWater(const char *iFluidName);
    MWater();

    void init() override;

    void    CvCp_T_o(const double &T,double &cvo,double &cpo) override;
    double	h_T_o(	const double &T) override;
    double	s_T_o(	const double &T) override;
    void	IdealState(const double &T,double &cpo,double &cvo,double &ho,double &so) override;
    void	CvCp_vT(const double &v,const double &T,double &cv,double &cp) override;

    // EOS partial diff.
    void	dPvdPT_vT(const double &v,const double &T,double &dPdv,double &dPdT) override;

    // thermal factors
    void	AvBpGt_vT(const double &v,const double &T,double &Av,double &Bp,double &Gt) override;


    double	P_vT(const double &v,const double &T) override;
    double	T_vP(const double &v,const double &P) override;
    void	v_PT(const double &P,const double &T,double &vl,double &vg) override;
    double	u_vT(const double &v,const double &T) override;
    double	s_vT(const double &v,const double &T) override;
    double	h_vT(const double &v,const double &T) override;
    double	h_vPT(const double &v,const double &P,const double &T) override;

    void	Prop_vT(const double &v,const double &T,double &P,double &u,double &h,double &s) override;

    // SATURATION LINE FUNCTIONS
    //double W_Tsat_v ( double& v) { return 0; }

    // Saturation
    void Prop_P_s(const double &P, double& vl,double& vg,double &T,double& ul,double& ug, double &hl, double &hg, double &sl, double &sg) override;
    void Prop_T_s(const double &T, double& vl,double& vg,double &P,double& ul,double& ug, double &hl, double &hg, double &sl, double &sg) override;
    double T_v_s (const double& v);
    double P_T_s(const double &T) override;
    double T_P_s(const double &P) override;
    void P_T_s(const double& T, double& P, double& DPDT);
    void T_P_s(const double& P, double& T, double& DPDT);

    // Liquid
    void Prop_PT_l (const double&P,const double& T,double& v, double &u, double &h, double &s) override;
    void TransProp_PT_l (const double& P,const double& T,double& ViscDyn,double& ViscCin,double& Rho,double& ThermCond,double& Prandtl) override;
    void TransProp_vT_l (const double& P,const double& T,double& ViscDyn,double& ViscCin,double& Rho,double& ThermCond,double& Prandtl) override;
    void FlowProp_vT_l(const double &v,const double &T,double &cp,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) override;
    double v_PT_l(const double &P,const double &T) override;
    double v_T_l(const double &T);
    void CvCp_vT_l(const double &T,double &cv,double &cp);
    void CvCp_vT_l(const double &v,const double &T,double &cv,double &cp) override;


    // Gaz
    void Prop_PT_g (const double&P,const double& T,double& v,double& u,double& h,double& s) override;
    void Prop_PT_g (const double& P,const double& T,double& h,double& s);
    //void Prop_vT_g (const double& v, const double& T,double& P,double& u,double& h,double& s);
    double v_PT_g(const double &P,const double &T) override;
    void Prop_PT_g ( double& P, double& T,double& v,double& u,double& h,double& s,double& DPV,double& DPT,double& CV,double& CP,double& SV);
    void Prop_vs_g (const double&v, const double& s,double&P,double& T,double& u,double& h);
    void TransProp_PT_g (const double& P,const double& T,double& ViscDyn,double& ViscCin,double& Rho,double& ThermCond,double& Prandtl) override;
    void TransProp_vT_g (const double& v,const double& T,double& ViscDyn,double& ViscCin,double& Rho,double& ThermCond,double& Prandtl) override;
    void FlowProp_vT_g(const double &v,const double &T,double &cp,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) override;
    void CvCp_vT_g(const double &v,const double &T,double &cv,double &cp) override;
    void CvCp_vT_g(const double &T,double &cv,double &cp);
    //    void Prop_hs_g ( double&h, double& s,double&v,double u,double& P,double& T);

// all

    //void Prop_vT_a(const double&V,const double&T,double&P, double&U, double&H,double& S);
    void PWater_vT(const double&V, const double&T, double&P, double&U, double&H, double& S, double& DPV);
    void Prop_vT_a(const double&V,const double&T,double&P, double&U, double&H,double& S,double &x) override;
    void vT_Ph_a(const double &P,const double &h, double& v, double& T);
    double v_Ph_a(const double &P,const double &h) override;
    double T_Ph_a(const double &P,const double &h);
    void vT_Ps_a(const double &P,const double &s, double& v, double& T);
    double T_Ps_a(const double &P,const double &s) override;
    double T_vP_a (const double &v, const double &P) override;
    double v_PT_a(const double &P,const double &T) override;
    double v_Ps_a(const double &P,const double &s) override;
//    void Sonic_PT_g ( double& Psup, double& Tsup,double& vLaval,double& PLaval,double& TLaval,double& Csonic);

// MWATER1

    void SWater(const double &T, double &P, double &DPDT);
    void DWater(const double &T, double &Dl);
    void SatWater(double &T, double &p, double &dpdT, int option);

    // calcul
    double p_error(const double& p, const double& px);
    bool p_not_converged(const double& p, const double& px);
    bool s_not_converged(const double& s, const double& sx);
    bool h_not_converged(const double &h, const double &hx, const double& T);

    double alter_v(const double& v, const double& vmin=0, const double& vmax=1e30);
    void vregulate(const double& v,double &dv);
    void Tregulate(const double& T,double &dT);
    double x_fraction(const double &var, const double &liq, const double &gaz);
    double x_fraction(const double &T, const double& P,const double& var,char choice);

    double TX(const double& T);
    double F1(const double& Tx);
    double F2(const double& Tx);
    double P_f1(const double& T, const double& f1);
    double DPDT(const double &P, const double& T, const double& f1, const double& f2);
    double DPDV(const double&V,const double&T);

    double Rho_PT_g(const double &P, const double &T);
    double Cp_T_g(const double &T);
    double ViscDyn_T_g(const double &T);
    double ThermCond_T_g(const double &T);
    double Rho_T_l(const double &T);
    double Cp_T_l(const double &T);
    double ViscDyn_T_l(const double &T);
    double ThermCond_T_l(const double &T);




private:
    double h_cst,s_cst;
    //double aa,bb,mm;
    // Cpo coefficients
    double c0,c1,c2,c3,c4;
};


#endif




//#endif //FLUIDSMAIN_CPP_WATER_H


/*class Water {
private:
    double T,v, Tsat, vLiq, Gaz, hLiq, hGaz, sLiq, sGaz, h, s, P, x, ViscDyn, ViscCin,ThermCond, Prandtl, Rho, vLaval, PLaval, TLaval, Csonic,Cv, Cp, Sv, Av, Bp, Gt,DPV, DPT;
    void LiqGazW ( double& P, double& T, double& x,double& h,double& s);

    void W_Prop_Px_b ( double& P, double& x,double& v,double& Tsat,double& h,double& s);
    void W_Prop_Ps_bg ( double& P, double& s,double& v,double& T,double& h,double& x);
    void W_Prop_Ph_bg ( double& P, double& h,double& v,double& T,double& s,double& x);
    void W_Prop_PTx_bg ( double& P,double& T,double& v,double& h,double& s,double& x);
    void W_Prop_vh_bg ( double& v, double& h,double& P,double& T,double& s,double& x);
    void W_Prop_vs_bg ( double& v, double& s,double& P,double& T,double& h,double& x);
    void W_Prop_vP_bg ( double& v, double& P,double& T,double& h,double& s,double& x);



    void W_GTH_PT_g ( double& P, double& T,double& Cv,double& Cp,double& Sv,double& Av,double& Bp,double& Gt);
    void W_GTH_PT_l ( double& P, double& T,double& Cv,double& Cp,double& Sv,double& Av,double& Bp,double& Gt);
public:
    Water();
    Water(double& P,double& T,double& v,double& h,double& s,double& x);
    Water(TState State);
    void T_P_s (double& P);
    void P_T_s (double& T);
    void Prop_PT_l (double& P, double& T);
    void PThsvGazW (double& P,double& T);
    void PThsGazW (double& P,double& T);
    void LiqGazW (double& P, double& T, double& x);
    void PhTLiqW (double& P, double& h);
    void PhTGazW (double& P, double& h);
    void PsThLiqW (double& P, double& s);
    void PsThGazW (double& P, double& s);

    void PvThsGazW (double& P,double& v);
    void PvThsLiqW (double& P,double& v);
    void vTPhsGazW (double& v,double& T);
    void vTPhsLiqW (double& v,double& T);

    void W_Prop_Px_b (double& P,double& x);
    void W_Prop_Ps_bg (double& P,double& s);
    void W_Prop_Ph_bg (double& P,double& h);
    void W_Prop_PTx_bg (double& P);
    void W_Prop_vh_bg (double& v,double& h);
    void W_Prop_vs_bg (double& v,double& s);
    void W_Prop_vP_bg (double& v,double& P);
    void W_Prop_vT_bg (double& v,double& T);

    void W_Prop_PT_g (double& P,double& T);
    void W_TransProp_PT_g (double& P,double& T);
    void W_TransProp_PT_l (double& P,double& T);
    void W_Sonic_PT_g (double& Psup,double& Tsup);
    void W_GTH_PT_g (double& P,double& T);
    void W_GTH_PT_l (double& P,double& T);

    void W_Prop_hs_g (double& h,double& s);
    void W_Prop_vh_g (double& v,double& h);
    void W_Prop_vs_g (double& v,double& s);

    double W_Psat_T (double& T);
    double W_Tsat_P (double& P);
    double W_Tsat_v (double& v);

    void W_Prop_P_s (double& P);


};


*/
