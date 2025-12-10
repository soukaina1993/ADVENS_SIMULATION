//
// Created by cornelia.blanke on 21.02.2024.
//

#ifndef FLUIDS_EXTCOOLPROP_H
#define FLUIDS_EXTCOOLPROP_H
#include "FluidCore.h"
#include "../externals/CoolProp/include/CoolProp.h"
using namespace std;
using namespace CoolProp;


class ExtCoolProp : public FluidCore {
public:
    ExtCoolProp()=default;
    explicit ExtCoolProp(const char *iFluidName);
    ~ExtCoolProp() override= default;

    //------------------------------------
    // BASE FUNCTIONS (explicit functions)

    // ideal state @ T,P=0

    void    init() override;
    void	CvCp_T_o(const double &T,double &cvo,double &cpo) override;
    double	h_T_o(	const double &T) override;
    double	s_T_o(	const double &T) override;
    void	IdealState(const double &T,double &cpo,double &cvo,double &ho,double &so) override;

    // EOS partial diff.
    void	dPvdPT_vT(const double &v,const double &T,double &dPdv,double &dPdT) override;          //TODO

    // thermal factors
    void	CvCp_vT(const double &v,const double &T,double &cv,double &cp) override;
    void	AvBpGt_vT(const double &v,const double &T,double &Av,double &Bp,double &Gt) override;   //TODO

    // vPT function (calcul de la pression, du volume massique ou de la temperature)
    double	P_vT(const double &v,const double &T) override;
    double	T_vP(const double &v,const double &P) override;
    void	v_PT(const double &P,const double &T,double &vl,double &vg) override;

    // u h s function
    double	u_vT(const double &v,const double &T) override;
    double	s_vT(const double &v,const double &T) override;
    double	h_vT(const double &v,const double &T) override;
    double	h_vPT(const double &v,const double &P,const double &T) override;

    void	Prop_vT(const double &v,const double &T,double &P,double &u,double &h,double &s) override;

    // sonic velocity
    double sv_vT(	const double &v,const double &T) override;


    //------------------------------------
    // SATURATION LINE FUNCTIONS
    double	T_P_s(	const double &P) override;
    double	P_T_s(	const double &T) override;

    void	Prop_Tx_s(	const double &T,const double &x,double &v,double &P,double &u,double &h,double &s) override;

    double v_PT_a(const double &P,const double &T) override;
    double v_Tx_a(const double &T, const double &x) override;

    double SurfaceTension(const double &T) override;

    //------------------------------------
    // COMMON FUNCTIONS (to all zones)
    void Prop_vT_a(const double &v,const double &T,double &P,double &u,double &h,double &s,double &x) override;
    double	v_Ph_a(	const double &P,const double &h) override;
    double	v_Ps_a(const double &P,const double &s) override;
    double	T_Ps_a(	const double &P,const double &s) override;

    double	P_vT_a(	const double &v,const double &T) override;
    double	u_vT_a(	const double &v,const double &T) override;
    double	h_vT_a(	const double &v,const double &T) override;
    double	s_vT_a(	const double &v,const double &T) override;
    double	T_vP_a(	const double &v,const double &P) override;
    double	h_vP_a(	const double &v,const double &P) override;
    double	v_Ts_a(	const double &T,const double &s) override;
    double	P_Ts_a(	const double &T,const double &s) override;
    double	s_Ph_a(	const double &P,const double &h) override;
    double	h_Ps(	const double &P,const double &s) override;
    double	T_hs_a(	const double &h,const double &s) override;

    void	FlowProp_vT_a(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) override;
    void	TransProp_vT_a(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) override;


    //------------------------------------
    // VAPOR ZONE FUNCTIONS
    void TransProp_vT_g(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) override;
    void TransProp_PT_g( const double &P,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) override;
    void FlowProp_vT_g(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) override;
    void CvCp_vT_g(	const double &v,const double &T,double &cv,double &cp) override;
    double v_PT_g(	const double &P,const double &T) override;
    void Prop_PT_g(	const double &P,const double &T,double &v,double &u,double &h,double &s) override;

    double	P_vT_g(	const double &v,const double &T) override;
    double	u_PT_g(	const double &P,const double &T) override;
    double	h_PT_g(	const double &P,const double &T) override;
    double	s_PT_g(	const double &P,const double &T) override;
    double	v_Ph_g(	const double &P,const double &h) override;
    double	T_hs_g(	const double &h,const double &s) override;
    double	T_vP_g(	const double &v,const double &P) override;



    //------------------------------------
    // LIQUID ZONE FUNCTIONS
    void TransProp_vT_l(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) override;
    void TransProp_PT_l(const double &P,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) override;
    void FlowProp_vT_l(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) override;
    void CvCp_vT_l(const double &v,const double &T,double &cv,double &cp) override;
    double v_PT_l(const double &P,const double &T) override;
    void Prop_PT_l(const double &P,const double &T,double &v,double &u,double &h,double &s) override;
    double	P_vT_l(	const double &v,const double &T) override;
    double	T_vP_l(	const double &v,const double &P) override;
    double	u_PT_l(	const double &P,const double &T) override;
    double	h_PT_l(	const double &P,const double &T) override;
    double	s_PT_l(	const double &P,const double &T) override;
    double	v_Ps_l(	const double &P,const double &s) override;
    double	T_Ps_l(	const double &P,const double &s) override;
};


#endif //FLUIDS_EXTCOOLPROP_H
