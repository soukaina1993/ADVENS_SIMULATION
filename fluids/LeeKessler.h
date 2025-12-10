//
// Created by Loan on 23.05.2022.
//

#ifndef MAIN_CPP_LEEKESSLER_H
#define MAIN_CPP_LEEKESSLER_H


#include <string>
#include "FluidCore.h"
#include "LK.h"

using namespace std;


class LeeKessler:public LK4{



public:
    struct LKFluidRec{
        string fluidName, fluidSymbole, fuildFormula, fluidSource;
        double Tfp, Tnb, PC, VC, TC, MW, ALPHA, BETA, GAMMA, DELTA, roL, TL;
        double wPitzer, Dipm;
        double ADA0, ADA1, TREF, WP, AF, ZC, HREFO, SREFO;
        double SC, HC;
    };
    LeeKessler(const char *iFluidName,const char* file_in );
    LeeKessler(const char *iFluidName);
    LeeKessler(string iFluidName);
    LeeKessler();
    LeeKessler(char *iFluidName,const double &ivc,const double &iPc,const double &iTc,const double &iTnb,const double &iMw,const double &iw,
               const double &iDipM,const double &ic0,const double &ic1,const double &ic2,const double &ic3,const double &ic4,const double &iTL,const double &iRoc,const double &iTfp);



    void init();

    void	CvCp_T_o(const double &T,double &cvo,double &cpo);

    double	h_T_o(	const double &T);
    double	s_T_o(	const double &T);
    void	IdealState(const double &T,double &cpo,double &cvo,double &ho,double &so);

    // EOS partial diff.
    void	dPvdPT_vT(const double &v,const double &T,double &dPdv,double &dPdT);

    // thermal factors
    void	CvCp_vT(const double &v,const double &T,double &cv,double &cp);
    void	AvBpGt_vT(const double &v,const double &T,double &Av,double &Bp,double &Gt);

    void CvCp_PT_g(const double &P,const double &T,double &cv,double &cp);
    void CvCp_PT_l(const double &P,const double &T,double &cv,double &cp);
    void CvCp_PT(const double &P,const double &T,double &cv,double &cp);
    // vPT function (calcul de la pression, du volume massique ou de la temperature)
    double	P_vT(const double &v,const double &T);

    double	P_vT_l(const double &v,const double &T);
    double	P_vT_g(const double &v,const double &T);

    double	T_vP(const double &v,const double &P);
    void	v_PT(const double &P,const double &T,double &vl,double &vg);
    double v_PT_l (const double & P,const double & T);
    // u h s function
    double	u_vT(const double &v,const double &T);
    double	s_vT(const double &v,const double &T);
    double	h_vT(const double &v,const double &T);
    double	h_vPT(const double &v,const double &P,const double &T);

    void	Prop_vT(const double &v,const double &T,double &P,double &u,double &h,double &s);


    // sonic velocity
//    double	sv_vT(	const double &v,const double &T);

    // SATURATION LINE FUNCTIONS
    double	T_P_s(	const double &P);
    double	P_T_s(	const double &T);
    void Prop_T_s(const double &T, double& vl,double& vg,double &P,double& ul,double& ug, double &hl, double &hg, double &sl, double &sg);
    void Prop_T_s(const double &T, double& vl,double& vg,double &P, double &hl, double &hg, double &sl, double &sg);
    void Prop_Tx_s(	const double &T,const double &x,double &v, double &P,double &u,double &h, double &s);
    void Prop_P_s(const double &P, double& vl,double& vg,double &T,double& ul,double& ug, double &hl, double &hg, double &sl, double &sg);


    // VAPOR ZONE FUNCTIONS
    void	TransProp_vT_g(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl);
    void    FlowProp_vT_g(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl);
    void Prop_PT_g (const double &P, const double & T,double&  v,double&  u, double& h, double& s);

    // LIQUID ZONE FUNCTIONS

    void    FlowProp_vT_l(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl);
    void	TransProp_vT_l(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl);
    void TransProp_vT_a(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl);
    void Prop_PT_l (const double &P, const double & T,double&  v,double&  u, double& h, double& s);

    // ALL
//    double	u_Ph_a(	const double &P);

};

/*
 *	----------------- Function InitLKFluid ----------------- 																			}
*	REFERENCE ENTHALPIES HREFO,HREF AND ENTROPIES SREFO AND SREF AS DEFINED IN
*	REPORT REFERRED TO ABOVE;INPUT DATA(APART FROM COMMON STATEMENTS) IS TREF.
*	AS THE PROGRAM IS WRITTEN, THE DATUM OF ENTROPY  AND ENTHALPY IS THE
*	SATURATED LIQUID AT 0 ¡C. AT THIS STATE ENTHALPY IS SET AT 200 KJ/KG AND
*	ENTROPY IS SET AT 1 KJ/KG K.																		}
*    CALCULATES 	AND ACENTRIC FACTOR WP

*	TNB K=NORMAL BOILING POINT
*
*	INPUT:		TREF	K
*	OUTPUT:	AF				PITZER ACENTRIC FACTOR RATIO FLUID
*				ZC		CRITICAL COMPRESSIBILITY FACTOR, ZC=PC * VC /(R*TC) , R=8314/MW
*				HREFO	KJ/KG	ENTHALPY AT TREF AND ZERO PRESSURE
*				SREFO KJ/KG K	ENTROPY AT TREF AND ZERO PRESSURE
*				VREF	M3/KG	SATURATED VAPOUR PROPERTIES AT TREF K
*				PREF	KN/M2
*				TREF	KJ/KG
*				HREF	KJ/KG
*				SREF	KJ/KG K
*    CALLS SUBROUTINES Ps_T,Ts_P,VRro_PTr_g	}

 */

#endif //MAIN_CPP_LEEKESSLER_H
