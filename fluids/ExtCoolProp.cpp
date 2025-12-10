//
// Created by cornelia.blanke on 21.02.2024.
//

#include "ExtCoolProp.h"

ExtCoolProp::ExtCoolProp(const char *iFluidName){
    FluidName = (char *) iFluidName;
    // get critical point
    vc = 1.0 / Props1SI(FluidName, "rhocrit");
    Pc = Props1SI(FluidName, "pcrit");
    Tc = Props1SI(FluidName, "Tcrit");
    Mw = Props1SI(FluidName, "molemass") * 1000;
    w = Props1SI(FluidName, "accentric");       //TODO  returns inf ??
    DipM = Props1SI(FluidName, "dipole_moment");//TODO  returns inf ??
    ExtCoolProp::init();
}


//------------------------------------
// BASE FUNCTIONS (explicit functions)

void ExtCoolProp::init() {
    double P;
    Prop_vT(vc,Tc,P,uc,hc,sc);
    Tref = 273.15;	// [K]
    href = 200.0E3;	// [J/kg]
    sref = 1.0E3;	// [J/(kg K)]
    r = Props1SI(FluidName, "gas_constant")*1e3 / Mw;
    //Zc = Pc*vc/r/Tc;
    Zc = PropsSI("Z", "T", Tc, "P", Pc, FluidName); // compressibility factor
    Tnb = T_P_s(101325);     // boiling point at atmospheric pressure
}

void ExtCoolProp::CvCp_T_o(const double &T, double &cvo, double &cpo) {
    cpo = PropsSI("C0", "T", T, "D", 1, FluidName);     //TODO returns inf ??
    cvo = cpo-r;
}
double ExtCoolProp::h_T_o(const double &T) {
    cerr << "Not implemented" << endl;
    exit(1);
}
double ExtCoolProp::s_T_o(const double &T) {
    cerr << "Not implemented" << endl;
    exit(1);
}
void ExtCoolProp::IdealState(const double &T, double &cpo, double &cvo, double &ho, double &so) {
    CvCp_T_o(T,cvo,cpo);
    ho = h_T_o(T);
    so = s_T_o(T);
}
void ExtCoolProp::dPvdPT_vT(const double &v, const double &T, double &dPdv, double &dPdT) {
    cerr << "Not implemented" << endl;
    exit(1);
}
void ExtCoolProp::CvCp_vT(const double &v, const double &T, double &cv, double &cp) {
    double rho = 1.0/v;
    cp = PropsSI("C", "T", T, "D", rho, FluidName);
    cv = PropsSI("O", "T", T, "D", rho, FluidName);
}
void ExtCoolProp::AvBpGt_vT(const double &v, const double &T, double &Av, double &Bp, double &Gt) {
    cerr << "Not implemented" << endl;
    exit(1);
}

double ExtCoolProp::P_vT(const double &v, const double &T) {
    return PropsSI("P", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::T_vP(const double &v, const double &P) {
    return PropsSI("T", "D", 1.0/v, "P", P, FluidName);
}
void ExtCoolProp::v_PT(const double &P, const double &T, double &vl, double &vg) {
    vl = 1.0 / PropsSI("D", "Q", 0.0, "P", P, FluidName);
    vg = 1.0 / PropsSI("D", "Q", 1.0, "P", P, FluidName);
}
double ExtCoolProp::u_vT(const double &v, const double &T) {
    return PropsSI("U", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::s_vT(const double &v, const double &T) {
    return PropsSI("S", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::h_vT(const double &v, const double &T) {
    return PropsSI("H", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::h_vPT(const double &v, const double &P, const double &T) {
    return h_vT(v, T);
}
void ExtCoolProp::Prop_vT(const double &v, const double &T, double &P, double &u, double &h, double &s) {
    double rho = 1.0/v;
    P = PropsSI("P", "D", rho, "T", T, FluidName);
    u = PropsSI("U", "D", rho, "T", T, FluidName);
    h = PropsSI("H", "D", rho, "T", T, FluidName);
    s = PropsSI("S", "D", rho, "T", T, FluidName);
}
double ExtCoolProp::sv_vT(	const double &v,const double &T){
    return PropsSI("speed_of_sound", "D", 1.0/v, "T", T, FluidName);
}

//------------------------------------
// SATURATION LINE FUNCTIONS
double	ExtCoolProp::T_P_s(	const double &P) {
    return PropsSI("T","P",P,"Q",0, FluidName);
}
double	ExtCoolProp::P_T_s(	const double &T) {
    return PropsSI("P","T",T,"Q",0, FluidName);
}
void ExtCoolProp::Prop_Tx_s(const double &T,const double &x,double &v,double &P,double &u,double &h,double &s){
    if (T<=Tc){
        v = 1.0 / PropsSI("D","T",T,"Q",x, FluidName);
        P = PropsSI("P","T",T,"Q",x, FluidName);
        u = PropsSI("U","T",T,"Q",x, FluidName);
        h = PropsSI("H","T",T,"Q",x, FluidName);
        s = PropsSI("S","T",T,"Q",x, FluidName);
    }
    else
        cerr <<"Warning: call to Prop_Tx_s with T > Tc ...."<<endl;
}
double ExtCoolProp::v_PT_a(const double &P,const double &T){
    double v;
    double Ps = P_T_s(T);
    if (Ps > P)
        v = v_PT_g(P,T);
    else if (Ps < P)
        v = v_PT_l(P,T);
    else {
        cerr << "Biphasic state: cannot compute from (P,T)";
        exit(1);
    }
    return v;
}
double ExtCoolProp::v_Tx_a(const double &T, const double &x){
    return 1.0 / PropsSI("D","T",T,"Q",x, FluidName);
}

double ExtCoolProp::SurfaceTension(const double &T){
    return PropsSI("SurfaceTension", "T", T, "Q", 0, FluidName);
}

//------------------------------------
// COMMON FUNCTIONS (to all zones)
void ExtCoolProp::Prop_vT_a(const double &v, const double &T,
              double &P, double &u, double &h, double &s, double &x){
    if (T<Tc)
    {
        Prop_vT(v, T, P, u, h, s);
        double vl, vg;
        v_PT(P, T, vl, vg);
        if ((vl < v) && (v < vg))
            x = PropsSI("Q", "D", 1.0/v, "T", T, FluidName);
    }
    else
        cerr << "T >= Tc" << endl;
}
double ExtCoolProp::v_Ph_a(const double &P,const double &h){
    return 1.0 / PropsSI("D", "P", P, "H", h, FluidName);
}
double ExtCoolProp::v_Ps_a(const double &P,const double &s){
    return 1.0 / PropsSI("D", "P", P, "S", s, FluidName);
}
double ExtCoolProp::T_Ps_a(const double &P,const double &s){
    return PropsSI("T", "P", P, "S", s, FluidName);
}
double ExtCoolProp::P_vT_a(	const double &v,const double &T){
    return PropsSI("P", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::u_vT_a(	const double &v,const double &T){
    return PropsSI("U", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::h_vT_a(	const double &v,const double &T){
    return PropsSI("H", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::s_vT_a(	const double &v,const double &T){
    return PropsSI("S", "D", 1.0/v, "T", T, FluidName);
}
double ExtCoolProp::T_vP_a(const double &v,const double &P) {
    return PropsSI("T", "D", 1.0/v, "P", P, FluidName);
}
double ExtCoolProp::h_vP_a(const double &v,const double &P){
    return PropsSI("H", "D", 1.0/v, "P", P, FluidName);
}
double ExtCoolProp::v_Ts_a(const double &T,const double &s){
    return 1.0 / PropsSI("D", "T", T, "S", s, FluidName);
}
double ExtCoolProp::P_Ts_a(const double &T,const double &s){
    return PropsSI("P", "T", T, "S", s, FluidName);
}
double ExtCoolProp::s_Ph_a(const double &P,const double &h){
    return PropsSI("S", "H", h, "P", P, FluidName);
}
double ExtCoolProp::h_Ps(const double &P,const double &s){
    return PropsSI("H", "S", s, "P", P, FluidName);
}
double ExtCoolProp::T_hs_a(const double &h,const double &s){
    return PropsSI("T", "H", h, "S", s, FluidName);
}
void ExtCoolProp::FlowProp_vT_a(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl){
    Ro = 1.0/v;
    CP = PropsSI("C", "T", T, "D", Ro, FluidName);
    ViscDyn = PropsSI("V", "T", T, "D", Ro, FluidName);
    ThermCond = PropsSI("L", "T", T, "D", Ro, FluidName);
    Prandtl = PropsSI("PRANDTL", "T", T, "D", Ro, FluidName);
}
void ExtCoolProp::TransProp_vT_a(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl){
    Ro = 1.0/v;
    ViscDyn = PropsSI("V", "T", T, "D", Ro, FluidName);
    ViscCin = ViscDyn/Ro;
    ThermCond = PropsSI("L", "T", T, "D", Ro, FluidName);
    Prandtl = PropsSI("PRANDTL", "T", T, "D", Ro, FluidName);
}


//------------------------------------
// VAPOR ZONE FUNCTIONS
void ExtCoolProp::TransProp_vT_g(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) {
    Ro = 1.0/v;
    ViscDyn = PropsSI("V", "T|gas", T, "D", Ro, FluidName);
    ViscCin = ViscDyn/Ro;
    ThermCond = PropsSI("L", "T|gas", T, "D", Ro, FluidName);
    Prandtl = PropsSI("PRANDTL", "T|gas", T, "D", Ro, FluidName);
}
void ExtCoolProp::TransProp_PT_g( const double &P,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) {
    Ro = PropsSI("D", "T|gas", T, "P", P, FluidName);
    ViscDyn = PropsSI("V", "T|gas", T, "P", P, FluidName);
    ViscCin = ViscDyn/Ro;
    ThermCond = PropsSI("L", "T|gas", T, "P", P, FluidName);
    Prandtl = PropsSI("PRANDTL", "T|gas", T, "P", P, FluidName);
}
void ExtCoolProp::FlowProp_vT_g(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) {
    Ro = 1.0/v;
    ViscDyn = PropsSI("V", "T|gas", T, "D", Ro, FluidName);
    ThermCond = PropsSI("L", "T|gas", T, "D", Ro, FluidName);
    CP = PropsSI("C", "T|gas", T, "D", Ro, FluidName);
    Prandtl = PropsSI("PRANDTL", "T|gas", T, "D", Ro, FluidName);
}
void ExtCoolProp::CvCp_vT_g(const double &v,const double &T,double &cv,double &cp){
    double rho = 1.0/v;
    cp = PropsSI("C", "T|gas", T, "D", rho, FluidName);
    cv = PropsSI("O", "T|gas", T, "D", rho, FluidName);
}
double ExtCoolProp::v_PT_g(const double &P,const double &T){
    return 1.0 / PropsSI("D", "T|gas", T, "P", P, FluidName);
}
void ExtCoolProp::Prop_PT_g(const double &P,const double &T,double &v,double &u,double &h,double &s){
    v = 1.0 / PropsSI("D", "P|gas", P, "T", T, FluidName);
    u = PropsSI("U", "P|gas", P, "T", T, FluidName);
    h = PropsSI("H", "P|gas", P, "T", T, FluidName);
    s = PropsSI("S", "P|gas", P, "T", T, FluidName);
}

double	ExtCoolProp::P_vT_g(const double &v, const double &T){
    double rho = 1.0/v;
    return PropsSI("P", "D|gas", rho, "T", T, FluidName);
}
double	ExtCoolProp::u_PT_g(const double &P, const double &T){
    return PropsSI("U", "P|gas", P, "T", T, FluidName);
}
double	ExtCoolProp::h_PT_g(const double &P, const double &T){
    return PropsSI("H", "P|gas", P, "T", T, FluidName);
}
double	ExtCoolProp::s_PT_g(const double &P, const double &T){
    return PropsSI("S", "P|gas", P, "T", T, FluidName);
}
double	ExtCoolProp::v_Ph_g(const double &P, const double &h){
    return 1.0 / PropsSI("D", "P|gas", P, "H", h, FluidName);
}
double	ExtCoolProp::T_hs_g(const double &h, const double &s){
    return PropsSI("T", "H|gas", h, "S", s, FluidName);
}
double	ExtCoolProp::T_vP_g(const double &v, const double &P){
    double rho = 1.0/v;
    return PropsSI("T", "D|gas", rho, "P", P, FluidName);
}


//------------------------------------
// LIQUID ZONE FUNCTIONS
void ExtCoolProp::TransProp_vT_l(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) {
    Ro = 1.0/v;
    ViscDyn = PropsSI("V", "T|liquid", T, "D", Ro, FluidName);
    ViscCin = ViscDyn/Ro;
    ThermCond = PropsSI("L", "T|liquid", T, "D", Ro, FluidName);
    Prandtl = PropsSI("PRANDTL", "T|liquid", T, "D", Ro, FluidName);
}
void ExtCoolProp::TransProp_PT_l(const double &P,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) {
    Ro = PropsSI("D", "T|liquid", T, "P", P, FluidName);
    ViscDyn = PropsSI("V", "T|liquid", T, "P", P, FluidName);
    ViscCin = ViscDyn/Ro;
    ThermCond = PropsSI("L", "T|liquid", T, "P", P, FluidName);
    Prandtl = PropsSI("PRANDTL", "T|liquid", T, "P", P, FluidName);
}
void ExtCoolProp::FlowProp_vT_l(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) {
    Ro = 1.0/v;
    ViscDyn = PropsSI("V", "T|liquid", T, "D", Ro, FluidName);
    ThermCond = PropsSI("L", "T|liquid", T, "D", Ro, FluidName);
    CP = PropsSI("C", "T|liquid", T, "D", Ro, FluidName);
    Prandtl = PropsSI("PRANDTL", "T|liquid", T, "D", Ro, FluidName);
}
void ExtCoolProp::CvCp_vT_l(const double &v,const double &T,double &cv,double &cp){
    double rho = 1.0/v;
    cp = PropsSI("C", "T|liquid", T, "D", rho, FluidName);
    cv = PropsSI("O", "T|liquid", T, "D", rho, FluidName);
}
double ExtCoolProp::v_PT_l(const double &P,const double &T){
    return 1.0 / PropsSI("D", "T|liquid", T, "P", P, FluidName);
}
void ExtCoolProp::Prop_PT_l(const double &P,const double &T,double &v,double &u,double &h,double &s){
    v = 1.0 / PropsSI("D", "P|liquid", P, "T", T, FluidName);
    u = PropsSI("U", "P|liquid", P, "T", T, FluidName);
    h = PropsSI("H", "P|liquid", P, "T", T, FluidName);
    s = PropsSI("S", "P|liquid", P, "T", T, FluidName);
}
double	ExtCoolProp::P_vT_l(const double &v,const double &T){
    double rho = 1.0/v;
    return PropsSI("P", "T|liquid", T, "D", rho, FluidName);
}
double	ExtCoolProp::T_vP_l(const double &v,const double &P){
    double rho = 1.0/v;
    return PropsSI("T", "P|liquid", P, "D", rho, FluidName);
}
double	ExtCoolProp::u_PT_l(const double &P,const double &T){
    return PropsSI("U", "P|liquid", P, "T", T, FluidName);
}
double	ExtCoolProp::h_PT_l(const double &P,const double &T){
    return PropsSI("H", "P|liquid", P, "T", T, FluidName);
}
double	ExtCoolProp::s_PT_l(const double &P,const double &T){
    return PropsSI("S", "P|liquid", P, "T", T, FluidName);
}
double	ExtCoolProp::v_Ps_l(const double &P,const double &s){
    return 1.0 / PropsSI("D", "P|liquid", P, "S", s, FluidName);
}
double	ExtCoolProp::T_Ps_l(const double &P,const double &s){
    return PropsSI("T", "P|liquid", P, "S", s, FluidName);
}