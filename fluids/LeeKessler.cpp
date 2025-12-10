//
// Created by Loan on 23.05.2022.
//

#include "LeeKessler.h"


LeeKessler::LeeKessler(const char *iFluidName,const char* file_in ) {
    const char sep = '\t';
    const int nbTitleLines = 2; //nombre de lignes avec des titres
    const int nbTitleColumns = 1; //nombre de colonnes avec des titres ie sans chiffre ie

    TInputFile  InputFile((char*)file_in, nbTitleLines, nbTitleColumns, sep);
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
    w = PtrSeq->Get(5);
    DipM = PtrSeq->Get(6);

    // cpo coefficients
    ALPHA = PtrSeq->Get(7)/1000;
    BETA = PtrSeq->Get(8)/1000;
    GAMMA = PtrSeq->Get(9)/1000;
    DELTA = PtrSeq->Get(10)/1000;
    roL = PtrSeq->Get(11); //pas de 5eme coefficent cp dans le fichier d'entrée
    TL = PtrSeq->Get(12);
    Roc = PtrSeq->Get(13);
    Tfp = PtrSeq->Get(14);
    // reference state default initialisation
    // liquid line
    Tref = 273.15;	// [K]
    href = 200.0E3;	// [J/kg]
    sref = 1.0E3;	// [J/(kg K)]
    init();
};
LeeKessler::LeeKessler(const char *iFluidName) {
    Readfluid(iFluidName);
    init();
};

LeeKessler::LeeKessler(string iFluidName){
    int n = iFluidName.length();
    char iFluidNamechar[n + 1];
    strcpy(iFluidNamechar, iFluidName.c_str());

    Readfluid(iFluidNamechar);
    init();
}

LeeKessler::LeeKessler() {

};

LeeKessler::LeeKessler(char *iFluidName,const double &ivc,const double &iPc,const double &iTc,const double &iTnb,const double &iMw,const double &iw,
                       const double &iDipM,const double &ic0,const double &ic1,const double &ic2,const double &ic3,const double &ic4,const double &iTL,const double &iRoc,const double &iTfp){
    FluidName = strdup(iFluidName);
    vc = ivc;
    Pc = iPc;
    kPc=iPc/1e3;// unitŽ de Pc, kN/m2
    Tc = iTc;
    Tnb = iTnb;
    Mw = iMw;
    w = iw;
    DipM = iDipM;

    // cpo coefficients
    ALPHA = ic0/1000;
    BETA = ic1/1000;
    GAMMA = ic2/1000;
    DELTA = ic3/1000;
    roL = ic4;

    TL = iTL;
    Roc =iRoc;
    Tfp = iTfp;
    // reference state default initialisation
    // liquid line
    Tref = 273.15;	// [K]
    href = 200.0E3;	// [J/kg]
    sref = 1.0E3;	// [J/(kg K)]
    init();
}

void LeeKessler::init() {
    Init_LKConst();
    InitLKFluid();
};

void	LeeKessler::CvCp_T_o(const double &TR,double &cvo,double &cpo) {
    double T = TR * Tc;
    cpo = ((DELTA * T + GAMMA) * T + BETA) * T + ALPHA;
    cvo = cpo - LKconst.RBAR / Mw;
};

double	LeeKessler::h_T_o(	const double &TR) {
    double T = TR * Tc;
    return T*(ALPHA+T*(BETA/2.0+T*(GAMMA/3.0+DELTA/4.0*T)))-Tref*(ALPHA+Tref*(BETA/2.0+Tref*(GAMMA/3.0+DELTA/4.0*Tref)))+ HREFO;

};
double	LeeKessler::s_T_o(	const double &TR) {
    double T = TR * Tc;
    return (sqrt(T) * T - sqr(Tref) * Tref) * DELTA / 3.0 + (sqrt(T) - sqrt(Tref)) * GAMMA / 2.0 + ln(T / Tref) * ALPHA + (T - Tref) * BETA + SREFO;

};
void	LeeKessler::IdealState(const double &T,double &cpo,double &cvo,double &ho,double &so) {
    CvCp_T_o(T,cvo,cpo);
    ho = h_T_o(T);
    so = s_T_o(T);
};

// EOS partial diff.
void	LeeKessler::dPvdPT_vT(const double &v,const double &T,double &dPdv,double &dPdT) {
    cout<<"Not implement"<<endl;
};

// thermal factors
void	LeeKessler::CvCp_vT(const double &v,const double &T,double &cv,double &cp) {
    double vl, vg, Psat, ul, ug, hl, hg, sl, sg;
    Prop_T_s(T, vl, vg, Psat, hl, hg, sl, sg);
    double P = P_vT(v, T);
    if (v <= vl) {
        //CvCp_PT_l(v, T, cv, cp);
        CvCp_PT_l(P, T, cv, cp);
    } else if (v >= vg) {
        //CvCp_PT_g(v, T, cv, cp);
        CvCp_PT_g(P, T, cv, cp);

    } else {
        double cv_g, cp_g, cv_l, cp_l, x;
        x = (v - vl) / (vg - vl);
        CvCp_vT_g(vg, T, cv_g, cp_g);
        CvCp_vT_l(vl, T, cv_l, cp_l);
        cv = cv_g * x + cv_l * (1 - x);
        cp = cp_g * x + cp_l * (1 - x);
    }
}


void LeeKessler::CvCp_PT_g(const double &P,const double &T,double &Cv,double &Cp){
    double Sv,Av,Bp,Gt;
    GTH_PT_g(P, T, Cv, Cp, Sv, Av, Bp, Gt);
}

void LeeKessler::CvCp_PT_l(const double &P,const double &T,double &Cv,double &Cp){
    double Sv,Av,Bp,Gt;
    GTH_PT_l(P, T, Cv, Cp, Sv, Av, Bp, Gt);
}


void LeeKessler::CvCp_PT(const double &P,const double &T,double &cv,double &cp){
    double Psat= P_T_s(T);

    if(Psat<P){
        CvCp_PT_l(P,T,cv,cp);
    }
    else if(Psat>P){
        CvCp_PT_g(P,T,cv,cp);
    }
    else
        cout<<"biphasique"<<endl;
}

void	LeeKessler::AvBpGt_vT(const double &v,const double &T,double &Av,double &Bp,double &Gt) {
    double vl, vg, Psat, ul, ug, hl, hg, sl, sg,Cv,Cp,Sv;
    Prop_T_s(T, vl, vg, Psat, hl, hg, sl, sg);
    double P = P_vT(v, T);
    if (v < vl) {
        GTH_PT_l(P, T, Cv, Cp, Sv, Av, Bp, Gt);
    } else if (v > vg) {
        GTH_PT_g(P, T, Cv, Cp, Sv, Av, Bp, Gt);

    } else {
        cout<<"Not implement"<<endl;
    };
};



// vPT function (calcul de la pression, du volume massique ou de la temperature)
double	LeeKessler::P_vT(const double &v,const double &T) {
    double vl, vg, Tsat, ul, ug, hl, hg, sl, sg,Psat;
    Ps_T(Psat,T);
    Psat=Psat*1000;

    LK4::Prop_P_s(Psat, vl, vg, Tsat, hl, hg, sl, sg);
    if (v < vl) {
        return P_vT_l(v,T);
    } else if (v > vg) {
        return P_vT_g(v,T);
    } else {
        return Psat;
    };
};
double	LeeKessler::T_vP(const double &v,const double &P) {
    double vl, vg,Tsat, ul, ug, hl, hg, sl, sg,DPT;
    LK4::Prop_P_s(P, vl, vg, Tsat, hl, hg, sl, sg);
    if (v < vl) {
        return T_vP_l(v,Tsat);
    } else if (v > vg) {
        return T_vP_g(v,Tsat);
    } else {
        return Tsat;
    };
};
void	LeeKessler::v_PT(const double &P,const double &T,double &vl,double &vg) {
    double Tsat;
    double ul,ug,hl,hg,sl,sg;
    LK4::Prop_P_s(P, vl,vg , Tsat, hl, hg, sl, sg);

};

// u h s function
double	LeeKessler::u_vT(const double &v,const double &T) {
    double P,u,h,s,x;
    Prop_vT_a(v,T,P, u, h,s,x);
    return u;

};
double	LeeKessler::s_vT(const double &v,const double &T) {
    double P,u,h,s,x;
    Prop_vT_a(v,T,P, u, h,s,x);
    return s;
};
double	LeeKessler::h_vT(const double &v,const double &T) {
    double P,u,h,s,x;
    Prop_vT_a(v,T,P, u, h,s,x);
    return h;
};
double	LeeKessler::h_vPT(const double &v,const double &P,const double &T) {
    double u,h,s,x;
    double Psat=P;
    Prop_vT_a(v,T,Psat, u, h,s,x);
    return h;
};

void	LeeKessler::Prop_vT(const double &v,const double &T,double &P,double &u,double &h,double &s) {
    double x;
        Prop_vT_a(v,T,P, u, h,s,x);
};


// sonic velocity


// SATURATION LINE FUNCTIONS
double	LeeKessler::T_P_s(	const double &P) {
    return Tsat_P (P);
};
double	LeeKessler::P_T_s(	const double &T) {
    return Psat_T(T);
};


void LeeKessler::Prop_T_s(const double &Tsat, double& vl,double& vg,double &P, double &hl, double &hg, double &sl, double &sg){
    double T;
    P=Psat_T(Tsat);
    LK4::Prop_P_s(P, vl,vg , T, hl, hg, sl, sg);
}

void LeeKessler::Prop_T_s(const double &T, double& vl,double& vg,double &P,double& ul,double& ug, double &hl, double &hg, double &sl, double &sg){
    Prop_T_s(T, vl,vg , P, hl, hg, sl, sg);
    ug=u_vP_g (vg,P);
    ul=hl - vl * P;
}

void LeeKessler::Prop_P_s(const double &P, double& vl,double& vg,double &T,double& ul,double& ug, double &hl, double &hg, double &sl, double &sg){
    LK4::Prop_P_s(P, vl,vg , T, hl, hg, sl, sg);
    ug=u_vP_g (vg,P);
    ul=hl - vl * P;
}

void LeeKessler::Prop_Tx_s(	const double &T,const double &x,double &v, double &P,double &u,double &h, double &s){
    if (T<=Tc)
    {	double vl,vg,ul,ug,hl,hg,sl,sg;
        Prop_T_s(T, vl,vg , P,ul,ug, hl, hg, sl, sg);
        v = vl + x*(vg-vl);
        u = ul + x*(ug-ul);
        h = hl + x*(hg-hl);
        s = sl + x*(sg-sl);
    }
    else	cerr <<"Warning: call to Prop_Tx_s with T > Tc ...."<<endl;
}


// VAPOR ZONE FUNCTIONS


void	LeeKessler::TransProp_vT_g(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) {
    double P;
    P= P_vT(v,T);
    TransProp_PT_g( P, T,  ViscDyn,  ViscCin,  Ro,  ThermCond,  Prandtl);
};

void    LeeKessler::FlowProp_vT_g(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) {
    double P,ViscCin,CV;
    P= P_vT(v,T);
    CvCp_vT(v,T,CV,CP);
    TransProp_PT_g( P, T,  ViscDyn,  ViscCin,  Ro,  ThermCond,  Prandtl);
};

void LeeKessler::Prop_PT_g (const double &P, const double & T,double&  v,double&  u, double& h, double& s){
    LK4::Prop_PT_g(P,T,v,h,s);
    u=h - v * P;
}


// LIQUID ZONE FUNCTIONS

void    LeeKessler::FlowProp_vT_l(const double &v,const double &T,double &CP,double &ViscDyn,double &Ro,double &ThermCond,double &Prandtl) {
    double P,ViscCin,CV;
    P= P_vT_l(v,T);
    CvCp_vT(v,T,CV,CP);
    TransProp_PT_l( P, T,  ViscDyn,  ViscCin,  Ro,  ThermCond,  Prandtl);
};
void	LeeKessler::TransProp_vT_l(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl) {
    LK_TransProp_vT_l( v, T,  ViscDyn,  ViscCin,  Ro,  ThermCond,  Prandtl);
};



void LeeKessler::TransProp_vT_a(const double &v,const double &T,double &ViscDyn,double &ViscCin,double &Ro,double &ThermCond,double &Prandtl){

    double P = P_T_s(T);
    double vl,vg;
    v_PT(P,T,vl,vg);
    if (v<vl){
        TransProp_vT_l(v,T,ViscDyn,ViscCin,Ro,ThermCond,Prandtl);
    }
    else if (v>vg){
        TransProp_vT_g(v,T,ViscDyn,ViscCin,Ro,ThermCond,Prandtl);
    }
    else{
        double ViscDyng,ViscCing,Rog,ThermCondg,Prandtlg,ViscDynl,ViscCinl,Rol,ThermCondl,Prandtll;

        TransProp_vT_g(P,T,ViscDyng,ViscCing,Rog,ThermCondg,Prandtlg);
        TransProp_vT_l(P,T,ViscDynl,ViscCinl,Rol,ThermCondl,Prandtll);

        double x=Mfraction(vl,vg,v);
        ViscDyn=ViscDynl*(1-x)+ViscDyng*(x);
        ViscCin=ViscCinl*(1-x)+ViscCing*(x);
        Ro=Rol*(1-x)+Rog*(x);
        ThermCond=ThermCondl*(1-x)+ThermCondg*(x);
        Prandtl=Prandtll*(1-x)+Prandtlg*(x);}
}

void LeeKessler::Prop_PT_l (const double &P, const double & T,double&  v,double&  u, double& h, double& s){
    LK4::Prop_PT_l(P,T,v,h,s);
    u=h - v * P;
}


double	u_Phv_a(const double &P,const double &h,const double &v){
    return h - v * P;
}



double LeeKessler::P_vT_l(const double &v, const double &T) {

    double PR;
    double TR = T / Tc;
    double VR = v/ vc *Zc;
    Pr_vTr_l(VR,PR,TR);

    return PR*kPc*1000;
}

double LeeKessler::P_vT_g(const double &v, const double &T) {

    double PR;
    double TR = T / Tc;
    double VR = v/ vc *Zc;
    Pr_vTr_g(VR,PR,TR);

    return PR*kPc*1000;
}


double LeeKessler::v_PT_l (const double & P,const double & T) {
// TODO il y a un probleme avec cette fonction
double v,u,h,s;
    Prop_PT_l(P,T,v,u,h,s);
    return v;
}