#include "HPipe.h"
#include <cmath>
#include <iostream>


using namespace std;

// CONSTRUCTEUR


HPipe::HPipe(const char *iPipeName, TInputFile &InputFile) {
    InputFile.open();

    int nbRecordsConduite = InputFile.nRecords();
    TSeqsArray RecordsConduite(InputFile.nRecords());// création du tableau VIDE avec le même nombre de Lignes que de nombre de records de InputFile
    InputFile.config(Tiofile::onrcds); //on configure le TSeqsArray en mode Records (ie ligne de résultats et pas colonnes)
    RecordsConduite = *InputFile.GetRecords(); //on ajoute les valeurs des Records de l'InputFile

    //Initialisation de la conduite
    char* ConduiteName = (char*)iPipeName;
    TSequence* PtrSeqPipe = Get(RecordsConduite, ConduiteName);
    if (!(PtrSeqPipe == NULL))
        cout << *PtrSeqPipe;

    PipeNumber = PtrSeqPipe->Get(0);
    Length= PtrSeqPipe->Get(1);
    Diam= PtrSeqPipe->Get(2);
    RugAbs= PtrSeqPipe->Get(3);
    flow=new Flow;                      // allocate new flow
    pStateIn=new TPhysicState;          // allocate new physical state
    flow_created = state_created = true;
    *flow->Massflow = PtrSeqPipe->Get(4);
}

HPipe::HPipe() {
    flow=new Flow;
    pStateIn=new TPhysicState;
    flow_created = state_created = true;
}

HPipe::HPipe(Flow *iflow) {
    flow=iflow;
    pStateIn=new TPhysicState;          // allocate new physical state
    state_created = true;
}

HPipe::HPipe(Flow_physic *iflow) {
    flow=iflow;
    pStateIn=iflow->pStateIn;           // pointer to physical state
}

HPipe::~HPipe() {
    if (flow_created) delete flow;
    if (state_created) delete pStateIn;
}

void HPipe::init() {}

BasicES &HPipe::performance() {
    return *this;
}

float HPipe::objective(TSequence &seq) {
    return 0;
}

void HPipe::ThermoPhysicCalculation() {
    pStateIn->P=flow->PIn();
    pStateIn->T=flow->TIn();
    pStateIn->x=flow->xIn();
    if (pStateIn->x > 0)
        flow->fluid->TransProp_PT_g(flow->getP(), flow->getT(), pStateIn->ViscDyn, pStateIn->ViscCin, pStateIn->Ro, pStateIn->ThermCond,pStateIn->Prandtl);
    else
        flow->fluid->TransProp_PT_l(flow->getP(), flow->getT(), pStateIn->ViscDyn, pStateIn->ViscCin, pStateIn->Ro, pStateIn->ThermCond, pStateIn->Prandtl);
}


double HPipe::MoodyFactor (double Reynolds, double roughness){
    if (Reynolds <= 0.0)
        return 1.0e30;
    else if (Reynolds <= 2320)
        return 64.0 / Reynolds;
    else {
        const double log_10 = log(10);      // log() is natural logarithm
        //double a = 2 / log(10);
        //double b = RugRel / 3.7;
        double bd = roughness * Reynolds * log_10 / 18.874;     // b * d
        //double d = Re * log(10) / 5.02;
        double log_d = log(Reynolds * log_10 / 5.02);        // log(d)
        double s = bd + log_d;
        //double q = pow(s, (s / (s + 1)));
        double log_q = log(pow(s, (s / (s + 1))));         // log(q)
        double g = bd + log_d - log_q;
        double z = log_q - log(g);
        //double dla = z * g / (g + 1);
        double dcf = z * g / (g + 1) * (1 + ((z / 2) / (pow((g + 1), 2) + ((z / 3) * (2 * g - 1)))));

        return pow(1 / (2 / log_10 * (log_d - log_q + dcf)), 2);
    }
}

double HPipe::MoodyFactor() {
    //double NewVarAux, VarAux, PrecVarAux, ErrVarAux, bb;
    Re = 4.0 * *flow->Massflow / (pStateIn->ViscDyn * M_PI * Diam);
    if (RugRel == 0.0)
        RugRel = RugAbs / Diam;
    return MoodyFactor(Re, RugRel);
}

double HPipe::Diameter() {
    double NewDiametre, Precdiametre, ErrDiametre;

    //{Initialisation de l'Algorithme :	}
    Diam = 0.3;
    Precdiametre = 0.005;

   // {Algorithme methode de la bissectrice }

    //{Calcul du coefficient de perte de charge :}

    do{
        F_Moody = MoodyFactor();
        // {calcul du  nouveau Diameter, NewDiametre:}
        NewDiametre = exp(0.2 * ln(8.0 *  Length * F_Moody * sqr(*flow->Massflow / M_PI) /  DeltaPr / pStateIn->Ro));
        /* {On regarde si la valeur de newDiametre est proche de celle de Diameter :}
         */
        ErrDiametre = abs(NewDiametre - Diam);
        Diam = NewDiametre;
    } while (ErrDiametre > Precdiametre);

    return  Diam;
}

void HPipe::InitConduiteData(double &iLength, double &iDiam, double &iRuAbs) {
    Diam = iDiam;
    Length = iLength;
    RugAbs = iRuAbs;
}

void HPipe::InitConduiteData(double &iMFlow, double &iLength, double &iDiam, double &iRuAbs, double &iDeltaP) {
    Diam = iDiam;
    Length = iLength;
    RugAbs = iRuAbs;
    DeltaPr = iDeltaP;
    *flow->Massflow = iMFlow;
}

void HPipe::Initflow(FluidCore*&iFluid, double iTin, double iPin, double iMFlow){
    flow->fluid=iFluid;
    flow->StateIn->T=iTin;
    flow->StateIn->P=iPin;
    *flow->Massflow=iMFlow;
}

void HPipe::getdataBase(string& Type,string & DN){
    string file_datapipe;
    file_datapipe = "dataBase/Pipe_";
    file_datapipe += Type+".txt";

    TInputFile  InputFile((char*)file_datapipe.c_str(), 2, 1, '\t');
    InputFile.open();
    InputFile.config(Tiofile::onrcds);
    TSeqsArray Records(InputFile.nRecords());
    Records = *InputFile.GetRecords();

    TSequence* PtrSeq = Get(Records, (char*) DN.c_str());
    Diam = PtrSeq->Get(0);
    RugAbs = PtrSeq->Get(1);
    U_pipe=PtrSeq->Get(2);

    double Diam_ext = PtrSeq->Get(3);
    double Ep_tube	 = PtrSeq->Get(4);
    if(Diam==0)
        Diam=Diam_ext-2*Ep_tube;

    if(U_pipe==0) {
        if (Diam_ext == 0)
            Diam_ext = Diam + Ep_tube * 2;

        double Diam_ext_insu = PtrSeq->Get(6);
        if (Diam_ext_insu == 0) {
            double Ep_insu = PtrSeq->Get(7);
            Diam_ext_insu = Diam_ext + Ep_insu * 2;
        }

        double U_i = PtrSeq->Get(8);
        double U_p = PtrSeq->Get(5);
        if (Ep_tube == 0)
            Ep_tube = (Diam_ext - Diam) / 2;

        U_pipe = Heattransfer_mat(Diam, Diam_ext, Diam_ext_insu, U_p, U_i); // calcul de heat transfer du material
    }
}







void HPipe::getdataBase(TSeqsArray& Records, string& DN){
        TSequence* PtrSeq = Get(Records, (char*) DN.c_str());
        Diam = PtrSeq->Get(0);
        RugAbs = PtrSeq->Get(1);
        U_pipe=PtrSeq->Get(2);

        double Diam_ext = PtrSeq->Get(3);
        double Ep_tube	 = PtrSeq->Get(4);
        if(Diam==0)     Diam=Diam_ext-2*Ep_tube;

        if(U_pipe==0){
            if(Diam_ext==0) Diam_ext=Diam+Ep_tube*2;

            double Diam_ext_insu = PtrSeq->Get(6);
            if(Diam_ext_insu==0){
                double Ep_insu = PtrSeq->Get(7);
                Diam_ext_insu=Diam_ext+Ep_insu*2;}

            double U_i  = PtrSeq->Get(8);
            double U_p  = PtrSeq->Get(5);
            if(Ep_tube==0) Ep_tube=(Diam_ext-Diam)/2;

            U_pipe=Heattransfer_mat(Diam,Diam_ext,Diam_ext_insu,U_p,U_i);} // calcul de heat transfer du material
}




void HPipe::BendData(double &F_Moody) {
    /* QUE FAIT CETTE FONCTION ??*/
    double BendR, BendA;
    double md;

    md =F_Moody;
    std::cout<<"--- BENDS DATA ---"<<std::endl;
    std::cout<<"Radius [mm]"<<std::endl;
    //readln(BendR);
    //write('Angle [� ]');
    //readln(BendA);
}

double HPipe::Debit() {
    /*RQ Pascal: La différence avec le while .. do réside dans le fait que le repeat ... until exécute toujours au moins une fois le bloc d'instructions avant d'évaluer l'expression booléenne alors que le while ...
     * do évalue immédiatement son expression booléenne avant d'exécuter le bloc d'instructions.
     */
    double newDebit;
    F_Moody = 0.02;
    newDebit = sqrt(DeltaPr * sqr(sqr(Diam)) * Diam * pStateIn->Ro * sqr(M_PI) / (8 * Length * F_Moody));
    do {
        *flow->Massflow = newDebit;
        F_Moody = MoodyFactor();
        newDebit = sqrt(DeltaPr * sqr(sqr(Diam)) * Diam * pStateIn->Ro * sqr(M_PI) / (8 * Length * F_Moody));
    } while (abs(newDebit - *flow->Massflow) / newDebit < 0.001);
    return newDebit;
}

void HPipe::PressureDropCalculation() {
    F_Moody = MoodyFactor();      // includes computation of F_Moody and Re
    // Speed = *flow->Massflow / pStateIn->Ro / (Pi * sqr(Diam) / 4);
    // DeltaPr = F_Moody / pStateIn->Ro * Length * sqr(*flow->Massflow / (Pi * sqr(Diam) / 4)) / Diam / 2;
    // simplified:
    Speed = Re * pStateIn->ViscCin / Diam;
    DeltaPr_m = F_Moody /Diam * (pStateIn->Ro /2) * sqr(Speed);
    DeltaPr = DeltaPr_m * Length;
}

void HPipe::updateFlow_P() {
    double Pout=flow->PIn()-DeltaPr;
    flow->setPout(Pout);
}


void HPipe::DiameterCalculation() {

    Diam = Diameter();
    Speed = *flow->Massflow / pStateIn->Ro / (M_PI * sqr(Diam) / 4);
}

double HPipe::Heattransfer_mat(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i){

    // return 1/(((Diam/2)*log((Diam_ext_insu/2)/(Diam_ext/2))/U_i) + ((Diam/2)*log((Diam_ext/2)/(Diam/2))/U_p));
    // simplified:
    return 2 / (Diam * (log(Diam_ext_insu/Diam_ext)/U_i + log(Diam_ext/Diam)/U_p) );
}

double HPipe::Heattransfer_mat(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &Diam_ext_mant,
                               double &U_p,double &U_i, double &U_m){
// three layers pipe including mantle
    return 2 / (Diam * (log(Diam_ext_mant/Diam_ext_insu)/U_m + log(Diam_ext_insu/Diam_ext)/U_i + log(Diam_ext/Diam)/U_p) );
}

double HPipe::Nusselt_lam(double Reynolds){
    // VDI-Waermeatlas eq. (6), laminar flow
    double Nu_m = pow(1.615 * (Reynolds * pStateIn->Prandtl * Diam / Length), 1.0 / 3.0);
    return pow(49.370896 + pow(Nu_m - 0.7, 3), 1.0 / 3.0);
}

double HPipe::Nusselt_turb(double Reynolds){
    // Gnielinski correlation
    return ((F_Moody/8) * (Reynolds-1000) * pStateIn->Prandtl) / (1 + 12.7*sqrt(F_Moody/8) * (pow(pStateIn->Prandtl,(2.0/3.0))-1));
}

double HPipe::Heattransfer_conv(){
    if (Length == 0)
        return 0.0;
    double Nu;
    F_Moody = MoodyFactor();      // includes computation of F_Moody and Re
    if (Re < 2300)
        Nu = Nusselt_lam(Re);
    else if (Re > 10000)
        Nu = Nusselt_turb(Re);
    else {
        double gamma = (Re - 2300)/7700;
        Nu = (1-gamma) * Nusselt_lam(2300) + gamma * Nusselt_turb(10000);
    }
    return Nu * pStateIn->ThermCond / Diam;     // heat transfer coeff in [W/(m2*K)]
}

double HPipe::Heattransfer_total(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i){
    double h_conv = Heattransfer_conv();
    double h_mat  = Heattransfer_mat(Diam,Diam_ext,Diam_ext_insu,U_p,U_i);
    return h_conv*h_mat / (h_conv+h_mat);
}

double HPipe::Heattransfer_total(double &Diam,double &Diam_ext,double &Diam_ext_insu, double &Diam_ext_mant,
                          double &U_p,double &U_i, double &U_m){
    double h_conv = Heattransfer_conv();
    double h_mat  = Heattransfer_mat(Diam,Diam_ext,Diam_ext_insu,Diam_ext_mant,U_p,U_i,U_m);
    return h_conv*h_mat / (h_conv+h_mat);
}

double HPipe::Heattransfer_total(){
    double h_conv = Heattransfer_conv();
    return h_conv*U_pipe / (h_conv+U_pipe);
}

void HPipe::tauCalculation(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i){
    U_total = Heattransfer_total(Diam, Diam_ext, Diam_ext_insu, U_p, U_i);
    if (flow->Cp == 0){
        cerr << "Cp is zero!" << endl;
        exit(1);
    }
    if (*flow->Massflow == 0)   // Massflow = 0 --> set tau to a very big value
        tau = 1.e+30;
    else
        tau = (U_total * Length * Diam * M_PI) / (*flow->Massflow * flow->Cp);
}


void HPipe::tauCalculation(){
    if (flow->Cp == 0){
        cerr << "Cp is zero!" << endl;
        exit(1);
    }
    if (*flow->Massflow == 0)   // Massflow = 0 --> set tau to a very big value
        tau = 1.e+30;
    else {
        U_total = Heattransfer_total();
        tau = (U_total * Length * Diam * M_PI) / (*flow->Massflow * flow->Cp);
    }
}

void HPipe::HeatlossCalculation(double Tambient) {
    tauCalculation();
    DeltaT = (flow->TIn() - Tambient) * (1 - exp(-tau));
}

void HPipe::updateFlow_T(){
    double Tout=flow->TIn()-DeltaT;
    flow->setTout(Tout);
}

void HPipe::LossCalculation(double Tambient){
    PressureDropCalculation();
    updateFlow_P();
    HeatlossCalculation(Tambient);
    updateFlow_T();
}


void HPipe::Display(){

    cout<<"Diametre                       [mm]   :       "<<Diam*1000<<endl;
    cout<<"PipeLength                     [mm]   :       "<<Length*1E3<<endl;
    cout<<"Rugosite absolue               [mm]   :       "<<RugAbs*1000<<endl;

    cout<<"Mass flow rate               [kg/s]   :       "<<*flow->Massflow<<endl;
    cout<<"Vol. flow                    [m3/h]   :       "<<*flow->Massflow/pStateIn->Ro * 3.6E3<<endl;
    cout<<"Masse volumique             [kg/m3]   :       "<<pStateIn->Ro<<endl;
    cout<<"Viscosite dynamique         [mPa s]   :       "<<pStateIn->ViscDyn * 1E3<<endl;

    cout<<"Reynolds number                 [-]   :       "<<Re<<endl;
    cout<<"Facteur de pertes de charge     [-]   :       "<<F_Moody<<endl;
    cout<<"Average speed                 [m/s]   :       "<<Speed<<endl;
    cout<<"Pressure drop                 [bar]   :       "<<DeltaPr / 1E5<<endl;

}

double HPipe::getLength() {
    return Length;
}

double HPipe::getDiameter() {
    return Diam;
}

double HPipe::getRugAbs() {
    return RugAbs;
}

double HPipe::getDeltaPr() {
    return DeltaPr;
}

double HPipe::getMassflow() {
    return *flow->Massflow;
}

double HPipe::getRe() {
    return Re;
}

double HPipe::getPrandtl() {
    return pStateIn->Prandtl;
}

double HPipe::getSpeed() {
    return Speed;
}

double HPipe::getViscDyn() {
    return pStateIn->ViscDyn;
}

double HPipe::getViscCin() {
    return pStateIn->ViscCin;
}

double HPipe::getConduct() {
    return pStateIn->ThermCond;
}

double HPipe::getRo() {
    return pStateIn->Ro;
}

double HPipe::getF_Moody() {
    return F_Moody;
}

void HPipe::setLength(double &length) {
    Length=length;

}

void HPipe::setDiameter(double &diameter) {
    Diam=diameter;
}

void HPipe::setRugAbs(double &rugabs) {
    RugAbs=rugabs;
}

void HPipe::setMassflow(double massflow) {
    *flow->Massflow=massflow;
}

void HPipe::setViscDyn(double &viscdyn) {
    pStateIn->ViscDyn=viscdyn;
}

void HPipe::setViscCin(double &visccin) {
    pStateIn->ViscCin=visccin;
}

void HPipe::setConduct(double &conduct) {
    pStateIn->ThermCond=conduct;
}

void HPipe::setRo(double &iro) {
pStateIn->Ro=iro;
}


double HPipe::getT(){
    return flow->getTlog();

}
double HPipe::getP(){
    return flow->getP();
}
