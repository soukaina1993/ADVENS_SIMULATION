#include "HPipe.h"
#include "Hydraulic.h"
#include <cmath>
#include <iostream>

using namespace std;

const double Pi = 3.14159;
const double g = 9.81;

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
    MassFlow = PtrSeqPipe->Get(4);
}



void HPipe::init() {}

BasicES &HPipe::performance() {
    return *this;
}

float HPipe::objective(TSequence &seq) {
    return 0;
}

void HPipe::ThermoPhysicCalculation(const TState &Point, PRfluid& fluid) {
    if (Point.x > 0) {
        fluid.TransProp_PT_g(Point.P, Point.T, ViscDyn, ViscCin, Ro, Conduct,Prandtl);
    }
    else
    {
        fluid.TransProp_PT_l(Point.P, Point.T, ViscDyn, ViscCin, Ro, Conduct, Prandtl);
    }

}

double HPipe::MoodyFactor() {
    double NewVarAux, RugRel, VarAux, PrecVarAux, ErrVarAux, bb;
       Re = 4.0 * MassFlow / (ViscDyn * Pi * Diam);
        if (Re <= 2320)
            {return 64 / Re;}
        else
            RugRel = RugAbs / Diam;
        //Initialisation algorithme
            VarAux = 8.0;
            //{La precision est donn�e sur VarAux : }
            PrecVarAux = 0.01;

           /* {Algorithme Methode de la bissectrice}
            {Calcul de NewVarAux :	}*/
            ErrVarAux = 1.0;

                do
                {
                    bb = exp(1.03 * ln(RugRel / 3.71)) + exp(1.1 * ln(4.26 / Re)) * exp(1.1 * ln(VarAux));
                    NewVarAux = -1.94 * ln(bb) / ln(10.0);
                    //{On regarde si la valeur de NewVarAux est proche de celle de VaeAux : }
                    ErrVarAux = abs(NewVarAux - VarAux);
                    VarAux = NewVarAux;
                }while (ErrVarAux > PrecVarAux);

            return 1.0 / sqrt(NewVarAux);
}

double HPipe::Diameter() {
    double NewDiametre, CoeffPr, Precdiametre, ErrDiametre, ee;

    //{Initialisation de l'Algorithme :	}
    Diam = 0.3;
        Precdiametre = 0.005;

       // {Algorithme methode de la bissectrice }

        //{Calcul du coefficient de perte de charge :}

        ErrDiametre = 1.0;
            do{
                CoeffPr = this->MoodyFactor();
                // {calcul du  nouveau Diameter, NewDiametre:}
                ee = 1.0 / 5.0;
                NewDiametre = exp(ee * ln(8.0 *  Length * CoeffPr * sqrt(MassFlow / Pi) /  DeltaPr / Ro));
                /* {On regarde si la valeur de newDiametre est proche de celle de Diameter :}
                 */
                ErrDiametre = abs(NewDiametre - Diam);
                Diam = NewDiametre;

            } while (ErrDiametre > Precdiametre);
        F_Moody = CoeffPr;
        return  Diam;
}

void HPipe::InitConduiteData(const double &iMFlow, const double &iLenght, const double &iDiam,
                            const double &iRuAbs, const double &iDeltaP) {
    Diam = iDiam;
    Length = iLenght;
    RugAbs = iRuAbs;
    DeltaPr = iDeltaP;
    MassFlow = iMFlow;

}

void HPipe::getdataBase(string& Type,string & DN){
        string file_datapipe="dataBase/Pipe_" + Type+".txt";

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

        U_pipe=Heattransfert_mat(Diam,Diam_ext,Diam_ext_insu,U_p,U_i);} // calcul de heat transfer du material
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

            U_pipe=Heattransfert_mat(Diam,Diam_ext,Diam_ext_insu,U_p,U_i);} // calcul de heat transfer du material
}




void HPipe::BendData(const double &F_Moody) {
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
    newDebit = sqrt(DeltaPr * sqrt(sqrt(Diam)) * Diam * Ro * sqrt(Pi) / (8 * Length * F_Moody));
    do {
        MassFlow = newDebit;
        F_Moody = MoodyFactor();
        newDebit = sqrt(DeltaPr * sqrt(sqrt(Diam)) * Diam * Ro * sqrt(Pi) / (8 * Length * F_Moody));
    } while (abs(newDebit - MassFlow) / newDebit < 0.001);
    return newDebit;
}

void HPipe::PressureDropCalculation(const TState &Point, PRfluid& fluid) {
    ThermoPhysicCalculation (Point, fluid);
    F_Moody = MoodyFactor();
   Speed = MassFlow / Ro / (Pi * sqrt(Diam) / 4);
    DeltaPr = F_Moody / Ro * Length * sqr(MassFlow / (Pi * sqr(Diam) / 4)) / Diam / 2;

}

void HPipe::DiameterCalculation(const TState &Point, PRfluid& fluid) {
    ThermoPhysicCalculation (Point, fluid);
    Diam = Diameter();
    Speed = MassFlow / Ro / (Pi * sqrt(Diam) / 4);
}

double HPipe::Heattransfert_mat(double &Diam,double &Diam_ext,double &Diam_ext_insu,double &U_p,double &U_i){

     return 1/(((Diam/2)*log((Diam_ext_insu/2)/(Diam_ext/2))/U_i) + ((Diam/2)*log((Diam_ext/2)/(Diam/2))/U_p));
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
    return MassFlow;
}

double HPipe::getRe() {
    return Re;
}

double HPipe::getPrandtl() {
    return Prandtl;
}

double HPipe::getSpeed() {
    return Speed;
}

double HPipe::getViscDyn() {
    return ViscDyn;
}

double HPipe::getViscCin() {
    return ViscCin;
}

double HPipe::getConduct() {
    return Conduct;
}

double HPipe::getRo() {
    return Ro;
}

double HPipe::getF_Moody() {
    return F_Moody;
}

void HPipe::setLength(const double &length) {
    Length=length;

}

void HPipe::setDiameter(const double &diameter) {
    Diam=diameter;
}

void HPipe::setRugAbs(const double &rugabs) {
RugAbs=rugabs;
}

void HPipe::setMassflow(const double &massflow) {
MassFlow=massflow;
}

void HPipe::setViscDyn(const double &viscdyn) {
ViscDyn=viscdyn;
}

void HPipe::setViscCin(const double &visccin) {
ViscCin=visccin;
}

void HPipe::setConduct(const double &conduct) {
Conduct=conduct;
}

void HPipe::setRo(const double &ro) {
Ro=ro;
}



