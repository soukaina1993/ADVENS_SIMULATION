//
// Created by robin.lemaire on 21.06.2021.
//

#include <cmath>
#include "TExchanger.h"

using namespace std;

gao 	 * TExchanger::ga 	   = new gao ("exchdata.txt", "exchrslts.txt");
/****************************************************************************/


/****************************************************************************/

TExchanger::TExchanger(TExchanger::Exchanger tp	, TExchanger:: Flowtype flw, TExchanger::Type sp)
{
    flow (flw);
    type (sp);
    ex(tp);

    bcga		= TExchanger::ga;	//default ga
    ivbles	= new TSequence (TExchanger::ga->data->genome());
    idxtp	= (unsigned int) tp;
    idxsp	= (unsigned int) sp;
    //On lit les parametres constructeur

init();
}

TExchanger::TExchanger(FluidCore *&fluidin,TExchanger::Exchanger tp	, TExchanger:: Flowtype flw, TExchanger::Type sp):fluidinside(fluidin)
{
    flow (flw);
    type (sp);

    bcga		= TExchanger::ga;	//default ga
    ivbles	= new TSequence (TExchanger::ga->data->genome());
    idxtp	= (unsigned int) tp;
    idxsp	= (unsigned int) sp;
    //On lit les parametres constructeur

    init();
}


void  TExchanger::optimize  (TExchanger::Criteria ob, gao::type GA) {
    obj = ob;
    if (obj != TExchanger::COST)
        bcga->maximize(*this, GA);
    else
        bcga->minimize(*this, GA);
}


/****************************************************************************/

void TExchanger::copy (TExchanger & src)
{
    Hot		=	src.Hot;
    Cold		=	src.Cold;
    Pinch		=	src.Pinch;
    Cost		=	src.Cost;
    Q		    =	src.Q;
    Area		=	src.Area;
    Cr	=	src.Cr;
    NUT		=	src.NUT;
    Efficiency    =   src.Efficiency;
}

/****************************************************************************/

void TExchanger::deltaT() {
    if (flw == 0) {
        LeftDeltaT = abs(Hot.StateIn->T - Cold.StateIn->T); //co-courant
        RightDeltaT = abs(Hot.StateOut->T - Cold.StateOut->T); //co-courant
    } else if (flw == 1) {
        LeftDeltaT = abs(Hot.StateIn->T - Cold.StateOut->T); //contre-courant
        RightDeltaT = abs(Hot.StateOut->T - Cold.StateIn->T); //contre-courant
    }
}


void TExchanger::hexregister (){
    switch (exchtp){
        case TExchanger::ANYTYPE:   {
            DefaultParams();
            break;
        }
        case TExchanger::P1:	{
            RegisterParams("P1");
            break;
        }
        case TExchanger::EV:	{
            RegisterParams("EV");
            break;
        }
    }
}

/****************************************************************************/

double TExchanger::getPinch() {
    deltaT();
//    if (LeftDeltaT<RightDeltaT)
//        Pinch	= LeftDeltaT;
//    else
//        Pinch	= RightDeltaT;
    Pinch = min(LeftDeltaT, RightDeltaT);
    return Pinch;
}

/****************************************************************************/

double TExchanger::deltaTlm(){

    deltaT();
    if(LeftDeltaT==RightDeltaT)
        DeltaTlog=LeftDeltaT;
    else
        DeltaTlog=(LeftDeltaT-RightDeltaT)/ln(LeftDeltaT/RightDeltaT);

    if(!(flw==0 || flw ==1))
    {
        cout<<"Echangeur à courant non parallèle, facteur de correction F appliqué"<<endl;
            if(Q==0)
                getPower();
            if(NUT==0)
                getNTU();
            if(Cr==0)
                getCr();

        double F;
        if (Cr==1)
            F=Q/(NUT*(1-Q));
        else
            F=1/(NUT*(1-Cr))*log((1-Cr*Q)/(1-Q));
        DeltaTlog *= F;
    }
    if (DeltaTlog==0){
        cout<<"Avertissement, deltaTlm=0 "<< endl;
    }
    return DeltaTlog;
}

/****************************************************************************/

double TExchanger::getSurface()
{
    getPower();
    deltaTlm();
    if (flw==0 || flw==1)
    {
        if (k != 0)
            Area=Q/k/DeltaTlog;
        else
        {
            double Ch;
            double Cc;
            Ch = *Hot.Massflow * Hot.SpecHeat;
            Cc = *Cold.Massflow * Cold.SpecHeat;
            //TODO
        }
    }
    else if (flw==2)
        getSurfaceByNTU();
    return Area;
}

/****************************************************************************/

double TExchanger::getPower()
{
    if(Cold.StateOut->T==0 && Hot.StateOut->T==0)
    {
        cerr << "Temperatures de sortie non definies, impossible de calculer la puissance" << endl;
        exit(1);
    }
    else {
        Q = *Cold.Massflow * Cold.SpecHeat * abs(Cold.StateIn->T - Cold.StateOut->T);

        if (Q == 0) {
            Q = *Hot.Massflow * Hot.SpecHeat * abs(Hot.StateIn->T - Hot.StateOut->T);

            if (Q == 0) {
                //    cout << "Not enough data" << endl;
            }
        }
    }
    return Q;
}

/****************************************************************************/

void TExchanger::getSpecHeat()
{
    Q=getPower();
    if (Q==0)
    {
        cerr << "Puissance non calculee" <<endl;
    }

    Hot.SpecHeat=Q/ *Hot.Massflow/(Hot.StateIn->T-Hot.StateOut->T);
    Cold.SpecHeat=Q/ *Cold.Massflow/(Cold.StateOut->T-Cold.StateIn->T);
}

/****************************************************************************/

void TExchanger::getMassFlow()
{
    if(Q==0){
        Q=getPower();
        if (Q==0)
        {
            cerr << "Puissance non calculee, impossible de calculer le débit" << endl;
        }
    }
    if (!(Hot.StateOut->T==0 && Cold.StateOut->T==0))
    {
        if (*Hot.Massflow==0)
        {
            *Hot.Massflow=abs(Q/Hot.SpecHeat/(Hot.StateIn->T-Hot.StateOut->T));

        }
        if (*Cold.Massflow==0)
        {
            *Cold.Massflow = abs(Q / Cold.SpecHeat / (Cold.StateOut->T - Cold.StateIn->T));
        }
    }
    else
    {
        cerr << "Températures de sortie non définies, impossible de calculer le débit"<< endl;
    }
}

/****************************************************************************/


double TExchanger::getThi()
{
    Q=getPower();
    if (Q==0)
    {
        cout << "Puissance non calculee" <<endl;
    }
    Hot.StateIn->T = Q/ *Hot.Massflow/Hot.SpecHeat+Hot.StateOut->T;
    return Hot.StateIn->T;
}

/****************************************************************************/

double TExchanger::getTho()
{
    if(*Hot.Massflow==0 && Area !=0)
    {
        Hot.StateOut->T=Hot.StateIn->T;
        double var=0.0;

        while(abs(Q-var)>1){
            Hot.StateOut->T=Hot.StateOut->T-0.1;
            deltaTlm();
            var=k*Area*DeltaTlog;
        }
    }
    else {
        if(Q==0){
            Q = getPower();
            if (Q == 0) {
                cout << "Puissance non calculee" << endl;
            }}
        Hot.StateOut->T = -(Q / *Hot.Massflow / Hot.SpecHeat - Hot.StateIn->T);
    }
    return Hot.StateOut->T;
}

/****************************************************************************/

double TExchanger::getTci()
{
    if(Q==0){
        Q=getPower();
        if (Q==0)
        {
            cout << "Puissance non calculee" <<endl;
        }}
    Cold.StateIn->T = -(Q/ *Cold.Massflow/Cold.SpecHeat-Cold.StateOut->T);
    return Cold.StateIn->T;
}
/****************************************************************************/

double TExchanger::getTco()
{
    if(*Cold.Massflow==0 && Area !=0)
    {
        Cold.StateOut->T=Cold.StateIn->T;
        double var;
        var=0;
        /*cout <<"Thi= "<<Hot.StateIn->T<<endl;
        cout<<"Tho= "<<Hot.StateOut->T<<endl;
        cout<<"Tci= "<<Cold.StateIn->T<<endl;
        cout <<"Q= "<<Q<<endl;
        cout <<"k= "<<k<<endl;
        cout <<"Area= "<<Area<<endl;*/
        while(abs(Q-var)>1){
            Cold.StateOut->T=Cold.StateOut->T+0.1;
            //cout <<"StateOut->T= "<<Cold.StateOut->T<<endl;
            deltaTlm();
            //cout <<"dTlm= "<<DeltaTlog<<endl;
            var=k*Area*DeltaTlog;
            //cout <<"var= "<<var<<endl;
        }
    }
    else {
        if(Q==0){
            Q = getPower();
            if (Q == 0) {
                cout << "Puissance non calculee" << endl;
            }}
        else
            Cold.StateOut->T = Q / *Cold.Massflow / Cold.SpecHeat + Cold.StateIn->T;
    }
    return Cold.StateOut->T;
}

/****************************************************************************/

double TExchanger::getCr()
{
    double Ch=*Hot.Massflow*Hot.SpecHeat;
    double Cc=*Cold.Massflow*Cold.SpecHeat;

    if(Ch<Cc)
        Cr=Ch/Cc;
    else
        Cr=Cc/Ch;
    return Cr;
}


/****************************************************************************/

double TExchanger::efficiency() {
    double Ch = *Hot.Massflow * Hot.SpecHeat;
    double Cc = *Cold.Massflow * Cold.SpecHeat;

    if (Cc >= Ch)
        Efficiency = (Hot.StateIn->T - Hot.StateOut->T) / (Hot.StateIn->T - Cold.StateIn->T);
    else
        Efficiency = (Ch * (Hot.StateIn->T - Hot.StateOut->T)) / (Cc * (Hot.StateIn->T - Cold.StateIn->T));

    return Efficiency;
}

/****************************************************************************/

double TExchanger::efficiencyByNUT(){
    getCr();
    if(flw == 0) //co-courant
        Efficiency = (1-exp(-NUT*(1+Cr))) / (1+Cr);

    else if(flw == 1) {  //contre-courant
        if (Cr != 0) {
            double expNUT = exp(NUT * (Cr - 1));
            Efficiency = (1 - expNUT) / (1 - Cr * expNUT);
        }
        else
            Efficiency = NUT / (1 + NUT);
    }

    else if(flw == 2)
    {
        int typeAutre;
        cout << "Votre e        changeur est :"<<endl;
        cout << "a courant croises, deux fluides brasses (0)"<< endl;
        cout << "a courant croises, Cmin non brasse (1)"<<endl;
        cout << "a courant croises, Cmax non brasse (2)"<<endl;
        cout << "a courant croises, deux fluides non brasses (3)"<<endl;
        cout << "Autre (4)"<< endl;
        cin >> typeAutre;

        switch (typeAutre) {
            case 0:
                Efficiency=pow((1/(1-exp(-NUT))+Cr/(1-exp(-NUT*Cr))-1/NUT),-1);//Croises, deux fluides brasses
                break;
            case 1:
                Efficiency=1/Cr*(1-exp(-Cr*(1-exp(-NUT)))); //Croisés, Cmin non mixé
                break;
            case 2:
                Efficiency=1-exp(-(1/Cr)*(1-exp(-Cr*NUT)));//Croises, Cmax non mixe
                break;
            case 3:
                Efficiency=1-exp((exp(-Cr*pow(NUT,0.78))-1)/(Cr*pow(NUT,-0.22)));//Croises, deux fluides non brasses
                break;
            case 4:
                cout << "Type d echangeur pas encore modelise"<< endl;
                break;
            default:
                cout << "Type invalide pour calcul NUT" << endl;
                break;
        }
    }
    return Efficiency;
}

/****************************************************************************/

double TExchanger::getNTU()
{
    getCr();
    if(flw==0) //co-courant
    {
        NUT=(-log(1-(1+Cr)*Efficiency))/(1+Cr);
    }
    else if(flw==1)//contre-courant
    {
        if (Cr != 1)
            NUT=(-log((1-Efficiency)/(1-Cr*Efficiency)))/(1-Cr);
        else
            NUT = Efficiency / (1 - Efficiency);
    }
    else if(flw==2)
    {
        int typeAutre;
        cout << "Votre echangeur est :"<<endl;
        cout << "a courant croises, deux fluides brasses (0)"<< endl;
        cout << "a courant croises, Cmin non brasse (1)"<<endl;
        cout << "a courant croises, Cmax non brasse (2)"<<endl;
        cout << "tubulaire, 1 passe cote coque et 2n passes cote tubes (3)"<<endl;
        cout << "Autre (4)"<<endl;
        cin >> typeAutre;

        switch (typeAutre) {
            case 0:
                NUT=-log(1+log(1-Efficiency* Cr)/Cr);
                break;
            case 1:
                NUT=-log(1+(log(1-Efficiency*Cr))/Cr); //Croisés, Cmin non mixé
                break;
            case 2:
                NUT=-log(1+Cr*log(1-Efficiency))/Cr;
                break;
            case 3:
                double y;
                y=(pow(Efficiency/2,-1)-1-Cr)/pow(1+Cr*Cr,0.5);
                NUT=-log((y-1)/(y+1))/pow(1+Cr*Cr,0.5);
                break;
            case 4:
                cout << "Type d echangeur pas encore modelise"<<endl;
                break;
            default:
                cout << "Type invalide pour calcul NUT" << endl;
                break;
        }
    }
    return NUT;
}

/****************************************************************************/

double TExchanger::getSurfaceByNTU()
{
    double Ch=*Hot.Massflow*Hot.SpecHeat;
    double Cc=*Cold.Massflow*Cold.SpecHeat;

    getNTU();

    if(Ch<Cc)
        Area=NUT*Ch/k;
    else
        Area=NUT*Cc/k;
    return Area;
}

/****************************************************************************/

void TExchanger::outletTemp()
{
    double Ch=*Hot.Massflow*Hot.SpecHeat;
    double Cc=*Cold.Massflow*Cold.SpecHeat;

    if(Ch<Cc)
    {
        Hot.StateOut->T=-(Efficiency*(Hot.StateIn->T-Cold.StateIn->T)-Hot.StateIn->T);
        Cold.StateOut->T=Ch/Cc*(Hot.StateIn->T-Hot.StateOut->T)+Cold.StateIn->T;
    }
    else
    {
        Cold.StateOut->T=Efficiency*(Hot.StateIn->T-Cold.StateIn->T)+Cold.StateIn->T;
        Hot.StateOut->T=-(Cc/Ch*(Cold.StateOut->T-Cold.StateIn->T)-Hot.StateIn->T);
    }
}
/****************************************************************************/

void TExchanger::init(){
    hexregister();
    variable  ();
}


TSequence & TExchanger::variable    () {

    /*
    Pin 	= ivbles->Get(TExpander::ga->data->idRecord("Pin"))  ;
    Tin 	= ivbles->Get(TExpander::ga->data->idRecord("Tin"))  ;
    Xin 	= ivbles->Get(TExpander::ga->data->idRecord("Xin"))  ;
    Pout 	= ivbles->Get(TExpander::ga->data->idRecord("Pout")) ;
    Mflow = ivbles->Get(TExpander::ga->data->idRecord("Mflow"));
    Nrot 	= ivbles->Get(TExpander::ga->data->idRecord("Nrot")) ;
    */

    Cold.StateIn->T 	= ivbles->Get(0);
    Cold.StateIn->P 	= ivbles->Get(1);
    Cold.StateIn->x 	= ivbles->Get(2);
    *Cold.Massflow = ivbles->Get(3);
    Hot.StateIn->T 	= ivbles->Get(4);
    Hot.StateIn->P 	= ivbles->Get(5);
    Hot.StateIn->x 	= ivbles->Get(6);
    *Hot.Massflow = ivbles->Get(7);
    Pinch= ivbles->Get(8);
    Area=ivbles->Get(9);
    Q=ivbles->Get(10);
    Cold.StateOut->T 	= ivbles->Get(11);
    Cold.StateOut->P 	= ivbles->Get(12);
    return * ivbles;
}

TSequence & TExchanger::variable    (TSequence& seq) {

    ivbles->copy (seq);
    return  variable();
}
/*
void TExchanger::init_old()
{
    //OUVERTURE ET LECTURE DU FICHIER DONNEES D ENTREES DE L ECHANGEUR
    const char sep = '\t';
    const int nbTitleLines = 1; //nombre de lignes avec des titres
    const int nbTitleColumns = 1; //nombre de colonnes avec des titres ie sans chiffre ie
    const char* file_inDATA = "BasicEnergySystem/exchanger/exchdata.txt";
    TInputFile InputFileDATA((char*)file_inDATA, nbTitleLines, nbTitleColumns, sep);
    InputFileDATA.open();

    TSeqsArray RecordsDATA(InputFileDATA.nRecords());// creation du tableau VIDE avec le même nombre de Lignes que de nombre de records de InputFile
    InputFileDATA.config(Tiofile::onrcds); //on configure le TSeqsArray en mode Records (ie ligne de résultats et pas colonnes)
    RecordsDATA = *InputFileDATA.GetRecords(); //on ajoute les valeurs des Records de l'InputFile

    //On lit les valeurs pour les états d'entrée chaud et froid
    TSequence* PtrSeqTCin = Get(RecordsDATA, "TCin");
    Cold.StateIn->SetT(PtrSeqTCin->Get(4));

    TSequence* PtrSeqPCin = Get(RecordsDATA, "PCin");
    Cold.StateIn->SetP(PtrSeqPCin->Get(4));

    TSequence* PtrSeqxCin = Get(RecordsDATA, "xCin");
    Cold.StateIn->SetX(PtrSeqxCin->Get(4));

    TSequence* PtrSeqmC = Get(RecordsDATA, "mC");
    if (PtrSeqmC->Get(0)==1)//Si la premiere colonne du tableau = 1, on prend la valeur en 4e colonne, sinon non
    {*Cold.Massflow=PtrSeqmC->Get(4);}
    else {*Cold.Massflow=0;}


    TSequence* PtrSeqTHin = Get(RecordsDATA, "THin");
    Hot.StateIn->SetT(PtrSeqTHin->Get(4));

    TSequence* PtrSeqPHin = Get(RecordsDATA, "PHin");
    Hot.StateIn->SetP(PtrSeqPHin->Get(4));

    TSequence* PtrSeqxHin = Get(RecordsDATA, "xHin");
    Hot.StateIn->SetX(PtrSeqxHin->Get(4));

    TSequence* PtrSeqmH = Get(RecordsDATA, "mH");
    if (PtrSeqmH->Get(0)==1)
    {*Hot.Massflow=PtrSeqmH->Get(4);}
    else {*Hot.Massflow=0;}


    TSequence* PtrSeqdTmin = Get(RecordsDATA, "dTmin");
    if (PtrSeqdTmin->Get(0)==1)
    {Pinch=PtrSeqdTmin->Get(4);}
    else {Pinch=0;}

    TSequence* PtrSeqS = Get(RecordsDATA, "S");
    if (PtrSeqS->Get(0)==1)
    {Area=PtrSeqS->Get(4);}
    else {Area=0;}

    TSequence* PtrSeqQ = Get(RecordsDATA, "P");
    if (PtrSeqQ->Get(0)==1)
    {Q=PtrSeqQ->Get(4)*1000;}
    else {Q=0;}

    TSequence* PtrSeqTCout = Get(RecordsDATA, "TCout");
    if (PtrSeqTCout->Get(0)==1)
    {Cold.StateOut->SetT(PtrSeqTCout->Get(4));
        Cold.StateOut->SetT(PtrSeqTCout->Get(4));
        Cold.StateOut->SetP(PtrSeqPCin->Get(4));}
    else {Cold.StateOut->T=0;}

    TSequence* PtrSeqTHout = Get(RecordsDATA, "THout");
    if (PtrSeqTHout->Get(0)==1)
    {Hot.StateOut->SetT(PtrSeqTHout->Get(4));
        Hot.StateOut->SetT(PtrSeqTHout->Get(4));
        Hot.StateOut->SetP(PtrSeqPHin->Get(4));}
    else {Hot.StateOut->T=0;}

    Cold.pStateWall->SetP(PtrSeqPCin->Get(4));
    Hot.pStateWall->SetP(PtrSeqPHin->Get(4));

    double Ta;
    Ta=((Cold.StateIn->GetT()+Cold.StateOut->GetT())/2+(Hot.StateIn->GetT()+Hot.StateOut->GetT())/2)/2;
    Cold.pStateWall->SetT(Ta);
    Hot.pStateWall->SetT(Ta);

    if (!(flw==0 || flw==1 || flw==2))
    {
        cout << "Mode incorrect pour initialisation" << endl;
    }
}
*/



/****************************************************************************/

BasicES & TExchanger::performance() {
    if (flw == 0 || flw == 1 || flw == 2) {

        if (*Hot.Massflow==0 || *Cold.Massflow==0)
        {
            getMassFlow();
        }

        if (!(Hot.StateOut->T == 0 && Cold.StateOut->T == 0) && Area==0) {
            design();
        }
        if (k == 0) {
            kCalculation();
        }

        double Ch;
        double Cc;
        Ch = *Hot.Massflow * Hot.SpecHeat;
        Cc = *Cold.Massflow * Cold.SpecHeat;

        if (Area != 0) {
            //CALCUL DU NUT
            if (Ch < Cc) {
                NUT = Area * k / Ch;
            } else {
                NUT = Area * k / Cc;
            }
            //CALCUL DE L EFFICACITE
            efficiencyByNUT();

            //CALCUL DES TEMPERATURES DE SORTIE
            if(Hot.StateOut->T==0){
                getTho();
            }
            if(Cold.StateOut->T==0){
                getTco();
            }

            //CALCUL DU PINCH
            getPinch();
        }
        else if (Pinch != 0) {
            if (flw == 0) {
                Cold.StateOut->T = (Ch * (Pinch - Hot.StateIn->T) - Cc * Cold.StateIn->T) / (-Ch - Cc);
                Hot.StateOut->T = Cold.StateOut->T + Pinch;
                getSurface();
                efficiency();
            }
            else if (flw == 1) {
                Cold.StateOut->T = Hot.StateIn->T - Pinch;
                Hot.StateOut->T = -(Cc / Ch * (Cold.StateOut->T - Cold.StateIn->T) - Hot.StateIn->T);
                getSurface();
                efficiency();
            }
            else {
                cout << "Mode incorrect pour le pincement: Veuillez rentrer la surface ou changer de type de flux"<< endl;
            }
        }
        else if (Q != 0) {
            Hot.StateOut->T = -(Q / (*Hot.Massflow * Hot.SpecHeat) - Hot.StateIn->T);
            Cold.StateOut->T = Q / (*Cold.Massflow * Cold.SpecHeat) + Cold.StateIn->T;
            getPinch();
            getSurface();
            efficiency();
        }
        else {
            cout << "Veuillez rentrer soit le pincement soit la surface soit la puissance de l'echangeur" << endl;
        }
    }
    else
    {
        cout << "Mode incorrect pour performance" << endl;
    }

/*    Hot.StateOut->SetT(Hot.StateOut->T);
    Hot.tStateOut.SetT(Hot.StateOut->T);
    Cold.StateOut->SetT(Cold.StateOut->T);
    Cold.tStateOut.SetT(Cold.StateOut->T);

    Hot.StateOut->SetP(Hot.StateIn->P);
    Cold.StateOut->SetP(Cold.StateIn->P);*/
}

/****************************************************************************/

float TExchanger::objective(TSequence &seq)
{
   /* if (){
        mxgao mxga ("mxdata.txt", "mxxrslts.txt"); //11-05-2020 le fichier joue le rôle de génome (avec des allèles)
        mxga.minimize(fObjective2, gao::simple);
        cout<<"That's it!"<<endl;
    }
    // A FAIRE suivant le enum
 */
    return 0;
}

/****************************************************************************/

void TExchanger::display()
{
    cout <<"/****************************************************************************/" << endl;
    cout<< endl;
    cout<<"                         Exchanger                           "<<endl;
    cout<<"Pinch Temperature               [C] :       "<<Pinch<<endl;
    cout<<"Heat power                     [kW] :       "<<Q/1E3<<endl;
    cout<<"Total Area                     [m2] :       "<<Area<<endl;
    cout<<"Cost                          [CHF] :       "<<Cost<<endl;
    cout<<"Efficiency                      [-] :       "<<Efficiency<<endl;
    cout<<endl;

    cout<<"                         Hot fluid "<<endl;
    cout<<"Fluid working pressure        [bar] :       "<<Hot.StateIn->P/1E5<<endl;
    cout<<"Inlet Temperature               [C] :       "<<Hot.StateIn->T<<endl;
    cout<<"Outlet Temperature              [C] :       "<<Hot.StateOut->T<<endl;
    cout<<"Mass flow rate               [kg/s] :       "<<*Hot.Massflow<<endl;
    cout<<"Specific heat              [J/kg/K] :       "<<Hot.SpecHeat<<endl;
    cout<<endl;

    cout<<"                         Cold fluid"<<endl;
    cout<<"Fluid working pressure        [bar] :       "<<Cold.StateIn->P/1E5<<endl;
    cout<<"Inlet Temperature               [C] :       "<<Cold.StateIn->T<<endl;
    cout<<"Outlet Temperature              [C] :       "<<Cold.StateOut->T<<endl;
    cout<<"Mass flow rate               [kg/s] :       "<<*Cold.Massflow<<endl;
    cout<<"Specific heat              [J/kg/K] :       "<<Cold.SpecHeat<<endl;
    cout<<endl;
    cout <<"/****************************************************************************/" << endl;
    cout << endl;
}


/****************************************************************************/

void TExchanger::CpCalculation(){
    ////Cp=dH/dT
    Hot.SpecHeat=(Hot.StateIn->GetH()-Hot.StateOut->GetH())/(Hot.StateIn->GetT()-Hot.StateOut->GetT());
    Cold.SpecHeat=(Cold.StateIn->GetH()-Cold.StateOut->GetH())/(Cold.StateIn->GetT()-Cold.StateOut->GetT());
}

/****************************************************************************/

void TExchanger::DefaultParams() {
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    e = 0;
    f = 0;
    k = 0;
    g = 0;
    h = 0;
}

/****************************************************************************/

void TExchanger::RegisterParams(const char *exchname) {
    a = params->Geta(exchname);
    b = params->Getb(exchname);
    c = params->Getc(exchname);
    d = params->Getd(exchname);
    e = params->Gete(exchname);
    f = params->Getf(exchname);
    k = params->Getk(exchname);
    g = params->Getg(exchname);
    h = params->Geth(exchname);
}

/****************************************************************************/


