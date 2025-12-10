#include "Hydraulic.h"
#include <cmath>
#include <iostream>

using namespace std;

const double Pi = 3.14159;
const double g = 9.81;

// CONSTRUCTEUR
Hydraulic::Hydraulic(void)
{

}
//=================================================================================================
// ECOULEMENT
//=================================================================================================
double
Hydraulic::Calcul_Vitesse(double Qv, double diametre)
{
    double vitesse;
    vitesse = 4*Qv/(Pi*diametre*diametre);
    return vitesse;
}
//=================================================================================================
double
Hydraulic::Calcul_Reynolds(double Ro, double vitesse, double diametre, double ViscDyn)
{
    double Reynolds;
    Reynolds = Ro*vitesse*diametre/ViscDyn;
    return Reynolds;
}
//=================================================================================================
// COEF DE PERTE DE CHARGE REGULIERE
//=================================================================================================
double
Hydraulic::Calcul_MoodyFactor(double diametre,
                              double RugAbs,
                              double Reynolds)
{
    double MoodyFactor, RugRel;
    if (Reynolds < 2000)
    {
        MoodyFactor = 64 / Reynolds;
    }

    else
    {
        double var = 10;
        double precvar = 0.01;
        double errvar = 0.1;
        double newvar;

        RugRel = RugAbs/diametre;

        while (errvar > precvar)
        {
            //exp(1.03*log(RugRel/3.71))+exp(1.1 * log(4.26/Reynolds))*exp(1.1*log(var));
            newvar = (-1.94 * log(exp(1.03*log(RugRel/3.71))+ exp(1.1 * log(4.26/Reynolds))*exp(1.1*log(var)))/log(10));
            errvar = fabs(newvar - var);
            var = newvar;
        }
    MoodyFactor = 1/sqrt(newvar);
    }
    cout << "le coef de perte de charge (Moody) est de : " << MoodyFactor << endl;

    return MoodyFactor;
}
//=================================================================================================
// Revoir la formule

//=================================================================================================
// COEF DE PERTE DE CHARGE SINGULIERE
//=================================================================================================

//=================================================================================================


//=================================================================================================
// FONCTIONS THERMIQUES
//=================================================================================================

//=================================================================================================
// OPTIMISATION
//=================================================================================================
double
Hydraulic::Calcul_diametre_int(double DeltaPr,
                               double ViscDyn,
                               double Ro,
                               double Qv)
{
    double diametre = 1000;
    double precdiam = 0;
    double newdiam;
    double precision = 0.01;
    double vitesse, Reynolds, MoodyFactor, RugAbs, longueur;

        while (fabs(precdiam - diametre) > precision)
            {
                precdiam = diametre;
                vitesse = 4*Qv/(Pi*diametre*diametre);
                Reynolds = Ro*vitesse*diametre/ViscDyn;

                MoodyFactor = Calcul_MoodyFactor(diametre, RugAbs, Reynolds);
                newdiam = exp(log(8*longueur*MoodyFactor*(Qv/Pi)*(Qv/Pi)/(DeltaPr *Ro))/5);
                diametre = newdiam;
            }

        cout << "Le diametre correspondant est " << diametre << " metres " << endl;

    return diametre;
}
//================================= FIN ======================================



