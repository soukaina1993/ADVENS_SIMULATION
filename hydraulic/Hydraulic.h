#ifndef DEF_HYDRAULIQUE
#define DEF_HYDRAULIQUE

#include "../BasicEnergySystem/BasicES.h"
/*
===============================================================================================================
DESCRIPTION
Dans ce code on int�gre l'ensemble des fonctions que l'on va utiliser pour faire les calculs sur les pertes de
charges, pertes thermiques, ... pour l'ensemble des composants hydraulique du r�seau

PDC = Perte de charge
PT = Perte thermique
lin = lin�aire
K = coef de perte de charge singuliere
sing = singuli�re
T� = Coude � 90 degr�s
===============================================================================================================
*/

class Hydraulic:public BasicES// HydrauliC dérive de BASICES
{
    public:
    Hydraulic(void);
// Calcul de la vitesse
    double Calcul_Vitesse (double Qv,
                           double diametre);
// Calcul du nombre de Reynolds
    double Calcul_Reynolds(double Ro,
                           double vitesse,
                           double diametre,
                           double ViscDyn);

// Calcul du coefficient de perte de charge
    double Calcul_MoodyFactor(double diametre,
                              double RugAbs,
                              double Reynolds);




// dans une jonction � s�paration sym�trique




/*
================================================================================================================
On integre les fonctions de dimensionnement
================================================================================================================
*/


// Calcul du diam�tre en fonction des pertes de charge
    double Calcul_diametre_int( double DeltaPr,
                                double ViscDyn,
                                double Ro,
                                double Qv);

};
#endif

// YOLAINE: Les méthodes enlevées
// Calcul des pertes de charge lin�aires --> maintenant dans Pipe
/*double Calcul_PDC_lin(  double longueur,
                        double diametre,
                        double Qv,
                        double MoodyFactor);

Calcul des pertes thermiques lin�aires --> dans Pipe
   double Calcul_PT_lin (  double &T_out,
                            double &Perte_thermique,
                            double longueur,
                            double T_in,
                            double T_amb,
                            double diametre_int,
                            double diametre_ext,
                            double Qv,
                            double Ro,
                            double cp);

// Calcul du coef de perte de charges singuliere :
// d'un coude brusque
    double Calcul_K_coude_brusque (double alpha);

// d'un coude arrondi
    double Calcul_K_coude_arrondi (double alpha,
                                        double diametre,
                                        double rayon_courbure);

 double Calcul_K_divergent (double Qin,
                               double Qout);

 double Calcul_K_convergent(double Qin,
                               double Qout);

                        */
