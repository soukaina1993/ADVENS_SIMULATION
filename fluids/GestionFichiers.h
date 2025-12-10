#include <iostream>
#include <ctype.h>#include <math.h>
#include <tgmath.h>#include <stdlib.h>#include <string.h>

void LireFluide(string &nom,
                double &vc,
                double &Pc,
                double &Tc,
                double &Tnb,
                double &Mw,
                double &w,
                double &DipM,
                double &c0,
                double &c1,
                double &c2,
                double &c3,
                double &c4);

void AjouterFluide();

void EcrireResultat( string nom,
                     double v,
                     double P,
                     double T,
                     double u,
                     double h,
                     double s,
                     double cp,
                     double cv,
                     double ViscDyn,
                     double ThermCond,
                     double Ro,
                     double Prandtl);


