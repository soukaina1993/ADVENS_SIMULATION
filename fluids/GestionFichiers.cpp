#include <iostream>
#include <fstream>

using namespace std;

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
                double &c4)
//il faut ajouter la fonction qui choisi le fluide adapt�
{
    ifstream ChercherGaz ("fluids/banque.txt");
    if (ChercherGaz)
    {
        string nom1;

        cout << "quel fluide voulez-vous utliser ?" << endl;
        cin >> nom1;
        cout << "vous desirez utiliser le :" << nom1 << endl;
        while (nom1.compare(nom) != 0) //si c'est egal a 0, ca veut dire que les deux string sont egaux //tant que les string ne sont pas egaux, tu fais la boucle
        {
            ChercherGaz >> nom;
            ChercherGaz >> vc;
            ChercherGaz >> Pc;
            ChercherGaz >> Tc;
            ChercherGaz >> Tnb;
            ChercherGaz >> Mw;
            ChercherGaz >> w;
            ChercherGaz >> DipM;
            ChercherGaz >> c0;
            ChercherGaz >> c1;
            ChercherGaz >> c2;
            ChercherGaz >> c3;
            ChercherGaz >> c4;
            //PROBLEME ICI ! si tu ne trouves pas le fluide, boucle infinie!!
        }
    }
     else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }

    return;
}

void AjouterFluide();
void AjouterFluide()

{
    string nom;
    double vc,Pc,Tc,Tnb,Mw,w,DipM,c0,c1,c2,c3,c4;

    ofstream AjouterGaz ("fluids/banque.txt",ios::app);
    if (AjouterGaz)
    {
        cout << "quel est son nom ? " << endl;
        cin  >> nom;
        cout << "quel est son vc ? " << endl;
        cin  >> vc;
        cout << "quel est son Pc ? " << endl;
        cin  >> Pc;
        cout << "quel est son Tc ? " << endl;
        cin  >> Tc;
        cout << "quel est son Tnb ? " << endl;
        cin  >> Tnb;
        cout << "quel est son Mw ? " << endl;
        cin  >> Mw;
        cout << "quel est son w ? " << endl;
        cin  >> w;
        cout << "quel est son DipM ? " << endl;
        cin  >> DipM;
        AjouterGaz << nom << " " << vc << " " << Pc << " " << Tc << " " << Tnb << " " << Mw << " " << w << " " << DipM << " " ;

        cout << "c0 ?, c1 ?, c2 ?, c3 ?, c4 ?" << endl;
        cin >> c0;
        cin >> c1;
        cin >> c2;
        cin >> c3;
        cin >> c4;
        AjouterGaz << c0 << " " << c1 << " " << c2 << " " << c3 << " " << c4 << endl;

    }
     else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }
return;
}

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
                     double Prandtl)
{
    ofstream AjouterResultat ("Resultats.txt",ios::app);
    if (AjouterResultat)
    {
        AjouterResultat << v << " " << P << " " << T << " " << u << " " << h << " " << s << " " << cp << " " << cv << " " ;
        AjouterResultat << ViscDyn << " " << ThermCond << " " << Ro << " " << Prandtl << endl;
    }
     else
    {
        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
    }
return;
}
