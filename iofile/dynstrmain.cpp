#ifdef _MSC_VER
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include "dynstring.h"
/*
int main()
{
	std::cout << "Hello World!\n";

	dynstring Texte = "Instance construite par texte = \"...\"";

	dynstring Chaine(Texte.length() + 4, '-', 50);
	cout << "cout << Texte << endl << Chaine << endl -->" << endl
		 << Texte << endl << Chaine << endl << endl;

	cout << "Saisie de chaine par cin >> chaine" << endl << endl;
	cin >> Chaine;

	cout << Chaine << endl << endl;


	Texte[0] = 'D';
	Texte[1] = 'i';

	cout << Texte << endl << endl;

	Texte.empty();
	Texte.copy(Chaine);
	cout << Texte << endl << endl;

	Texte = "je vais";
	cout << Texte << endl << endl;


	Chaine.empty();
	Chaine = "à l'école";
	cout << Chaine << endl << endl;

	dynstring * Pstr = Chaine.clone();
	cout <<*Pstr << endl << endl;
	cout <<Pstr->length() << endl << endl;
	cout <<Pstr->size() << endl << endl;


	char liste[] = " à l'école";
	dynstring Str = Texte + liste;
	cout <<Str << endl << endl;

	dynstring Str2 = Texte + liste + " demain ";
	cout <<Str2 << endl << endl;
	cout <<Str2.length() << endl << endl;


	char car = 'm';
	cout << Str2+car << endl << endl;


	dynstring Str3 = "mais sans stress ";
	dynstring Str4 = Str2+Str3;
	cout <<Str4 << endl << endl;
	cout <<Str4.length() << endl << endl;
	cout <<Str4.size() << endl << endl;


	dynstring Str5 = Str4+Chaine;
	cout <<Str5 << endl << endl;
	cout <<Str5.length() << endl << endl;

	Str4 = "C'est mon choix";
	Str5 = "C'est mon choix";

	if (Str4 != Str5)
	    cout << "les deux string sont différents" << endl;
	else
	    cout << "les deux string sont identiques" << endl;

	cout << Str4.length() << endl;
	cout << Str4.size() << endl;

	cout << Str5.length() << endl;
	cout << Str5.size() << endl;

	return 0;
}
*/