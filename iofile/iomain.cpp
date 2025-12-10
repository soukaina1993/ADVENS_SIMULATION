#include"iofile.h"
#include"dynstring.h"
#include <iostream>

using namespace std;  // introduces namespace std
/*
int main() {
    // separateur de colonnes
    const char sep = ' ';

    // Fichier d'extraction et d'insertion
    const char* file_in = "iofdata.txt";

    // Nombre de lignes / Colones de titres
    const int nbTitleLines = 1;
    const int nbTitleColumns = 1;

    // Debut du traitement
    cout << "Test Input Ouput file : " << file_in << endl << endl;

    TInputFile  InputFile((char*)file_in, nbTitleLines, nbTitleColumns, sep);
    InputFile.open();

    TOutputFile OutputFile((char*)"iofdata.out");
    OutputFile.open();

    cout << InputFile.firstline() << endl;
    cout << InputFile.firstcolum() << endl;

    //-------------------------------------------------
    int nbFields = InputFile.nFields();
    TSequence RecordSeq(nbFields, TSequence::line, sep);
    int nbLine = InputFile.nLines();
    int nbCol = InputFile.nColums();
    cout << "nombre de lignes  : " << nbLine << endl << endl;
    cout << "nombre de colonnes : " << nbCol << endl << endl;

    // FRO: GetRecord prendd un index 1-based, 2 == 2ieme colone
    // Donc si deux colones de titres, renvoie NULL
    if (InputFile.GetRecord(RecordSeq, nbTitleLines + 1) != NULL)
        cout << RecordSeq << endl;

    if (InputFile.NextRecord(RecordSeq))
        cout << RecordSeq << endl;

    // latchtech: si j'ai bien compris, pas possible de lire "atlas" vu que le tableau s'initialise a la colonne 2
    dynstring turbname = "Atlas";
    if (InputFile.GetRecord(RecordSeq, turbname) != NULL)
        cout << RecordSeq << endl;

    cerr << std::endl;

    //-------------------------------------------------
    int nbRecords = InputFile.nRecords();
    TSequence FieldSeq(nbRecords, TSequence::colum, sep);

    // latchtech: que represente l'index 2 ? (0-based ou 1-based)
    if (InputFile.GetField(FieldSeq, nbTitleLines + 1) != NULL)
    {
        cout << FieldSeq << endl;
    }

    if (InputFile.NextField(FieldSeq, nbTitleLines + 1))
    {
        cout << FieldSeq << endl;
    }

    if (InputFile.NextField(FieldSeq))
    {
        cout << FieldSeq << endl;
    }

    cout << InputFile.titleline() << endl;

    //-------------------------------------------------
    nbRecords = InputFile.nRecords();
    TSeqsArray Records(InputFile.nRecords());
    InputFile.config(Tiofile::onrcds);
    //InputFile>>Records;
    Records = *InputFile.GetRecords();

    //-------------------------------------------------
    nbFields = InputFile.nFields();
    TSeqsArray Fields(InputFile.nFields());
    InputFile.config(Tiofile::onflds);
    //InputFile>>Fields;
    Fields = *InputFile.GetFields();

    //dynstring Titre = InputFile.titleline();
    //cout<<InputFile.titleline(Titre)<<endl;

    //-----------------------------------------------------------

    OutputFile.SetFileName("iofrcd1.out");
    OutputFile << RecordSeq << Records;
    OutputFile.SetFileName("iofrcd2.out");
    OutputFile.titleline(InputFile.titleline());
    OutputFile << Records << Records;

    //-----------------------------------------------------------

    OutputFile.SetFileName("ioffld1.out");
    OutputFile << Fields << Fields;

    //-----------------------------------------------------------

    OutputFile.SetFileName("ioffld2.out");
    if (!OutputFile.Set2DArray(Records, Fields))
        cerr << "SetFields (); File Error : " << OutputFile.GetFileName() << endl << endl;

    InputFile.close();
    OutputFile.close();

    //-----------------------------------------------------------

    cout << endl << endl;
    //dynstring linename="Maneur", colname="Diam";

    // latchtech: si on suit le fil : le "Records" lu ici ne contient pas de titre vu qu'on commence la lecture a row=2, col=2
    char* temp = (char*)"Maneur";
    TSequence* PtrSeq = Get(Records, temp);
    if (!(PtrSeq == NULL))
        cout << *PtrSeq;
    PtrSeq = Get(Fields, (char*)"Diam");
    if (!(PtrSeq == NULL))
        cout << *PtrSeq;

    cout << endl << endl;
    unsigned int SeqID = 3;
    PtrSeq = Get(Records, SeqID);
    if (!(PtrSeq == NULL))
        cout << *PtrSeq;
    PtrSeq = Get(Fields, SeqID);
    if (!(PtrSeq == NULL))
        cout << *PtrSeq;

    if (PtrSeq != NULL) {
        cout << PtrSeq->Get(0) << endl;
        cout << PtrSeq->Set(0, 70) << endl;
        cout << *PtrSeq << endl;
    }


    if (isexist(Records, (char*) "Maneur"))
        cout << "oui" << endl;
    else
        cout << "non" << endl;
    unsigned int n1 = howoften(Records, (char*) "Maneur");
    cout << n1 << endl;

    return 0;
}

*/