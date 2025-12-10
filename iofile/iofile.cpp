/****************************************************************************/
#include <fstream>
#include <iostream>
#include <string>
#include "array.h"
#include "dynstring.h"
#include "iofile.h"
#include <sstream>

#include <stdlib.h>//ajout le 29-03-2020 pour exit
#include <math.h>

/****************************************************************************/
using namespace std;
//-------------------------------------------------
TSequence * Get (const TSeqsArray &Seqs, char * ch, unsigned int cptr)
{
   if (cptr >= Seqs.size())
      cptr=0;
   unsigned int cpt=cptr;
   dynstring st = ch;

	   while (st != *Seqs[cpt]->Title && cpt != Seqs.size() - 1) {
		cpt++;
		}

   if (st!=*Seqs[cpt]->Title)
      return NULL;
   return Seqs[cpt];
}

TSequence * Get (const TSeqsArray &Seqs, int d,char * ch, unsigned int cptr)
{
    if (cptr >= Seqs.size())
        cptr=0;
    unsigned int cpt=cptr;
    dynstring st = ch;

    while (st != *Seqs[cpt]->Title && cpt != Seqs.size() - 1) {
        cpt++;
    }

    if (st!=*Seqs[cpt]->Title)
        return NULL;
    return Seqs[cpt+d];
}
//-------------------------------------------------
TSequence * Get (const TSeqsArray &Seqs, const unsigned int Index)
{
    int a= Index-1;
    if (!(a < Seqs.size()+1))
      return NULL;
   return Seqs[a];
}

double GetNum (const TSeqsArray &Seqs, char * ch, unsigned int cptr)
{
    if (cptr >= Seqs.size())
        cptr=0;
    unsigned int cpt=cptr;
    dynstring st = ch;

    while (st != *Seqs[cpt]->Title && cpt != Seqs.size() - 1) {
        cpt++;
    }

    if (st!=*Seqs[cpt]->Title)
        return 0;
    return cpt;
}


//-------------------------------------------------
bool isexist (const TSeqsArray &Seqs, char * sname)
{
   unsigned int cpt=0;
   dynstring st=sname;
   while (st!=*Seqs[cpt]->Title && cpt!=Seqs.size()-1)
      cpt++;
   if (st!=*Seqs[cpt]->Title)
      return false;
   return true;
}
//-------------------------------------------------
unsigned int howoften (const TSeqsArray &Seqs, char * sname)
{
   unsigned int cpt=0;
   dynstring st=sname;

   for (unsigned int i=0; i < Seqs.size(); i++){
      if (st == *Seqs[i]->Title)
         cpt++;
   }
   return cpt;
}
/****************************************************************************/
/*                            Input file class                              */
/****************************************************************************/

//-------------------------------------------------
TInputFile::TInputFile (char * AFileName, unsigned int nTline, unsigned int nTcol, const char colDelimiter) : Tiofile(AFileName)
{
  //FFileName = AFileName;
  config(Tiofile::onrcds);
  nTleLn=nTline;
  nTleCo=nTcol;
  delimiter = colDelimiter;

  if (Dimension ()){
     FRecord = new TSequence (nFlds, TSequence::line, colDelimiter);
     FField  = new TSequence (nRcds, TSequence::colum, colDelimiter);
     FRcds = new Array<TSequence*> (nRecords());
     if (GetRecords (*FRcds)==NULL){
       cerr <<"GetRecords (); File Error : " << GetFileName() << endl;
       FRcds = NULL; exit(1);
     }
     FFlds = new Array<TSequence*> (nFields());
     if (GetFields (*FFlds)==NULL){
       cerr <<"GetRecords (); File Error : " << GetFileName() << endl;
       FFlds = NULL; exit(1);
     }
  }
  FFileID = new ifstream (FFileName);
}

    //-------------------------------------------------
    TInputFile::TInputFile (TInputFile & inFile): Tiofile(inFile.GetFileName()){
          *this=inFile;
    }


//-------------------------------------------------
bool
TInputFile::Dimension () {
      if (!((nLine = GetLineNumber())==0)){
         nRcds= nLine-nTleLn;
      }
      else return false;
      if (!((nCol = GetColumNumber ())==0)){
         nFlds= nCol-nTleCo;
      }
      else return false;
      return true;
   }

 //-------------------------------------------------
 //ligne de titre obtunue a partir des s�quences
 dynstring & TInputFile::titleline  (const char* Name) {
	 int nbFields = nFields();
      Array <TSequence*> temp(nbFields);
      temp=*FFlds;
      Tls->empty();
      *Tls = (char*)Name;
      for (unsigned int i=0; i<nFields(); i++){
         *Tls = *Tls + "	";
         *Tls = *Tls + *temp[i]->Title ;
      }
      return *Tls;
   }
  //-------------------------------------------------
  //ligne de titre directement lue au fichier
 dynstring & TInputFile::firstline  () {
      char Buffer[STREAM_LINE_MAX_SIZE];
      ifstream FileID (FFileName, ios::in);
      if (FileID.fail()) {
          cerr <<"Erreur � l'ouverture de " << FFileName << endl << endl;
          return *Tls="";
      }
      else{
          Tls->empty();
          FileID.getline(Buffer, sizeof(Buffer));
          *Tls=Buffer;
      }
      FileID.close();
      return *Tls;
   }
 //-------------------------------------------------
  //titre en colonne directement lue au fichier
 dynstring & TInputFile::firstcolum  () {
      dynstring ch;
      ifstream FileID (FFileName, ios::in);
      if (FileID.fail()) {
          cerr <<"Erreur � l'ouverture de " << FFileName << endl << endl;
          return *Tls="";
      }
      else{
          Tls->empty();
          FileID>>*Tls;
          for (unsigned int i=0;i<nLines()-1;i++){
              ignoreLine(&FileID, STREAM_LINE_MAX_SIZE, LINE_DEFAULT_DELIMITER);
              FileID>>ch;
              *Tls = *Tls  + "	"+ ch;
          }
      }
      FileID.close();
      return *Tls;
   }
//-------------------------------------------------
 void TInputFile::debug(const char* tag) {
	 // dbg latch
	// get current pos

	 int pos = (int)FFileID->tellg();
	 char test[256] = { '\0' };
	 char* reader = &test[0];
	 char c = '\0';
	 while (!FFileID->eof() && c != ' ') {
		 c = '\0';
		 FFileID->get(c);
		 *(reader++) = c;
	 }
	 *reader = '\0';
	 FFileID->seekg(pos, ios::beg);
	 cerr << tag << " - reader: " << pos << ", " << test << std::endl;
	 // end dbg latch
}

void
TInputFile::InitFilePtr(unsigned int i, unsigned int j)
   {
    // cout << "iofile.cpp\n";
      if (i<=nLine && j<=nCol) {
         FFileID->seekg(0, ios::beg);
         FRecord->noID=0;
         FField->noID=0;
		 DownWord(i);
		 NextWord(j);

      }else{
         cerr <<"Index de ligne ou de colonne incorrect"<< endl << endl;
         exit(1);
      }
   }
//-------------------------------------------------
unsigned int
TInputFile::GetLineNumber (){
    char Buffer[STREAM_LINE_MAX_SIZE];
    int nLine=0;
    ifstream FileID (FFileName, ios::in);
	FileID.clear();
	FileID.seekg(0, ios::beg);
    while (!FileID.eof()) {
       if (FileID.fail()) {
          cerr <<"GetLineNumber Erreur a l'ouverture de " << FFileName << endl << endl;
          return -1;
       }
       else if(FileID.getline(Buffer, sizeof(Buffer)))
           nLine++;
    }
    FileID.close();
    return nLine;
}

//-------------------------------------------------
unsigned int
TInputFile::GetWordNumber (){
    char Buffer[STREAM_COL_MAX_SIZE];
    int nWord=0;
    ifstream FileID (FFileName, ios::in);
	FileID.clear();
	FileID.seekg(0, ios::beg);
    while (!FileID.eof()) {
       if (FileID.fail()) {
          cerr <<"GetWordNumber Erreur a l'ouverture de " << FFileName << endl << endl;
          return -1;
       }
       else{
          FileID >> Buffer;
          nWord++;
       }
    }
    FileID.close();
    return nWord;
}

//-------------------------------------------------
unsigned int
TInputFile::GetColumNumber (){
    int nColum;
    nColum =  (int) (GetWordNumber() / GetLineNumber()) ; // (int / int) = int
    return nColum;
}

//-----------------------------------------------------------
//Ignore un word et incr�mente le compteur de colonne
//-----------------------------------------------------------
bool
TInputFile::NextWord (){
    if (FFileID->fail()){
		cerr <<"Erreur de traitement du fichier " << strerror(errno) << endl << endl;

		return false;
		exit(1);
     }else{
        ignoreLine(FFileID, STREAM_COL_MAX_SIZE, this->delimiter);
        FField->noID++;
     }
    return true;
}

//-----------------------------------------------------------
//Ignore un nombre donn� de word et incr�mente autant le compteur
//de colonne
//-----------------------------------------------------------
bool
TInputFile::NextWord (unsigned int Number){
	for (unsigned int i = 0; i < Number; i++) {
		if (! NextWord()) {
			return false;
		}
	}

    return true;

}
//-----------------------------------------------------------
//Ignore un nombre donn� de lignes en utilisant la m�thode normale ou le hack
// en fonction du #define USE_IGNORE_LINE_HACK
// delimiter attendu = '\n'
//-----------------------------------------------------------
void ignoreLine(ifstream* stream, unsigned int size, const char delimiter) {
#if USE_IGNORE_LINE_HACK
	ignoreLineHack(stream, size, delimiter);
#else
	stream->ignore(size, delimiter);
#endif
}
//-----------------------------------------------------------
//Ignore un nombre donn� de lignes en utilisant � l'aide d'une boucle while
// delimiter attendu = '\n'
//-----------------------------------------------------------
void ignoreLineHack(ifstream* stream, unsigned int size, const char delimiter) {
	char c = '\0';
	for (unsigned int i = 0; i < size; ++i) {
		stream->get(c);
		if (c == delimiter || stream->eof()) {
			break;
		}
	}
}
//-----------------------------------------------------------
bool
TInputFile::DownWord (){
    if (FFileID->fail()){
		cerr <<"Erreur de traitement du fichier " << strerror(errno) << endl << endl;
		return false;
		exit(1);
     }else{
        if (FRecord->noID<nLine-1){ //nRcds
		   ignoreLine(FFileID, STREAM_LINE_MAX_SIZE, LINE_DEFAULT_DELIMITER);
           FRecord->noID++;
        }else{
           FFileID->seekg(ios::beg);
           FRecord->noID = nTleLn-1;
           FField->noID++;
        }
        for (unsigned int i=0;i<FField->noID;i++)
              ignoreLine(FFileID, STREAM_COL_MAX_SIZE, this->delimiter);
     }
    return true;
}
//-----------------------------------------------------------
bool
TInputFile::DownWord (unsigned int Number){
    if (FFileID->fail()){
		cerr <<"Erreur de traitement de fichier " << endl << endl;
		return false;

    } else {
        if (!(Number==0)){
           for (unsigned int i=0;i<Number;i++){
              DownWord ();
           }
        }
    }
    return true;
}

//-----------------------------------------------------------
bool
TInputFile::NextRecord (TSequence &Sequence){
    if (FFileID->fail()){
	cerr <<"Erreur de traitement du fichier " << strerror(errno) << endl << endl;
	return false;
	exit(1);
     }else{
        InitFilePtr(FRecord->noID+1, nTleCo-1);
        GetRecord(Sequence);
     }
    return true;
}

//-----------------------------------------------------------
bool
TInputFile::NextRecord (TSequence &Sequence, unsigned int Number){
    if (FFileID->fail()){
	cerr <<"Erreur de traitement de fichier "  << endl << endl;
	return false;
	exit(1);
     }else{
        InitFilePtr(FRecord->noID+Number, nTleCo-1);
        GetRecord(Sequence);
    }
    return true;
}
//-----------------------------------------------------------
bool
 TInputFile::NextField (TSequence &Sequence){
    if (FFileID->fail()){
	cerr <<"Erreur de traitement du fichier " << strerror(errno) << endl << endl;
       // cout<<Sequence<<endl;
        cout<<Sequence.size()<<endl;

	return false;
	exit(1);
     }else{
        InitFilePtr(nTleLn-1, FField->noID+1);
        GetField(Sequence);
     }
    return true;
}

//-----------------------------------------------------------
bool
TInputFile::NextField (TSequence &Sequence, unsigned int Number){
    if (FFileID->fail()){
	cerr <<"Erreur de traitement de fichier "  << endl << endl;
	return false;
	exit(1);
     }else{
        InitFilePtr(nTleLn-1, FField->noID+Number);
        GetField(Sequence);
    }
    return true;
}
//-----------------------------------------------------------
//Lecture de la ligne (ou Record) courante et remise a jour des compteurs.
//Retourne un pointeur sur la structure ou NULL si probl�me
//-----------------------------------------------------------
TSequence * TInputFile::GetRecord (TSequence &Sequence){
     if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return NULL;
	exit(1);
     }
     else{
        InitFilePtr(FRecord->noID, nTleCo-1);
        Sequence.copy(*FRecord);
        *FFileID >> Sequence;
     }
  return &Sequence;
}
//-----------------------------------------------------------
TSequence * TInputFile::GetRecord (TSequence &Sequence, const dynstring &Str){
   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return NULL;
	exit(1);
   }
   else{
	InitCursor();
	do
	{
	   DownWord();
	   *FFileID>>*Sequence.Title;
	   InitFilePtr(FRecord->noID, FField->noID);
	} while (Str!=*Sequence.Title && FRecord->noID!=nLine-1);
	if (Str!=*Sequence.Title)
            return NULL;
   }
   return GetRecord (Sequence);
 }
//-----------------------------------------------------------
//Lecture du record Index� a partir de la valeur courante.
//Retourne un pointeur sur la structure ou NULL si probl�me
//-----------------------------------------------------------
TSequence * TInputFile::GetRecord (TSequence &Sequence, unsigned int Index){
  if (Index>=nTleLn && Index<GetLineNumber()) {
     if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return NULL;
	exit(1);
     }
     else{
	InitFilePtr(Index-1, nTleCo-1);
	Sequence.copy(*FRecord);
	*FFileID >> Sequence;
     }
  }else{
     cerr <<"L'index donn� correspond � une ligne de titre"<< endl << endl;
     cerr <<"ou est sup�rieur au nombre total de ligne"<< endl << endl;
     return NULL;
     exit(2);
  }
  return &Sequence;
}
//-----------------------------------------------------------
Array <TSequence*> * TInputFile::GetRecords (Array <TSequence*> &ARcds){
   TSequence TempSeq (FRecord->nElements, TSequence::line, this->delimiter);
   ifstream * OldFileID = FFileID;
   unsigned int nbRcds = nRcds;
   ifstream FileID (FFileName, ios::in);
   FFileID = &FileID;
   if (!(ARcds.size() <= nRcds)) nbRcds = nRcds;
   else nbRcds = ARcds.size();

   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return NULL;
	exit(1);
    }
    else{
	InitFilePtr(nTleLn-1, nTleCo-1);
	for (unsigned int i= 0; i<nbRcds ; i++){
		if (NextRecord(TempSeq)) {
			ARcds[i] = TempSeq.clone();
		}
		else {
			std::cerr << "Erreur d�tection d'un item : record = " << i << std::endl;
		}

	}
    }
   FFileID = OldFileID;
   FileID.close();
   return &ARcds;
 }
//-----------------------------------------------------------
//Lecture de la colonne (ou Field) courante et mis a jour des compteurs.
//Retourne un pointeur sur la structure ou NULL si probl�me
//-----------------------------------------------------------
TSequence * TInputFile::GetField (TSequence &Sequence){
     if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return NULL;
	exit(1);
     }
     else{
	InitFilePtr(nTleLn-1, FField->noID);
	Sequence =*FField;
	*FFileID >> Sequence;
     }
  return &Sequence;
}
//-----------------------------------------------------------
//Lecture du field Index� a partir de la valeur courante.
//Retourne un pointeur sur la structure ou NULL si probl�me
//-----------------------------------------------------------
TSequence * TInputFile::GetField (TSequence &Sequence, unsigned int Index){
	unsigned int colNumber = GetColumNumber();
	if (Index >= nTleCo && Index < colNumber) {
		if (FFileID->fail()){
			cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
			return NULL;
			exit(1);
		}
		else{
			InitFilePtr(nTleLn-1, Index-1);
			Sequence =*FField;
			*FFileID >> Sequence;
		}
	}else{
		cerr <<"L'index donn� correspond � une ligne de titre"<< endl << endl;
		cerr <<"ou est sup�rieur au nombre total de ligne"<< endl << endl;
		return NULL;
		exit(2);
	}
	return &Sequence;
}
//-----------------------------------------------------------
Array <TSequence*> * TInputFile::GetFields (Array <TSequence*> &AFlds){
   TSequence TempSeq (FField->nElements, TSequence::colum, this->delimiter);
   ifstream * OldFileID = FFileID;
   unsigned int nbFlds = nFlds;
   ifstream FileID (FFileName, ios::in);
   FFileID = &FileID;
   if (!(AFlds.size() <= nFlds)){
      nbFlds = nFlds;
   }else{
      nbFlds = AFlds.size();
   }
   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return NULL;
	exit(1);
   }
   else{
	InitFilePtr(nTleLn-1, nTleCo-1);
	for (unsigned int i= 0; i<nbFlds ; i++){
	   NextField(TempSeq);
	   AFlds[i] = TempSeq.clone();
	}
   }
   FFileID = OldFileID;
   FileID.close();
   return &AFlds;
 }

/****************************************************************************/
/*                             output file class                            */
/****************************************************************************/


//-----------------------------------------------------------
//constructeur d'un identificateur de fichier de sortie
//-----------------------------------------------------------
TOutputFile::TOutputFile (char * AFileName): Tiofile(AFileName)
{
      //FFileName = AFileName;
      config(Tiofile::onrcds);
      FFileID = new ofstream (FFileName, ios::out);
      TOutputFile::close();
}
//-----------------------------------------------------------
//Ecrit la s�quence indiqu�e dans le fichier de sortie
//Ouverture de fichier en mode ecriture simple
//-----------------------------------------------------------
bool TOutputFile::SetRecord (TSequence &Sequence){
   if (!isopened())
      TOutputFile::open();
   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return false;
	exit(1);
     }
     else{
	*FFileID<<Sequence;
     }
   return true;
 }

//-----------------------------------------------------------
//Additionne la s�quence indiqu�e dans le fichier de sortie
//Ouverture de fichier en mode ajout
//-----------------------------------------------------------
bool TOutputFile::AddRecord (TSequence &Sequence){
   if (!isopened())
      {FFileID->open(FFileName, ios::app);isopen = TOutputFile::opened;}
   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return false;
	exit(1);
     }
     else{
	*FFileID<<Sequence;
     }
   return true;
 }

//-----------------------------------------------------------
//Ecrit les s�quences en ligne dans le fichier de sortie
//Ouverture de fichier en mode �criture simple
//-----------------------------------------------------------
bool TOutputFile::SetRecords (Array<TSequence*> &Rcds){
   if (!isopened())
      {FFileID->open(FFileName, ios::out);isopen = TOutputFile::opened;}
   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return false;
	exit(1);
     }
     else{
	for (unsigned int i=0; i<Rcds.size(); i++)
	   *FFileID<<*Rcds[i];
     }
   return true;
 }
//-----------------------------------------------------------
//Ecrit les s�quences en ligne dans le fichier de sortie
//Ouverture de fichier en mode ajout
//-----------------------------------------------------------
bool TOutputFile::AddRecords (Array<TSequence*> &Rcds){
   if (!isopened())
      {FFileID->open(FFileName, ios::app);isopen = TOutputFile::opened;}
   if (!SetRecords(Rcds))
       return false;
   return true;
 }
//-----------------------------------------------------------
//�crit les s�quences en colonne dans le fichier de sortie
//Ouverture de fichier en mode �cciture simple
//-----------------------------------------------------------
bool TOutputFile::SetFields (Array<TSequence*> &Flds){
   if (!isopened()) {
      FFileID->open(FFileName, ios::out);
      isopen = TOutputFile::opened;
      //TOutputFile::titleline (Flds, "	");
      TOutputFile::titleline (Flds);
   }
   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return false;
	exit(1);
   }
   else{
	writevalues(Flds);
   }
   return true;
 }
//-----------------------------------------------------------
//Ecrit les s�quences en colonne dans le fichier de sortie
//Ouverture de fichier en mode ajout
//-----------------------------------------------------------
bool TOutputFile::AddFields (Array<TSequence*> &Flds){
   if (!isopened())
      {FFileID->open(FFileName, ios::app);isopen = TOutputFile::opened;}
   if (FFileID->fail()){
	cerr <<"Erreur de manipulation du fichier" << FFileName << endl << endl;
	return false;
	exit(1);
   }
   else{
	writevalues(Flds);
   }
   return true;
}

//-----------------------------------------------------------
//Ecrit la ligne de titre  dans le fichier de sortie
//Ouverture de fichier en mode �criture
//-----------------------------------------------------------
void TOutputFile::titleline (const Array<TSequence*> &Flds, const char* Name){
   if (!isopened())
      {FFileID->open(FFileName, ios::out);isopen = TOutputFile::opened;}
   if (strcmp(Name, "	") != 0)
      *FFileID<<Name;
   for (unsigned int i=0; i<Flds.size(); i++)
      *FFileID<<'	'<<*Flds[i]->Title;
   *FFileID<<endl;
}
//-----------------------------------------------------------
void TOutputFile::titleline (dynstring & tleline){
   if (!isopened())
      {FFileID->open(FFileName, ios::out);isopen = TOutputFile::opened;}
   *FFileID<<tleline<<endl;
}

//-----------------------------------------------------------
//Ecrit les lignes de donn�es  dans le fichier de sortie
//Ouverture de fichier en mode �criture
//-----------------------------------------------------------
void TOutputFile::writevalues  (Array<TSequence*> &Flds){
   if (!isopened())
      {FFileID->open(FFileName, ios::app);isopen = TOutputFile::opened;}
   for (unsigned int j=0; j<Flds[0U]->nElements; j++){
      for (unsigned int i=0; i<Flds.size(); i++){
          *FFileID<<"\t"<<Flds[i]->Elements->Get(j);
      }
      *FFileID<<endl;
   }
}
//-----------------------------------------------------------
//Reconstitue l'ens. du tableau � partir des s�quences lignes
//et colonnes. Ouverture du fichier en mode �criture
//-----------------------------------------------------------
bool TOutputFile::Set2DArray (Array<TSequence*> &Rcds, Array<TSequence*> &Flds){
   if (!isopened())
      {FFileID->open(FFileName, ios::out);isopen = TOutputFile::opened;}
   //SetTitles (Flds);
   TOutputFile::titleline (Flds);
   if (!AddRecords(Rcds))
       return false;
   return true;
 }

/****************************************************************************/
/*                  operateurs d'insertion et d'extraction                  */
/****************************************************************************/


//-----------------------------------------------------------
TInputFile   &operator>>(TInputFile  &InFile, TSequence &Seq){
   switch (Seq.RangeType){
   case TSequence::line : {
         if (InFile.GetRecord (Seq)==NULL)
            cerr <<"GetRecord (); File Error : " << InFile.GetFileName() << endl;
         break;
      }
    case TSequence::colum : {
         if (InFile.GetField (Seq)==NULL)
            cerr <<"GetField (); File Error : " << InFile.GetFileName() << endl;
         break;
      }
   }

   return InFile;
}
//-----------------------------------------------------------
TInputFile   &operator>>(TInputFile  &InFile, Array<TSequence*> &Seqs){
   switch(InFile.cfg){
   case Tiofile::onrcds :{
        if (InFile.GetRecords (Seqs)==NULL)
           cerr <<"GetRecords (); File Error : " << InFile.GetFileName() << endl;
        break;
     }
   case Tiofile::onflds:{
        if (InFile.GetFields (Seqs)==NULL)
           cerr <<"GetFields (); File Error : " << InFile.GetFileName() << endl;
        break;
     }
  }
  return InFile;
}
//-----------------------------------------------------------
TOutputFile &operator<<(TOutputFile &OutFile, TSequence &Seq){
   if (!OutFile.SetRecord(Seq))
       cerr <<"SetRecord (); File Error : " << OutFile.GetFileName() << endl;
   return OutFile;
}
//-----------------------------------------------------------
TOutputFile &operator<<(TOutputFile &OutFile, Array<TSequence*> &Seqs){
	unsigned int idx = 0;
	switch(Seqs[idx]->RangeType) {
		case TSequence::line : {
			if (!OutFile.SetRecords(Seqs))
				cerr <<"SetRecords (); File Error : " << OutFile.GetFileName() << endl;
			break;
		}
		case TSequence::colum: {
			if (!OutFile.SetFields(Seqs))
				cerr <<"SetFields (); File Error : " << OutFile.GetFileName() << endl;
			break;
		}
	}
	return OutFile;
}
//-----------------------------------------------------------
istream &operator>>(ifstream &FileID, TSequence &Seq)
{
  switch(Seq.RangeType){
  case TSequence::line :{
  	FileID >> *Seq.Title;
  	for (unsigned int i=0;i<Seq.nElements;i++){
     	   FileID>> Seq.Value;
    	   Seq.Elements->Set(i, Seq.Value);
     	}
  	break;
     }
  case TSequence::colum:{
  	FileID>> *Seq.Title;
  	for (unsigned int i=0;i<Seq.nElements;i++){
  		ignoreLine(&FileID, STREAM_LINE_MAX_SIZE, LINE_DEFAULT_DELIMITER);
		for (unsigned int i = 0; i < Seq.noID; i++) {
            ignoreLine(&FileID, STREAM_COL_MAX_SIZE, Seq.delimiter);
		}
     	FileID>> Seq.Value;
    	Seq.Elements->Set(i, Seq.Value);
  	}
  	break;
     }
  }
  return FileID;
}
//-----------------------------------------------------------
ostream &operator<<(ostream &FileID, TSequence &Seq)
{
  //FileID<<Seq.noID<<'	'<<Seq.Value<<'	'<<*Seq.Title;
  FileID << *Seq.Title;
  for (unsigned int i=0;i<Seq.nElements;i++){
     Seq.Value = Seq.Elements->Get(i);
     FileID<<'	'<<Seq.Value;
  }
  FileID<<endl;
  return FileID;
}

/****************************************************************************/
/*                                iodata class                              */
/****************************************************************************/
ioData::ioData (TInputFile & inFile, TOutputFile& outFile){
      ifile   = new TInputFile (inFile);
      ofile   = new TOutputFile(outFile);
      //Records = new TSeqsArray (inFile.nRecords());
      //Fields  = new TSeqsArray (inFile.nFields());
      Records = inFile.GetRecords ();
      Fields  = inFile.GetFields  ();
      //load();
}
//-----------------------------------------------------------
ioData::ioData (const char * ifname, const char * ofname, unsigned int nTln, unsigned int nTcol)
{
      ifile   = new TInputFile  ((char*)ifname, nTln, nTcol);
      Records = new TSeqsArray  (ifile->nRecords());
      Fields  = new TSeqsArray  (ifile->nFields());
      if (strcmp(ofname, "--") != 0)
         ofile   = new TOutputFile ((char*)ofname);
      else
         ofile   = NULL;
      load();
}
//-----------------------------------------------------------
ioData::ioData (TInputFile & inFile){
      ifile   = new TInputFile (inFile);
      Records = inFile.GetRecords ();
      Fields  = inFile.GetFields  ();
      ofile   = NULL;
}
 //-----------------------------------------------------------
void ioData::load(void) {
   	ifile->open();
   	ifile->config(Tiofile::onrcds);
   	Records = ifile->GetRecords();
   	ifile->config(Tiofile::onflds);
         Fields = ifile->GetFields();
         ifile->close();
}
//-----------------------------------------------------------
void ioData::write(void){
      if (!(ofile==NULL)){
         ofile->open();
         ofile->titleline(*Fields);
         *ofile<<*Records;
         ofile->close();
      }
}
/****************************************************************************/


void multi_data(TSeqsArray F, const char* name_var,int number,int Lines,vector<vector<double>> &Q){
    vector<double> Qtemp(Lines);
    Q.clear();

    for (int j(0); j < number; ++j) {
        TSequence *PtrSeq = Get(F,j, (char *) name_var);
        for (unsigned int i(0); i < Lines; ++i) {
            Qtemp[i] = PtrSeq->Get(i);
        }
        Q.push_back(Qtemp);
    }

}

void one_data(TSeqsArray F, const char* name_var,int Lines,vector<double> &Q) {

    Q.clear();


    TSequence *PtrSeq = Get(F,(char *) name_var);
    for (unsigned int i(0); i < Lines; ++i) {
        Q.push_back(PtrSeq->Get(i));
    }



}
void one_dataT(TSeqsArray F, const char* name_var,int Lines,vector<double> &Q) {
    Q.clear();
    TSequence *PtrSeq = Get(F,(char *) name_var);
    for (unsigned int i(0); i < Lines; ++i) {
        Q.push_back(PtrSeq->Get(i)+273.15);
    }
}


void Getdata(const char* file_in, vector<double> &A,int &Time,const char *dataName,int inTline,int inTcol, char isep) {
    TInputFile InputFile((char *) file_in, inTline, inTcol, isep);
    TSeqsArray Fields = *InputFile.GetFields();
    TSequence *PtrSeq = Get(Fields, (char *) dataName);
    A.clear();
    Time = InputFile.nRecords();
    if (PtrSeq == NULL)
        cout << " File Error";
    for (int i(0); i < Time; ++i) {
        A.push_back(PtrSeq->Get(i));
    }
}

void Getdata(const char* ifile_in, vector<double> &A,int &Time,const int &index,const int inTline,const int inTcol,const char isep) {
    TInputFile  InputFile((char*)ifile_in, inTline, inTcol, isep);
    TSeqsArray Fields = *InputFile.GetFields();
    TSequence* PtrSeq = Get(Fields, index);
    A.clear();
    Time=InputFile.nRecords();
    if (PtrSeq == NULL)
        cout << " File Error";
    for(int i(0); i<Time; ++i)
    {
        A.push_back(PtrSeq->Get(i));
    };
}
void Getdata_test(const char* ifile_in, vector<double> &A,int &Time,const int &index,const int inTline,const int inTcol,const char isep) {
    fstream fin;

    // Open an existing file
    fin.open(ifile_in, ios::in);

}

void SetOutput(vector <double> &print, char* ifile_out){
    TOutputFile File((char*) ifile_out);
    TSequence Fields(print.size());
    for (int j(0);j<print.size();++j) {
        Fields.Set(j, print[j]);}

    File.SetRecord(Fields);
    // File.AddRecord(Fields);

};

void SetOutput(vector <double> &print,vector <double>& print2, char* ifile_out){
    TOutputFile File((char*) ifile_out);
    TSequence Fields(print.size());
    for (int j(0);j<print.size();++j) {
        Fields.Set(j, print[j]);}

    File.SetRecord(Fields);
    for (int j(0);j<print2.size();++j) {
        Fields.Set(j, print2[j]);}
    File.AddRecord(Fields);

};

void SetOutput(vector<vector <double>> &print, char* ifile_out){
    TOutputFile File((char*) ifile_out);
    Array<TSequence*> Rklds(print.size());
    vector<TSequence> Fields={};

   for (int i(0);i<print.size();i++) {
       TSequence Field1(print[0].size());
       for (int j(0); j < print[0].size(); j++)
           Field1.Set(j, print[i][j]);
       Fields.push_back(Field1);
       Rklds.Set(i,&Fields[i]);
   }
    File.SetFields(Rklds);

};


void SetOutput2(vector<vector <double>> &print, char* ifile_out){
    TOutputFile File((char*) ifile_out);
    TSequence Fields(print[0].size());
    for (int j(0);j<print[0].size();++j) {
        Fields.Set(j, print[0][j]);}
    File.SetRecord(Fields);

    for (int i(1);i<print.size();i++) {
        for (int j(0); j < print[i].size(); ++j) {
            Fields.Set(j, print[i][j]);
        }
        File.AddRecord(Fields);
    }
};


vector<vector<string>> read_record(char* file_input){
//    cout<<file_input<<endl;
    // File pointer
    fstream fin;

    // Open an existing file
    fin.open(file_input, ios::in);

    if(!fin.is_open()) throw runtime_error("Could not open file");
    // Read the Data from the file
    // as String Vector
    vector<string> row;
    vector<vector<string>> table;
    string line, word, temp,all;
    while (fin >> temp) {
        // read an entire row and
        // store it in a string variable 'line'
        getline(fin,all,' ');
        // used for breaking words
        stringstream s(all);

        // read every column data of a row and
        // store it in a string variable, 'word'
        while (getline(s, line)) {
            stringstream s2(line);
            while (getline(s2, word,';')) {
            // add all the column data
            // of a row to a vector
            row.push_back(word);
            }
            if (row.size()>0){
                table.push_back(row);
            }
            row.clear();
        }

    }

    return table;
}


