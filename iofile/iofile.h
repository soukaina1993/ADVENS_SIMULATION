

#ifndef IOFILE_MK_H
#define IOFILE_MK_H
/****************************************************************************/

/*
REMARQUES: 
   TSequences
   * Enregistrement = Ligne de valeurs (hors titres)
   * Field = Colonne de valeurs (hors titres)
   * Par défaut, une TSequence prend les lignes ie les enregistrements !
QUESTIONS/
   1- Dans la classe TSequence, on définit Evaluator mais uniquement avec pour paramètre un TSequence. 
      Dans le cas où l'on souhaite avoir recours à un génôme (comme le MIT), ne faut-il pas définir un Evaluator avec comme paramètre un Génome ? 
      Ou tout est fait via des TSequences (à voir plus loin dans le code, 02-06-2020) ?
*/
#ifdef _MSC_VER
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#endif

#include <fstream>
#include <iostream>
#include <string>
#include "array.h"
#include "dynstring.h"
#include "../math/maths.h"

#define TAILLE 128
#define STREAM_COL_MAX_SIZE 32
#define STREAM_LINE_MAX_SIZE 1024
#define COL_DEFAULT_DELIMITER '\t'
#define LINE_DEFAULT_DELIMITER '\n'
#define USE_IGNORE_LINE_HACK false // pas utilis�: � laisser � false (corrige un bug de compilateur de code::blocks 20, on ne sait trop comment)

/****************************************************************************/
class TSequence;
class Tiofile;
class TInputFile;
class TOutputFile;
class ioData;
//class ioParams;

typedef Array<TSequence *> TSeqsArray;

void ignoreLineHack(ifstream* stream, unsigned int size, const char delimiter = LINE_DEFAULT_DELIMITER);	// use a while loop to next end of line
void ignoreLine(ifstream* stream, unsigned int size, const char delimiter = LINE_DEFAULT_DELIMITER);		// ignore line using "whileloop" or "stream->ignore" depending on USE_IGNORE_LINE_HACK

/****************************************************************************/
class TSequence {
public:
   typedef enum TRangeType {line=0, colum=1} TRangeType; //définir si le tableau est en lignes ou en colonnes 

   typedef float (*Evaluator)(TSequence& Seq); //définition de Evaluator ==> pointeur de fonction vers la fonction objectif qui prend en paramètre des TSequences

   //Constructeurs
   TSequence (unsigned int nElmts, TSequence::TRangeType RangeType = TSequence::line, const char colDelimiter = COL_DEFAULT_DELIMITER)
   {
      this->RangeType = RangeType;
      nElements 	= nElmts;
      Elements 	= new Array<float> (nElements);
      Title 	= new dynstring (STREAM_COL_MAX_SIZE);
      //Title 	= new dynstring (64);
      eval  	= (TSequence::Evaluator)0; 
	  delimiter = colDelimiter;
   }

   TSequence(const TSequence & orig)
   {
      Elements = new Array<float> (orig.nElements);
      Title = new dynstring (orig.Title->size());
      copy(orig);
   }

   //virtual ~TSequence  (void) 	{}
   virtual ~TSequence  () 	{delete Elements; delete Title;}

   //op�rateurs diverses
   friend   std::istream &operator>>(std::istream &FileID, TSequence &Seq);
   friend   std::ostream &operator<<(std::ostream &FileID, TSequence &Seq);
   float & operator[](unsigned int i) {return Value=Elements->Get(i);}
   void operator=(TSequence &Seq)	  {copy (Seq);}

   //fonctions diverses
   unsigned int size () 		{return nElements;}
   float Get(unsigned int pos){
      if (!(pos<nElements)) return 0;
      //return Value=Elements->Get(pos);
      return Elements->Get(pos);
   }

   float Set(unsigned int pos, float val){
      if (!(pos<nElements)) return 0;
      else Elements->Set(pos, val);
      //return Value=val;
      return val;
   }

   float Get(){return Value;}

   //copie la s�quence d'origine (� employer pour 2 seq identiques)
   void copy(const TSequence & orig){
      eval	= orig.eval;
      RangeType 	= orig.RangeType;
      noID 	= orig.noID;
      Value 	= orig.Value;
      Title->copy (*orig.Title);
      nElements 	= orig.nElements;
      for (unsigned int i=0;i<orig.nElements;i++)
         Elements->Set(i, orig.Elements->Get(i));
      delimiter = orig.delimiter;
  }

  //colle la s�quence d'origine (peut �tre employ� pour 2 seq diff�rentes)
   void paste(const TSequence & orig){
      eval	= orig.eval;
      RangeType 	= orig.RangeType;
      noID 	= orig.noID;
      Value 	= orig.Value;
      Title->copy (*orig.Title);
      unsigned int Size = orig.nElements;
      if (!(nElements>Size)) Size = nElements;
      for (unsigned int i=0; i<Size; i++)
         Elements->Set(i, orig.Elements->Get(i));
      delimiter = orig.delimiter;
  }

   //retourne un ptr sur un clone de l'objet
   TSequence * clone()	{return new TSequence(*this);}

   //Ptr de fonction d'�valuation de la s�quence
   TSequence::Evaluator evaluator() 			 {return eval;}
   TSequence::Evaluator evaluator(const TSequence::Evaluator f) {return eval=f;}
   float evaluate () {return this->Value = (*eval)(*this);} // ???

   //-------------------------------------------------
public:

   TSequence::Evaluator eval;	// Ptr de function vers la fonction objectif à évaluer pour chaque génôme (au sens MIT ie ensemble de gènes)
   TRangeType RangeType;		// type de serie (ligne ou colonne)
   Array<float> * Elements;		// Ptr sur le tabeau des donn�es de la sequence
   unsigned int nElements;		// nombre d'elements-donnees de la sequence
   unsigned int noID;			// num�ro identificateur de la s�quence (indice)
   float Value;					// valeur du donnee courante
   dynstring * Title;			// Ptr sur le string courant
   char delimiter;				// separateur de colonnes
};

//-------------------------------------------------
//retourne la s�quence sp�cifi�e par nom ou par indice
//a partir du compteur cptr sp�cifi�
TSequence * Get (const TSeqsArray &Seqs, char * ch, unsigned int cptr=0);
TSequence * Get (const TSeqsArray &Seqs,int d, char * ch, unsigned int cptr=0);
TSequence * Get (const TSeqsArray &Seqs, const unsigned int Index);

//how often is there
bool isexist (const TSeqsArray &Seqs, char * sname);
unsigned int howoften (const TSeqsArray &Seqs, char * sname);


/****************************************************************************/
class Tiofile{
public:
   typedef enum fisopen {opened=0, closed=1} fisopen;
   typedef enum fconfig {onrcds=0, onflds=1} fconfig;

   Tiofile(char * AFileName) {FFileName = AFileName; Tls=new dynstring (STREAM_LINE_MAX_SIZE);}
   virtual ~Tiofile  (void)  {delete Tls;}


   virtual void close (void) {isopen=Tiofile::closed;}
   virtual void open  (void) {isopen=Tiofile::opened;}

   //is it still open
   bool isopened (void) {
      switch(isopen){
      case Tiofile::opened :{return true;  break;}
      case Tiofile::closed :{return false; break;}
      }
      return true;
   }
   //file congiguration
   void config (const Tiofile::fconfig acfg) {cfg = acfg;}

   //Files
   char * GetFileName () {return FFileName;}
   virtual char * SetFileName (char * AFileName) {return FFileName = AFileName;}
   virtual dynstring & title (void)	    	    {return *Tls;}

protected:
   char * FFileName;			//Nom de fichier
   dynstring * Tls;			//ligne de titre
   fisopen isopen;
   fconfig cfg;
};

/****************************************************************************/
class TInputFile : public Tiofile {
public:

   //-------------------------------------------------
   TInputFile (char * AFileName, unsigned int nTline = 1, unsigned int nTcol = 1, const char colDelimiter = COL_DEFAULT_DELIMITER);
   TInputFile (TInputFile & inFile);
   virtual ~TInputFile(void) {
   	delete FFileID;
         delete FRecord; delete FField;
         delete FRcds; delete FFlds;
   }
    //close file
   virtual void close (void) {
    if (isopened()) {
         FFileID->close();
         Tiofile::close();
      }
   }
   //open file
   virtual void open  (void) {
      if (!isopened()) {
         FFileID->open(FFileName, ios::in);
         Tiofile::open();
      }
   }
   //-------------------------------------------------
   //controle et d�placement du curseur
   bool	NextWord ();						//curseur au mot suivant
   bool	NextWord (unsigned int Number);				//mot suivant distant de Nb
   bool	DownWord ();						//mot en dessous
   bool	DownWord (unsigned int Number);				//mot dessous distant de Nb

   //-------------------------------------------------
   //Recherche et extraction dans le fichier d'entr�e

   bool  NextRecord (TSequence &Sequence);				 //ligne suivante
   bool  NextRecord (TSequence &Sequence, unsigned int Number);	 //ligne distante de Nb
   bool  NextField  (TSequence &Sequence);				 //colonne suivant
   bool  NextField  (TSequence &Sequence, unsigned int Number);	 //colonne distante de Nb

   TSequence * GetRecord (TSequence &Sequence);			 //ligne courante
   TSequence * GetRecord (TSequence &Sequence, const dynstring &Str);	 //recherche par nom
   TSequence * GetRecord (TSequence &Sequence, const unsigned int Index);//recherche par index

   TSequence * GetField  (TSequence &Sequence);			 //colonne courante
   TSequence * GetField  (TSequence &Sequence, const dynstring &Str);	 //recherche par nom
   TSequence * GetField  (TSequence &Sequence, const unsigned int Index);//recherche par index

   TSeqsArray * GetRecords (TSeqsArray &ARcds);		 	 //extrait les lignes sp�cifi�es
   TSeqsArray * GetFields  (TSeqsArray &AFlds);	 	 	 //extrait les colonnes sp�cifi�es
   TSeqsArray * GetRecords () {return FRcds;}			 //extrait toutes les lignes
   TSeqsArray * GetFields  () {return FFlds;}			 //extrait toutes les colonnes

   //-------------------------------------------------
   //fonctions diverses
   unsigned int nRecords 	() {return nRcds;} //nombre de lignes (hors lignes de titre)
   unsigned int nFields 	() {return nFlds;} //nombre de colonnes (hors colonnes de titre)
   unsigned int nLines 	() {return nLine;}  //nombre de lignes totales
   unsigned int nColums 	() {return nCol;} //nombre de colonnes totales
   unsigned int nTLines 	() {return nTleLn;} //nombre de lignes de titre
   unsigned int nTColums 	() {return nTleCo;} //nombre de lignes de colonnes 

   char * SetFileName (char * AFileName) {
      FFileName = AFileName; TInputFile::close(); return FFileName;
   }
   void InitCursor() {InitFilePtr(nTleLn - 1, nTleCo - 1);}	//curseur au d�but des seqs.
   dynstring & titleline (const char* Name="Title");		//ligne obtenue a partir des seqs
   dynstring & firstline  ();				//ligne extrait dans le fichier
   dynstring & firstcolum ();				//ligne extrait dans le fichier


   friend TInputFile   &operator>>(TInputFile  &InFile, TSequence &Seq);
   friend TInputFile   &operator>>(TInputFile  &InFile, Array<TSequence*> &Seqs);

protected:
   //-------------------------------------------------

   void InitFilePtr(unsigned int i=0, unsigned int j=0);	//curseur au debut du fichier
   bool Dimension ();					//dimensions du fichier

   unsigned int	GetLineNumber	();			//nbre de ligne
   unsigned int	GetColumNumber	();			//nbre de colonne
   unsigned int	GetWordNumber	();			//nombre total de mot

   // latchtech
   void debug(const char* tag); // debug next word

   //-------------------------------------------------

private:
   ifstream * FFileID;			//Ptr sur l'identificateur de fichier
   TSequence * FRecord;			//Ptr sur un enregistrement
   TSequence * FField;			//Ptr sur un champs
   Array <TSequence*> * FRcds;		//Ptr sur un tableau de pointeurs sur les enregistrements (valeurs en ligne), 
   Array <TSequence*> * FFlds;		//Ptr sur un tableau de fields

   unsigned int nTleLn;			//nombre de lignes de titre
   unsigned int nTleCo;			//nombre de colonnes de titre
   unsigned int nRcds; 			//nombre d'enregistrement
   unsigned int nFlds;			//nombre de champs de donn�es
   unsigned int nLine;			//nombre de ligne
   unsigned int nCol;			//nombre de colonne

   char delimiter;				//delimiter par defaut utilis� pour l'extraction de donnn�es
};

/****************************************************************************/
template<class T>
TOutputFile &operator<<(TOutputFile &OutFile, T& a);

class TOutputFile : public Tiofile{
public:
   TOutputFile (char * AFileName);
   ~TOutputFile (void) {delete FFileID;}

   //close file
   void close (void) {
    if (isopened()) {
         FFileID->close();
         Tiofile::close();
      }
   }
   //open file
   void open  (void) {
      if (!isopened()) {
         FFileID->open(FFileName, ios::out);
         Tiofile::open();
      }
   }
   //Records
   bool SetRecord  (TSequence &Sequence);   		//�crit la seq. en format ligne
   bool AddRecord  (TSequence &Sequence); 		//ajoute la seq. en format ligne
   bool SetRecords (Array<TSequence*> &Rcds);	//�crit les seqs. en format ligne
   bool AddRecords (Array<TSequence*> &Rcds);	//ajoute les seqs. en format ligne

   //Fields
   bool SetFields  (Array<TSequence*> &Flds);	//�crit les seqs. en format colonne
   bool AddFields  (Array<TSequence*> &Flds);	//ajoute les seqs. en format colonne
   //void SetTitles  (Array<TSequence*> &Flds);
   //void SetValues  (Array<TSequence*> &Flds);

   bool Set2DArray (Array<TSequence*> &Rcds, Array<TSequence*> &Flds);

   //File name & title line
   char* SetFileName ( const char * AFileName) {
      FFileName = (char*) AFileName; TOutputFile::close(); return FFileName;
   }
   void titleline  (const Array<TSequence*> &Flds, const char* Name = "Title");
   void titleline  (dynstring & tleline);

   //FileID
   ofstream & FileID(void) {return * FFileID;}


   //friends
   friend TOutputFile &operator<<(TOutputFile &OutFile, Array<TSequence*> &Rcds);
   friend TOutputFile &operator<<(TOutputFile &OutFile, TSequence &Seq);


    void writevalues  (Array<TSequence*> &Flds);

private:
   ofstream * FFileID;			//Ptr sur l'identificateur de fichier


};
//sortie de la classe 29-03-2020
template<class T>
   TOutputFile &operator<<(TOutputFile &OutFile, T& a){
   	*OutFile.FFileID<<a;
   	return OutFile;
   }
/****************************************************************************/
TInputFile   &operator>>(TInputFile  &InFile, TSequence &Seq);
TInputFile   &operator>>(TInputFile  &InFile, Array<TSequence*> &Seqs);

TOutputFile  &operator<<(TOutputFile &OutFile, TSequence &Seq);
TOutputFile  &operator<<(TOutputFile &OutFile, Array<TSequence*> &Seqs);

std::istream &operator>>(std::ifstream &FileID, TSequence &Sequence);
std::ostream &operator<<(std::ostream &FileID, TSequence &Sequence);


/****************************************************************************/
class ioData{
public:
   ioData (TInputFile & inFile, TOutputFile& outFile);
   ioData (const char * ifname, const char * ofname = "--", unsigned int nTln=1, unsigned int nTcol=1);
   ioData (TInputFile & inFile);
   virtual ~ioData(void) {delete Records; delete Fields;}

   unsigned int count (void) {return ifile->nRecords();}
   void load(void);
   virtual void write(void);
   //-------------------------------------------------
   TSeqsArray & records  () 		   {return *Records;}
   TSeqsArray & fields   () 		   {return *Fields;}
   //-------------------------------------------------
   //get sequence
   TSequence & record    (char * rcdname) 	   {return *Get(*Records, rcdname);}
   TSequence & record    (unsigned int rcdidx)  {return *Get(*Records, rcdidx);}
   TSequence & field     (char * fldname) 	   {return *Get(*Fields , fldname);}
   TSequence & field     (const char * fldname)      {return field((char*)fldname);}

   TSequence & field     (unsigned int fldidx)  {return *Get(*Fields , fldidx);}
   //-------------------------------------------------
   //get sequence id
   unsigned int idRecord (char * rcdname){
   	return record(rcdname).noID-ifile->nTLines();
   }
   unsigned int idField  (char * fldname){
   	return field (fldname).noID-ifile->nTColums();
   }
   //-------------------------------------------------
   //get specified value
   float value (unsigned int i, unsigned int j){return record(i).Get(j);}
   float value(const char* rcdname, char* fldname) { return value((char*)rcdname, fldname); }
   float value (char * rcdname, char * fldname){
   	return record(rcdname).Get(idField(fldname));
   }
   //set specified value
   float value (char * rcdname, char * fldname, float val) {
   	return record(rcdname).Set(idField(fldname), val);
   }
   float value (unsigned int i, unsigned int j, float val) {
   	return record(i).Set(j,val);
   }
   //-------------------------------------------------
   char * ifilename () 		{return ifile->GetFileName();}
   char * ofilename () 		{return ofile->GetFileName();}
   char * ofilename (char * fn) 	{return ofile->SetFileName(fn);}
   TInputFile & ifileID  	() 	{return *ifile;}
   TOutputFile & ofileID   () 	{return *ofile;}

public:
   TInputFile  * ifile;
   TOutputFile * ofile;
   TSeqsArray  * Records;
   TSeqsArray  * Fields ;
};

void one_data(TSeqsArray F, const char* name_var,int Lines,vector<double> &Q);
void one_dataT(TSeqsArray F, const char* name_var,int Lines,vector<double> &Q);
void multi_data(TSeqsArray F, const char* name_var,int number,int Lines,vector<vector<double>> &Q);

void Getdata(const char* ifile_in, vector<double> &data,int &Time,const char *dataName,int nTline=1,int nTcol=1, char isep = '\t');
void Getdata(const char* ifile_in, vector<double> &data,int &Time,const int &index=1, int nTline=1,int nTcol=1, char isep = '\t');

void SetOutput(vector<double> &vector1, char* ifile_out=(char*)"Results.txt");
void SetOutput(vector<double> &vector1, vector <double> &vector2, char* ifile_out=(char*)"Results.txt");
void SetOutput(vector<vector <double>> &print, char* ifile_out=(char*)"Results.txt");
void SetOutput2(vector<vector <double>> &print, char* ifile_out=(char*)"Results.txt");
vector<vector<string>> read_record(char* file_input);

/****************************************************************************/




#endif

