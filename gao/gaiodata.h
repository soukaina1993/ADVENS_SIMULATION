/*---------------------------------------------------------------------------
 DESCRIPTION: (le 09-05-2020)
 * ifname: input file name. Il définit les allèles possibles d'un individu (les intervalles de valeurs des gènes/varibles indpendantes qu'un gène peut avoir). 
 * nTln: nombre de lignes de titre
 * nTcol: nombre de colonnes de titre
 * Ce fichier regroupe les classes permettant de gérer les données de l'algoruthme génétique. 
 * gaiodata permet la gestion des données d'entrées et de sorties de l'ogorithme génétique choisi. 
 * gaData: variables indépendentes, valeurs du génome
 * gaParam: Paramètres de l'algorithme (génération, population ...)
 * gaMixedData: gestion de variables indépendantes de "type" différents (continues, discrètes, binaires)

---------------------------------------------------------------------------- */
#ifndef GAIODATA_H
#define GAIODATA_H
#include <cmath>
#include "ga/ga.h"

#include <iostream>
#include "../iofile/dynstring.h"
#include "../iofile/iofile.h"
#include "../BasicEnergySystem/BasicES.h"

/****************************************************************************/

class ioData;  // 09-05-2020: Cette classe est définie dans iofile.
class gaData;  /* 10-05-2020: Données et informations sur le génôme des individus. Regroupe les informations, valeurs caractéristiques des variables indépendantes ie des gènes qu'un individu peut avoir. 
                              Les variables indépendantes constituent les gènes des chromosomes (c'est le génome de l'individu, son patrimoine génétique). 
                              On associe à chaque variable indépendante (gène) un allèle (plage valeurs que le gène peut prendre) et une valeur d'initalisation (ici nommée "genome").
                              A une ligne correspond un gène possible (avec l'allèle ie l'intervalle de valeur que peut prendre ce gène).
                              Un individu est représenté par l'ensemble des gènes du fichier dont il prend une valeur définie par l'algorithme génétique. 
                              Classe mère de "gaMixedData" */ 
class gaParams; // 10-05-2020: paramètre de l'algorithme d'optimisation
class gaMixedData; // 10-05-2020: Classe fille de gaData et MixedAlleleSet. Elle permet de gérer des variables indépendantes réelles, entières, binaires et/ou discrètes pour un même génôme. 
class gaMoData;  //10-05-2020: Gestion de plusieurs fonctions objectifs sur des données mixtes

/****************************************************************************/
class gaParams: public ioData
{ //09-05-2020: il ne manque pas un public ? FRO : ajoute
public:

   gaParams (const char * ifname, unsigned int nTln=1, unsigned int nTcol=1):
   ioData(ifname, "--", nTln, nTcol) {
   	params = new GAParameterList;
   }
   gaParams (TInputFile & ifile): ioData(ifile) {params = new GAParameterList;}

   virtual ~gaParams(void){delete params;}

   //Get position value
   int   zpop  (const char * ganame) {return (int)trunc (value("PopulationSize", (char*)ganame));}
   int   ngen  (const char * ganame) {return (int)trunc (value("nGeneration"   , (char*)ganame));}
   float pcros (const char * ganame) {return 	 value("pCrossover"    , (char*)ganame) ;}
   float pmut  (const char * ganame) {return 	 value("pMutation"     , (char*)ganame) ;}
   int   nrepl (const char * ganame) {return (int)trunc (value("nReplacement"  , (char*)ganame));}
   int   nmig  (const char * ganame) {return (int)trunc (value("nMigration"    , (char*)ganame));}
   int   npop  (const char * ganame) {return (int)trunc (value("nPopulation"   , (char*)ganame));}
   int   sfreq (const char * ganame) {return (int)trunc (value("ScoreFrequency", (char*)ganame));}
   int   ffreq (const char * ganame) {return (int)trunc (value("FlushFrequency", (char*)ganame));}

   //-------------------------------------------------
   GAParameterList * registered (const char * ganame);
public:
   GAParameterList * params;
};

/****************************************************************************/
class gaData:public ioData{ // 10-05-2020: informations, valeurs caractéristiques des variables indépendantes. Les variables indépendantes constituent les gènes des chromosomes (c'est le génome de l'individu, son patrimoine génétique). On associe à chaque variable indépendante (gène) un allèle (plage valeurs que le gène peut prendre) et une valeur d'initalisation (ici nommée "genome")
public:
   //-------------------------------------------------
   gaData (TInputFile & inFile, TOutputFile& outFile);
   gaData (TInputFile & inFile);
   gaData (const char * ifname, const char * ofname="--", unsigned int nTln=1, unsigned int nTcol=1); //fichier de sortie par défaut ("--")// nombre titre de ligne et nombre titre de colonne 

   virtual ~gaData(void)	{
   			   delete ConAlles;
   			   for (unsigned int i=0; i<nGnom; i++)
   			      if (!(memo->Get(i)==NULL)) delete memo->Get(i);
     			   delete memo;
     			}
   friend gaData & operator<<(gaData &Data, TSequence &Seq);
   //Get specified sequence
   TSequence &   isopt    (void) 		{return field("isopt") ;}
   TSequence &   lower    (void) 		{return field("lower") ;}
   TSequence &   upper    (void) 		{return field("upper") ;}
   TSequence &   weight   (void) 		{return field("weight");}
   TSequence &   genome   (void) 		{return field("genome");}


   //Get position value
   float isopt  (unsigned int i)  		{return isopt ().Get(i);}
   float lower  (unsigned int i)		{return lower ().Get(i);}
   float upper  (unsigned int i)		{return upper ().Get(i);}
   float weight (unsigned int i)		{return weight().Get(i);}
   float genome (unsigned int i)		{return genome().Get(i);}

   //Set position value
   float isopt  (unsigned int i, float f)    {return isopt ().Set(i,f);}
   float lower  (unsigned int i, float f)    {return lower ().Set(i,f);}
   float upper  (unsigned int i, float f)    {return upper ().Set(i,f);}
   float weight (unsigned int i, float f)    {return weight().Set(i,f);}
   float genome (unsigned int i, float f)    {return genome().Set(i,f);}

   //-------------------------------------------------

   //float evaluate ();
   TSequence::Evaluator evaluator(const TSequence::Evaluator f);
   BasicES * evaluator(BasicES & thecpnt);
   virtual float evaluate ();


   virtual gaData & store (GAGenome&  g);		//genome recordered
   virtual GARealAlleleSetArray * registered (); 	//define all the alleles

   //-------------------------------------------------
   //write data
   virtual void title (void);
   void write   (void)	{*this<<genome();}
   void display (void)	{std::cout<<genome().noID<<"\t"<<genome().Value<<"\t"<<genome();}


   //-------------------------------------------------
   //memorization of the best genome
   unsigned int nGenome (unsigned int gn);
   unsigned int nGenome (void)		{return nGnom;}
   TSequence *  memorize(void);

   //-------------------------------------------------
   //historic of the components
   bool historic (TOutputFile & ofile);			//for all the memorized genomes
   bool historic (TOutputFile & ofile, unsigned int noInd);	//for the specified genome


   //-------------------------------------------------
   //sensitivity analysis
   virtual
   bool sensitivity (TOutputFile & ofile, unsigned int noGene);//for the specified gene



public:

   GARealAlleleSetArray * ConAlles;       		//continious allele set
   TSequence::Evaluator   eval; 			//function pointer
   BasicES * 	        cpnt; 			//component pointer
   float score;					//objective value

   Array <TSequence*> * memo;			//Array of genomes
   unsigned int nGnom;				//Number of genome
   unsigned int iGnom;				//current genome




};
gaData  &operator<<(gaData  &Data, TSequence &Seq);	//write file one sequence

/****************************************************************************/
class gaMixedData:public gaData, public MixedAlleleSet { 
public:
   //-------------------------------------------------
   gaMixedData (TInputFile & inFile);
   gaMixedData (TInputFile & inFile, TOutputFile& outFile);
   gaMixedData (const char * ifname, const char * ofname= "--", unsigned int nTln=1, unsigned int nTcol=2);
   virtual ~gaMixedData(void) 	     {distroy();} //
   TSequence &   Continious	 ()  {return * scon;}
   TSequence &   Integer		 ()  {return * sint;}
   TSequence &   Discretized	 ()  {return * sdis;}
   TSequence &   Binary		 ()  {return * sbin;}
   //-------------------------------------------------
   gaData & store (GAGenome&  g);		//gene recordered
   GARealAlleleSetArray * registered (); 	//define all the alleles
   void title (void);
   MixedAlleleSet * mixedalleles ();

protected:
   void create 	();
   void distroy 	();
   unsigned int nUsed () ;
   unsigned int howoften (char * sname);
   unsigned int howoften(const char* sname) { return howoften((char*)sname); }
   unsigned int RealRegister ();
   unsigned int IntRegister  ();
   unsigned int DisRegister  ();
   unsigned int BinRegister  ();
   unsigned int RealStore (GAGenome&  g);
   unsigned int DisStore  (GAGenome&  g);
   unsigned int IntStore  (GAGenome&  g);
   unsigned int BinStore  (GAGenome&  g);
   TSequence  & sequence  (TSequence &seq, char * sname) ;
   TSequence& sequence(TSequence& seq, const char* sname) { return sequence(seq, (char*)sname); }
   virtual bool sensitivity (TOutputFile & ofile, unsigned int noGene);

public:
   TSequence * scon;			//sequence of continious values
   TSequence * sint;			//sequence of integer values
   TSequence * sdis;			//sequence of discretized values
   TSequence * sbin;			//sequence of binary values

};

/****************************************************************************/
struct TTarget {
   float min; 				//objective min value
   float max;				//objective max value
   float fOpt;				//objective opt value
   float Wght;				//relative weight for objective
   float Wnf;				//normalizing factor
   					//evaluator
   TSequence::Evaluator eval;
   unsigned int tgID;			//target ID
   char * name;				//Name of the objective
};

/****************************************************************************/
class gaMoData:public gaMixedData {
public:
   //-------------------------------------------------
   gaMoData (TInputFile & inFile, unsigned int nfunc=1);
   gaMoData (TInputFile & inFile, TOutputFile& outFile, unsigned int nfunc=1);
   gaMoData (const char * ifname, const char * ofname= "--", unsigned int nfunc=1,
   	    unsigned int nTln=1, unsigned int nTcol=2);

   virtual ~gaMoData(void)	{distroy ();}

   //-------------------------------------------------
   TSequence & vCurrent  (unsigned int i) {return *VdecAr->Get(i);}	//get ieme Decision variables
   TSequence & pCurrent  (unsigned int i) {return *PoffAr->Get(i);}	//get ieme function vector
   TSequence * payoff    (void); 					//all objectives for a same genome
   BasicES * Set (BasicES & thecpnt, char*objname);			//set current component
   TSequence::Evaluator Set (const TSequence::Evaluator f, char*objname); //set current evaluator
   TTarget 	& target	(unsigned int i) {return *tgets->Get(i);}	//get ieme target
   Array<TTarget*>& targets (void);					//get all target

   void  evaluator(unsigned int ii);				//set evaluator of ieme target
   float utility  (unsigned int ii);				//evaluate the ieme target
   float minimaxi (void);						//normalized function

   //-------------------------------------------------
   //write data
   void motitle (void);
   void mowrite (void);
   //void display (void)	{cout<<genome().noID<<"\t"<<genome().Value<<"\t"<<genome();}


protected:
   void create 	();
   void distroy 	();


public:

   Array <TSequence*> * VdecAr;			//Array of decision variable
   Array <TSequence*> * PoffAr;			//Payoff table
   unsigned int nPoff;				//size of payoff table
   unsigned int iPoff;				//Current element of payoff table

   Array <TTarget*> * tgets;			//array of targets
   unsigned int nObj;				//Number of Objective function
   unsigned int iObj;				//Current Objective function



};



/****************************************************************************/
#endif //GAIODATA_H