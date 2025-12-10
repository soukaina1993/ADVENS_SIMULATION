/*---------------------------------------------------------------------------
 DESCRIPTION:

---------------------------------------------------------------------------- */
#ifndef GAO_H
#define GAO_H

#ifdef _MSC_VER
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
#endif

#include "ga/ga.h"

#include <iostream>
#include "../iofile/dynstring.h"
#include "../iofile/iofile.h"
//#include "BasicES.h"
#include "gaiodata.h"

/****************************************************************************/
class gao; //11-05-2020: Classe mère (non instantiable?) pour la déclaration de la plateforme d'optimisation
class mxgao; //11-05-2020: Plateforme d'optimisation pour des valeurs d'entrées mixtes
class mogao; //11-05-2020: Gestion d'une optimisation multiobjectif (avec des varibles indépendantes mixtes)

/****************************************************************************/
float objective    (GAGenome& g);
float mxobjective  (GAGenome& g); 
//float edcomparator (const GAGenome& g1, const GAGenome& g2);

const int ndiplay = 10;

/****************************************************************************/
class gao{
public:

   typedef enum type 	{simple=0, sstate=1, stggle=2, inc=3, deme=4, complex=5, mobjective=6} type;
   typedef enum extramum 	{max=0, min=1} extramum;

   //--------------------------------------------------
   gao(const char * ifname = "gadata.txt", const char * ofname = "garesults.txt", type GA=gao::simple);
   gao(TInputFile & inFile, TOutputFile& outFile, type GA=gao::simple);
   gao(TInputFile & inFile, type GA=gao::simple);

   virtual ~gao(void){
      data->ofile->close();
      delete data;
      if (!(InitGen==NULL)) delete InitGen;
   }

   //-------------------------------------------------
   void maximize (BasicES & thecpnt, type GA=gao::simple);
   void maximize (GAGenome::Evaluator  fg=(GAGenome::Evaluator)0 , type GA=gao::simple);
   void maximize (TSequence::Evaluator fs=(TSequence::Evaluator)0, type GA=gao::simple);

   void minimize (BasicES & thecpnt, type GA=gao::simple);
   void minimize (GAGenome::Evaluator  fg=(GAGenome::Evaluator)0 , type GA=gao::simple);
   void minimize (TSequence::Evaluator fs=(TSequence::Evaluator)0, type GA=gao::simple);
   //-------------------------------------------------
   virtual void init  ();
   virtual void start ();
   virtual void stop  ();
   //-------------------------------------------------
   GAGeneticAlgorithm * selectga (GAGenome& g, type GA=gao::simple);
   extramum search (extramum optimum=gao::max)  {return opt = optimum;}
   gaData * userdata ()			   {return data;}

   //-------------------------------------------------
   friend float objective (GAGenome& g);		//friend float objective (GAGenome& g);

   virtual void write   (const GAGenome& g);		//write out a genome
   virtual void write   (const GAPopulation& p);	//write out whole population
   virtual void display (const GAGenome& g);		//display one genome
   virtual void display (const GAPopulation& p);	//display whole population


   virtual TSequence * memorize (const GAGenome& g);	//keep memory a genome
   bool historic 	  (const char * ofname="historic.txt");	//historic of the best genomes

   						//one gene sensitivity
   virtual bool sensitivity   (unsigned int noGene, const char * ofname="sensitivity.txt");
   						//the n-first genes sensitivity
   virtual bool sensitivities (unsigned int nbGene, const char * ofname="sensitivity.txt");
   unsigned int nSensitivity  (unsigned int nbGene)	{return nSen = nbGene;}
   unsigned int nSensitivity  (void)		{return nSen;}

   //-------------------------------------------------
   char * resultfile   ()		  {return data->ofile->GetFileName ();}
   char * resultfile   (const char* ofname) {return data->ofile->SetFileName (ofname); }
   //-------------------------------------------------
   TInputFile & ifileID   ()	  {return data->ifileID ();}
   TOutputFile & ofileID   () 	  {return data->ofileID ();}
   //-------------------------------------------------

   TSequence &  initgenome  (void) 	{return *InitGen;} //the genome used for the initialization

public:
   type ThisGA;
   extramum opt;
   gaData * data;
   GAGenome * chrom;
   GARealAlleleSetArray * alleles;
   GAParameterList * pm;
   GAGeneticAlgorithm * ga;
   GAGenome::Evaluator  obj; //10-05-2020: FONCTION OBJECTIF !! (définie dans la classe GAGenome)
   static gaParams * setting;
   unsigned int nSen;				//Number of sensitivity variables
   TSequence * InitGen;				//Init Genome Pointer

};

/****************************************************************************/

class mxgao: public gao{
public:

   //--------------------------------------------------
   mxgao(const char * ifname = "gamxdata.txt", const char * ofname = "garesults.txt", type GA=gao::complex);
   mxgao(TInputFile & inFile, TOutputFile& outFile, type GA=gao::complex);
   mxgao(TInputFile & inFile, type GA=gao::complex);
   virtual ~mxgao(void){
      mxdata->ofile->close();
   }
   //-------------------------------------------------
   virtual void init  ();
   virtual void start ();
   friend float mxobjective (GAGenome& g);
   void write (const GAGenome& g);			//write out a genome
   void display (const GAGenome& g);		//display one genome
   virtual void write (const GAPopulation& p);	//write out whole population
   virtual void display (const GAPopulation& p);	//display whole population
   virtual TSequence * memorize(const GAGenome& g);	//keep memory a genome
   //-------------------------------------------------

public:

   MixedAlleleSet * mxalleles;
   gaMixedData 	* mxdata;

};

/****************************************************************************/
struct gaotype {
   char * Name;				//Name of the objective
   gao::type GA;				//Modele of the used GA
   gao::extramum EX;			//The direction of the optimum (min or max)
   dynstring * rfname;			//result file name
   dynstring * hfname;			//historical file name
   dynstring * sfname;			//sensitivity file name

};
/****************************************************************************/
class mogao: public mxgao{
public:

   //--------------------------------------------------
   mogao(const char * ifname = "gamxdata.txt", const char * ofname = "garesults.txt", unsigned int nfunc=1);
   mogao(TInputFile & inFile, TOutputFile& outFile, unsigned int nfunc=1);
   mogao(TInputFile & inFile, unsigned int nfunc=1);

   virtual ~mogao(void){
      for (unsigned int i=0; i<nMdl; i++)
      	if (!(gamdl->Get(i)==NULL)) {
      	   delete gamdl->Get(i)->rfname;
      	   delete gamdl->Get(i)->hfname;
      	   delete gamdl->Get(i)->sfname;
      	   delete gamdl->Get(i);
      	}
      delete gamdl;
   }
   //-------------------------------------------------
   void Set (BasicES & thecpnt, char*objname, type GA=gao::sstate, extramum EX=gao::min);
   void Set (TSequence::Evaluator fs, char*objname, type GA=gao::sstate, extramum EX=gao::min);
   gaotype & modele (unsigned int i) {return *gamdl->Get(i);}
   gaotype & modele (unsigned int i, gaotype * mdl) {gamdl->Set(i, mdl); return *mdl;}

   void filename (unsigned int i);
   void optimize (unsigned int NoTarget);
   void optimize ();

   //-------------------------------------------------
   void init  ();
   void start ();
   //void stop  ();
   friend float moobjective (GAGenome& g);
   friend float moutility 	 (GAGenome& g);

   void write (const GAGenome& g);			//write out a genome
   void display (const GAGenome& g);		//display one genome
   virtual void write (const GAPopulation& p);	//write out whole population
   virtual void display (const GAPopulation& p);	//display whole population
   virtual TSequence * memorize(const GAGenome& g);	//keep memory a genome

   //-------------------------------------------------
public:

   gaMoData 	 * modata;
   Array<gaotype*> * gamdl;			//ga modele array
   unsigned int nMdl;				//size of ga modele array
   unsigned int iMdl;				//lengh of ga modele


};
/****************************************************************************/
#endif //GAO_H