//"GAMixedGenome Header v2.0 Malick kane 2000\n";
/* ----------------------------------------------------------------------------
  mixed.h
  akm 1April97
  Copyright (c) 1997 Adam K Molyneaux
                     all rights reserved

 DESCRIPTION:
 A mixed integer and real genome.

 MODIFICATIONS:
 modified by malick.kane@epfl.ch, 2000
 Extension: A mixed integer, real, binary and discretized real genome
---------------------------------------------------------------------------- */
#ifndef _GA_MIXED_GENOME_H_
#define _GA_MIXED_GENOME_H_

#include "GARealGenome.h"

/****************************************************************************/
class MixedAlleleSet{
public:
   MixedAlleleSet ()
   {
      realals	= NULL; ncon=0;
      disals	= NULL; ndis=0;
      binals	= NULL; nbin=0;
      intals	= NULL; nint=0;

   }
   virtual ~MixedAlleleSet(void) 	     {} //

   GARealAlleleSetArray & realalleles () {return *realals;}
   GARealAlleleSetArray & binalleles  () {return *binals ;}
   GARealAlleleSetArray & intalleles  () {return *intals ;}
   GARealAlleleSetArray & disalleles  () {return *disals ;}

   unsigned int  nContinious	  () {return ncon;}
   unsigned int  nInteger		  () {return nint;}
   unsigned int  nDiscretized	  () {return ndis;}
   unsigned int  nBinary		  () {return nbin;}

public:

   GARealAlleleSetArray * realals;       	//Continious allele set
   GARealAlleleSetArray * disals;       	//discretized real allele set
   GARealAlleleSetArray * binals;       	//binary allele set (0,1 enumerated)
   GARealAlleleSetArray * intals;       	//Integer allele set

   unsigned int ncon;
   unsigned int nint;
   unsigned int ndis;
   unsigned int nbin;

};
/****************************************************************************/

class MixedGenome : public GAGenome {
public:
  GADefineIdentity("MixedGenome", 245);

  static void  MixedInitializer (GAGenome&);
  static int   MixedMutator     (GAGenome&, float);
  static float MixedComparator  (const GAGenome&, const GAGenome&);
  static int   MixedCrossover   (const GAGenome&, const GAGenome&,
				   GAGenome*, GAGenome*);
public:

  MixedGenome(MixedAlleleSet & mxls,
	     GAGenome::Evaluator f=NULL, void* u=NULL);
  MixedGenome(const MixedGenome & orig);
  MixedGenome& operator=(const GAGenome& g);
  virtual ~MixedGenome();
  virtual GAGenome* clone(GAGenome::CloneMethod) const ;
  virtual void copy (const GAGenome & c);
  virtual int  equal(const GAGenome& g) const;
  virtual int  read(std::istream & is);
  virtual int  write(std::ostream & os) const;

  GARealGenome & realgenome() const {return *rgnom;}
  GARealGenome & intgenome () const {return *ignom;}
  GARealGenome & disgenome () const {return *dgnom;}
  GARealGenome & bingenome () const {return *bgnom;}


protected:
  GARealGenome *rgnom;
  GARealGenome *dgnom;
  GARealGenome *ignom;
  GARealGenome *bgnom;

  //MixedAlleleSet * mxls;
  static unsigned int nr;
  static unsigned int ni;
  static unsigned int nd;
  static unsigned int nb;

};
/****************************************************************************/


#endif

