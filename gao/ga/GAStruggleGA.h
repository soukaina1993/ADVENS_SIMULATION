/* ----------------------------------------------------------------------------
  Struggle.h
  mbwall jan96
  Copyright (c) 1995-1996 Massachusetts Institute of Technology
                          all rights reserved

  Modified.
  3April97 Adam Molyneaux 
  
  Modified.
  27July2000 malick.kane@epfl.ch 

---------------------------------------------------------------------------- */
#ifndef _StruggleHGAH_h_
#define _StruggleHGAH_h_

#include "ga.h"
#ifdef HOOPS
#include <hc.h>
#endif

typedef int (*DrawPopulation)(GAPopulation & pop);

#define UNIFORM_WITH_CORRECTION (0)
#define UNIFORM_MAKEPOP (1)
#define PROPORTIONAL  (2)

int DrawPopulationDefault(GAPopulation &pop);

class DistanceScaling : public GAScalingScheme {
public:
  GADefineIdentity("DistanceScaling", 252);
  DistanceScaling() {}
  DistanceScaling(const DistanceScaling & arg) { copy(arg); }
  DistanceScaling & operator=(const GAScalingScheme & arg)
    { copy(arg); return(*this); }
  virtual ~DistanceScaling() {}
  virtual GAScalingScheme* clone() const { return new DistanceScaling(*this); }
  virtual void evaluate(const GAPopulation & p);
  virtual void copy(const GAScalingScheme & arg)
    { if(&arg != this && sameClass(arg)) GAScalingScheme::copy(arg); }
  void assignInd(GAGenome *g) { genome=g; }
  GAGenome * getInd(void) { return(genome); }

protected:
  GAGenome *genome;
};


class StruggleHGA : public GAGeneticAlgorithm {
public:
  GADefineIdentity("StruggleHGA", 237);
  static GAParameterList& registerDefaultParameters(GAParameterList&);
  
  StruggleHGA(const GAGenome& g);
  StruggleHGA(const GAPopulation& p);
  virtual ~StruggleHGA() { delete child; }
  virtual void initialize(unsigned int seed=0);
  virtual void step();
  StruggleHGA & operator++() { step(); return *this; }
  inline int SelectMethod(int ss) {_SelectMethod=ss;return(_SelectMethod);}
  inline int SelectMethod(void) {return(_SelectMethod);}
  void Display(DrawPopulation p) {drawpop=p;}
  inline float LocalChoiceProbability(float p) {pLocal=p;return(pLocal);}
  inline float LocalChoiceProbability() {return(pLocal);}

protected:
  GAPopulation & MakeChildPop(GAPopulation &cpop);
  GAGenome *child;
  int _SelectMethod;
  float pLocal;    // prob of local selection
  DrawPopulation drawpop;	/* function to draw population */
};

/*#ifndef NO_STREAMS
inline ostream & operator<< (ostream & os, StruggleHGA & arg)
{arg.write(os); return(os);}
inline istream & operator>> (istream & is, StruggleHGA & arg)
{arg.read(is); return(is);}
#endif
*/
#endif


