#include <iostream>
#include "GARealGenome.h"
#include "GAMixedGenome.h"

unsigned int MixedGenome::nr;
unsigned int MixedGenome::ni;
unsigned int MixedGenome::nd;
unsigned int MixedGenome::nb;


int GAListBlendCrossover(const GAGenome& p1, const GAGenome& p2,
		     GAGenome* c1, GAGenome* c2);
//---------------------------------------------------------------
MixedGenome::MixedGenome(MixedAlleleSet & mxls, GAGenome::Evaluator f, void* u):
GAGenome(MixedInitializer,
	MixedMutator,
	MixedComparator)
{
    evaluator(f);
    userData(u);
    crossover(MixedCrossover);
    MixedGenome::nr	= mxls.nContinious();
    MixedGenome::ni	= mxls.nInteger();
    MixedGenome::nd	= mxls.nDiscretized();
    MixedGenome::nb	= mxls.nBinary();

    if(MixedGenome::nr){
       rgnom = new GARealGenome (mxls.realalleles(), f, u);
       rgnom->crossover(GARealUniformCrossover);
    }
    if(MixedGenome::ni){
       ignom = new GARealGenome (mxls.intalleles(),  f, u);
       ignom->crossover(GARealUniformCrossover);
       ignom->mutator(GARealSwapMutator);
    }
    if(MixedGenome::nd){
       dgnom = new GARealGenome (mxls.disalleles(),  f, u);
       dgnom->crossover(GARealUniformCrossover);
       dgnom->mutator(GARealSwapMutator);
    }
    if(MixedGenome::nb){
       bgnom = new GARealGenome (mxls.binalleles(),  f, u);
       bgnom->crossover(GARealUniformCrossover);
       bgnom->mutator(GARealSwapMutator);

    }

}
//---------------------------------------------------------------
MixedGenome::MixedGenome(const MixedGenome & orig) {

  if(MixedGenome::nr)	rgnom = new GARealGenome (orig.realgenome());
  if(MixedGenome::ni)	ignom = new GARealGenome (orig.intgenome());
  if(MixedGenome::nd)	dgnom = new GARealGenome (orig.disgenome());
  if(MixedGenome::nb)	bgnom = new GARealGenome (orig.bingenome());

  copy(orig);
}
//---------------------------------------------------------------
MixedGenome& MixedGenome::operator=(const GAGenome& g) { copy(g); return *this; }

MixedGenome::~MixedGenome() {

   if(MixedGenome::nr)	delete rgnom;
   if(MixedGenome::ni)	delete ignom;
   if(MixedGenome::nd)	delete dgnom;
   if(MixedGenome::nb)	delete bgnom;

}
//---------------------------------------------------------------
GAGenome* MixedGenome::clone(GAGenome::CloneMethod) const {
  return new MixedGenome(*this);
}

void MixedGenome::copy(const GAGenome & c){
  if(&c != this && sameClass(c)){
    GAGenome::copy(c);
    MixedGenome & bc = (MixedGenome &)c;
    /*
    if(MixedGenome::nr)	rgnom->copy(*(bc.rgnom));
    if(MixedGenome::ni)	ignom->copy(*(bc.ignom));
    if(MixedGenome::nd)	dgnom->copy(*(bc.dgnom));
    if(MixedGenome::nb)	bgnom->copy(*(bc.bgnom));
    */

    if(MixedGenome::nr){
       if(rgnom->length()==bc.rgnom->length())
          rgnom->copy(*(bc.rgnom));
       else{
          delete rgnom;
          rgnom = (GARealGenome*)bc.rgnom->clone();
       }
    }

    if(MixedGenome::ni){
       if(ignom->length()==bc.ignom->length())
          ignom->copy(*(bc.ignom));
       else{
          delete ignom;
          ignom = (GARealGenome*)bc.ignom->clone();
       }
    }

    if(MixedGenome::nd){
       if(dgnom->length()==bc.dgnom->length())
          dgnom->copy(*(bc.dgnom));
       else{
          delete dgnom;
          dgnom = (GARealGenome*)bc.dgnom->clone();
       }
    }

    if(MixedGenome::nb){
       if(bgnom->length()==bc.bgnom->length())
          bgnom->copy(*(bc.bgnom));
       else{
          delete bgnom;
          bgnom = (GARealGenome*)bc.bgnom->clone();
       }
    }

  }
}
//---------------------------------------------------------------
int MixedGenome::equal(const GAGenome& g) const {
  MixedGenome& genome = (MixedGenome&)g;
  int eql=0;
  if(MixedGenome::nr){
       GARealGenome &g1 = (GARealGenome &) (this->realgenome());
       GARealGenome &g2 = (GARealGenome &) (genome.realgenome());
       eql = g1.equal(g2);
  }
  if(MixedGenome::ni){
       GARealGenome &h1 = (GARealGenome &) (this->intgenome());
       GARealGenome &h2 = (GARealGenome &) (genome.intgenome());
       eql = (eql && h1.equal(h2));
  }
  if(MixedGenome::nd){
       GARealGenome &p1 = (GARealGenome &) (this->disgenome());
       GARealGenome &p2 = (GARealGenome &) (genome.disgenome());
       eql = (eql && p1.equal(p2));
  }
  if(MixedGenome::nb){
       GARealGenome &q1 = (GARealGenome &) (this->bingenome());
       GARealGenome &q2 = (GARealGenome &) (genome.bingenome());
       eql = (eql && q1.equal(q2));
  }
  return eql;

}
//---------------------------------------------------------------
int MixedGenome::read(std::istream & is) {
  float sc=0.0;
  if(MixedGenome::nr)	is>>*rgnom;
  if(MixedGenome::ni)	is>>*ignom;
  if(MixedGenome::nd)	is>>*dgnom;
  if(MixedGenome::nb)	is>>*bgnom;
  is>>sc;
  this->score(sc);// set the score
  return is.fail() ? 1 : 0;
}
//---------------------------------------------------------------
int MixedGenome::write(std::ostream & os) const {
  if(MixedGenome::nr)	os<<*rgnom;
  if(MixedGenome::ni)	os<<*ignom;
  if(MixedGenome::nd)	os<<*dgnom;
  if(MixedGenome::nb)	os<<*bgnom;
  os << _score;

  return os.fail() ? 1 : 0;
}
//---------------------------------------------------------------


// These are the default initialization, mutation, and comparator operators for
// this genome class.  They are defined as static functions of the composite
// genome class and they're defaults for the class.  But they can be overridden
// on any instance of the genome.

// The initializer just calls the initializer for each of the genomes that are
// in the composite genome.

// I would have used simply 'Initializer', 'Mutator', etc rather than
// 'CompositeInitializer' but old versions of g++ are brain-dead and don't
// get the encapsulation properly.

void MixedGenome::MixedInitializer(GAGenome & c) {
  MixedGenome & child = (MixedGenome &)c;

  if(MixedGenome::nr)	child.realgenome().initialize();
  if(MixedGenome::ni)	child.intgenome ().initialize();
  if(MixedGenome::nd)	child.disgenome ().initialize();
  if(MixedGenome::nb)	child.bingenome ().initialize();

  child._evaluated = gaFalse;
}
//---------------------------------------------------------------

// The mutator just calls the mutator for each of the component genomes.
int MixedGenome::MixedMutator(GAGenome & c, float pmut) {
  MixedGenome & child = (MixedGenome &)c;
  int nm=0;
  if(MixedGenome::nr)	nm = child.realgenome().mutate(pmut);
  if(MixedGenome::ni)	nm = nm + child.intgenome().mutate(pmut);
  if(MixedGenome::nd)	nm = nm + child.disgenome().mutate(pmut);
  if(MixedGenome::nb)	nm = nm + child.bingenome().mutate(pmut);
  if(nm) child._evaluated = gaFalse;
  return nm;
}

//---------------------------------------------------------------
// The comparator just calls the comparators for each of the component genomes,
// then adds the score using different weight (for example 1 for a real, 10 for
// a discretized real genome, 100 for integer and 1000 for a binary).
// Theese weights are used to maintaint the diversity.
float MixedGenome::MixedComparator(const GAGenome& a, const GAGenome& b) {

  MixedGenome& sis = (MixedGenome &)a;
  MixedGenome& bro = (MixedGenome &)b;
  float ff=0;

  if(MixedGenome::nr) ff = ff + 1    * NormSqrDisComparator(sis.realgenome(), bro.realgenome());
  if(MixedGenome::nd) ff = ff + 10   * NormSqrDisComparator(sis.disgenome (), bro.disgenome ());
  if(MixedGenome::ni) ff = ff + 100  * NormSqrDisComparator(sis.intgenome (), bro.intgenome ());
  if(MixedGenome::nb) ff = ff + 1000 * NormSqrDisComparator(sis.bingenome (), bro.bingenome ());
  return ff;
}

// The comparator just calls the comparators for each of the component genomes,
// then adds the score using different weight (for example 1-real, 10-discretized,
// for the real genome, 100 integer and 100-binary). to maintaint the diversity.
/*float MixedGenome::MixedComparator(const GAGenome& a, const GAGenome& b) {

  MixedGenome& sis = (MixedGenome &)a;
  MixedGenome& bro = (MixedGenome &)b;
  int lr	= sis.realgenome().length();
  int ld	= sis.disgenome().length();
  int li	= sis.intgenome().length();
  int lb	= sis.bingenome().length();
  int len= lr + ld + li+ lb;
  float ff=0;

  if(MixedGenome::nr)	ff = ff + lr * GARealElementComparator (sis.realgenome(), bro.realgenome());
  if(MixedGenome::nd)	ff = ff + ld * GARealElementComparator (sis.disgenome (), bro.disgenome ());
  if(MixedGenome::ni)	ff = ff + li * GARealElementComparator (sis.intgenome (), bro.intgenome ());
  if(MixedGenome::nb)	ff = ff + lb * GARealElementComparator (sis.bingenome (), bro.bingenome ());
  return ff/len;
}*/
//---------------------------------------------------------------

// The crossover operator invokes the crossover for each of the genomes in the
// composite genome.
int
MixedGenome::
MixedCrossover(const GAGenome& a, const GAGenome& b,
		   GAGenome* c, GAGenome* d){
  MixedGenome& mom = (MixedGenome&)a;
  MixedGenome& dad = (MixedGenome&)b;

  int n=0;
  if(c && d){
    MixedGenome& sis = (MixedGenome&)*c;
    MixedGenome& bro = (MixedGenome&)*d;
    if(MixedGenome::nr)
       GARealUniformCrossover  ( mom.realgenome(),  dad.realgenome(),
       		   &sis.realgenome(), &bro.realgenome());
    if(MixedGenome::ni)
       GARealUniformCrossover( mom.intgenome(),  dad.intgenome(),
       			   &sis.intgenome(), &bro.intgenome());
    if(MixedGenome::nd)
       GARealUniformCrossover( mom.disgenome(),  dad.disgenome(),
       			   &sis.disgenome(), &bro.disgenome());
    if(MixedGenome::nb)
       GARealUniformCrossover( mom.bingenome(),  dad.bingenome(),
       			   &sis.bingenome(), &bro.bingenome());
    n = 2;
    sis._evaluated= gaFalse;
    bro._evaluated= gaFalse;
  }
  else if(c){
    MixedGenome& sis = (MixedGenome&)*c;
    if(MixedGenome::nr)
       GARealUniformCrossover (mom.realgenome(), dad.realgenome(), &sis.realgenome(), 0);
    if(MixedGenome::ni)
       GARealUniformCrossover(mom.intgenome(), dad.intgenome (), &sis.intgenome (), 0);
    if(MixedGenome::nd)
       GARealUniformCrossover(mom.disgenome(), dad.disgenome (), &sis.disgenome (), 0);
    if(MixedGenome::nb)
       GARealUniformCrossover(mom.bingenome(), dad.bingenome (), &sis.bingenome (), 0);
    n = 1;
    sis._evaluated= gaFalse;
  }
  else if(d){
    MixedGenome& bro = (MixedGenome&)*d;
    if(MixedGenome::nr)
       GARealUniformCrossover (mom.realgenome(), dad.realgenome(), 0, &bro.realgenome());
    if(MixedGenome::ni)
       GARealUniformCrossover(mom.intgenome(), dad.intgenome (), 0, &bro.intgenome ());
    if(MixedGenome::nd)
       GARealUniformCrossover(mom.disgenome(), dad.disgenome (), 0, &bro.disgenome ());
    if(MixedGenome::nb)
       GARealUniformCrossover(mom.bingenome(), dad.bingenome (), 0, &bro.bingenome ());
    n = 1;
    bro._evaluated= gaFalse;
  }
  return n;
}


