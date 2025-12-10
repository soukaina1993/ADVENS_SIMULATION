// $Header: /nfs/dsi/cvs/galib/ga/GARealGenome.C,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  real.C
  mbwall 11nov95
  Copyright (c) 1995-1996 Massachusetts Institute of Technology
                          all rights reserved

 DESCRIPTION:
   Source file for the real number specialization of the array genome.

 Modified by malick.kane@epfl.ch Nov.2000
---------------------------------------------------------------------------- */

#include "GARealGenome.h"
//#include <string.h>

#define SQR(x) ((x)*(x))

#if !defined(NO_STREAMS) && defined(USE_STDH)
#include <iostream.h>
#elif !defined(NO_STREAMS)
#include <iostream>
#endif

// template <> const char *
// GA1DArrayAlleleGenome<float>::className() const {return "GARealGenome";}
//template <> int
//GA1DArrayAlleleGenome<float>::classID() const {return GAID::FloatGenome;}

// template <> GA1DArrayAlleleGenome<float>::
// GA1DArrayAlleleGenome(unsigned int length, const GAAlleleSet<float> & s,
// 		      GAGenome::Evaluator f, void * u) :
// GA1DArrayGenome<float>(length, f, u){
//   naset = 1;
//   aset = new GAAlleleSet<float>[1];
//   aset[0] = s;

//   initializer(DEFAULT_REAL_INITIALIZER);
//   mutator(DEFAULT_REAL_MUTATOR);
//   comparator(DEFAULT_REAL_COMPARATOR);
//   crossover(DEFAULT_REAL_CROSSOVER);
// }

// template <> GA1DArrayAlleleGenome<float>::
// GA1DArrayAlleleGenome(const GAAlleleSetArray<float> & sa,
// 		      GAGenome::Evaluator f, void * u) :
// GA1DArrayGenome<float>(sa.size(), f, u){
//   naset = sa.size();
//   aset = new GAAlleleSet<float>[naset];
//   for(int i=0; i<naset; i++)
//     aset[i] = sa.set(i);

//   initializer(DEFAULT_REAL_INITIALIZER);
//   mutator(DEFAULT_REAL_MUTATOR);
//   comparator(DEFAULT_REAL_COMPARATOR);
//   crossover(DEFAULT_REAL_CROSSOVER);
// }

// template <>
// GA1DArrayAlleleGenome<float>::~GA1DArrayAlleleGenome(){
//   delete [] aset;
// }



#ifndef NO_STREAMS
// The read specialization takes in each number and stuffs it into the array.
// template <> int
// GA1DArrayAlleleGenome<float>::read(istream & is) {
//   unsigned int i=0;
//   float val;
//   do{
//     is >> val;
//     if(!is.fail()) gene(i++, val);
//   } while(!is.fail() && !is.eof() && i < nx);

//   if(is.eof() && i < nx){
//     GAErr(GA_LOC, className(), "read", gaErrUnexpectedEOF);
//     is.clear(ios::badbit | is.rdstate());
//     return 1;
//   }
//   return 0;
// }

// No need to specialize the write method.
#endif







// We must also specialize the allele set so that the alleles are handled
// properly.  Be sure to handle bounds correctly whether we are discretized
// or continuous.  Handle the case where someone sets stupid bounds that
// might cause an infinite loop for exclusive bounds.
// template <> float
// GAAlleleSet<float>::allele() const {
//   float value = 0.0;
//   if(core->type == GAAllele::ENUMERATED)
//     value = core->a[GARandomInt(0, core->sz-1)];
//   else if(core->type == GAAllele::DISCRETIZED){
//     float n = (core->a[1] - core->a[0]) / core->a[2];
//     int m = (int)n;
//     if(core->lowerb == GAAllele::EXCLUSIVE) m -= 1;
//     if(core->upperb == GAAllele::EXCLUSIVE) m -= 1;
//     value = core->a[0] + GARandomInt(0,(int)m) * core->a[2];
//     if(core->lowerb == GAAllele::EXCLUSIVE) value += core->a[2];
//   }
//   else{
//     if(core->a[0] == core->a[1] &&
//        core->lowerb == GAAllele::EXCLUSIVE &&
//        core->upperb == GAAllele::EXCLUSIVE) {
//       value = core->a[0];
//     }
//     else {
//       do {
// 	value = GARandomFloat(core->a[0], core->a[1]);
//       } while ((core->lowerb == GAAllele::EXCLUSIVE && value == core->a[0]) ||
// 	       (core->upperb == GAAllele::EXCLUSIVE && value == core->a[1]));
//     }
//   }
//   return value;
// }

// If someone asks for a discretized item that is beyond the bounds, give them
// one of the bounds.  If they ask for allele item when there is no
// discretization or enumeration, then error and return lower bound.
template <> float
GAAlleleSet<float>::allele(unsigned int i) const {
  float value = 0.0;
  if(core->type == GAAllele::ENUMERATED)
    value = core->a[i % core->sz];
  else if(core->type == GAAllele::DISCRETIZED){
    float n = (core->a[1] - core->a[0])/core->a[2];
    unsigned int m = (unsigned int)n;            // what about bogus limits?
    if(core->lowerb == GAAllele::EXCLUSIVE) m -= 1;
    if(core->upperb == GAAllele::EXCLUSIVE) m -= 1;
    if(i > m) i = (int)m;
    value = core->a[0] + i*core->a[2];
    if(core->lowerb == GAAllele::EXCLUSIVE) value += core->a[2];
  }
  else{
    GAErr(GA_LOC, "GAAlleleSet", "allele", gaErrNoAlleleIndex);
    value = core->a[0];
  }
  return value;
}







/* ----------------------------------------------------------------------------
   Operator specializations
---------------------------------------------------------------------------- */
// The Gaussian mutator picks a new value based on a Gaussian distribution
// around the current value.  We respect the bounds (if any).
//*** need to figure out a way to make the stdev other than 1.0
int
GARealGaussianMutator(GAGenome& g, float pmut){
  GA1DArrayAlleleGenome<float> &child=
    DYN_CAST(GA1DArrayAlleleGenome<float> &, g);
//  register int n, i;
  int n, i;
  if(pmut <= 0.0) return(0);

  float nMut = pmut * (float)(child.length());
  int length = child.length()-1;
  if(nMut < 1.0){		// we have to do a flip test on each element
    nMut = 0;
    for(i=length; i>=0; i--){
      float value = child.gene(i);
      if(GAFlipCoin(pmut)){
	if(child.alleleset(i).type() == GAAllele::ENUMERATED ||
	   child.alleleset(i).type() == GAAllele::DISCRETIZED)
	  value = child.alleleset(i).allele();
	else if(child.alleleset(i).type() == GAAllele::BOUNDED){
	  value += (float)GAUnitGaussian();
	  value = GAMax(child.alleleset(i).lower(), value);
	  value = GAMin(child.alleleset(i).upper(), value);
	}
	child.gene(i, value);
	nMut++;
      }
    }
  }
  else{				// only mutate the ones we need to
    for(n=0; n<nMut; n++){
      int idx = GARandomInt(0,length);
      float value = child.gene(idx);
      if(child.alleleset(idx).type() == GAAllele::ENUMERATED ||
	 child.alleleset(idx).type() == GAAllele::DISCRETIZED)
	value = child.alleleset(idx).allele();
      else if(child.alleleset(idx).type() == GAAllele::BOUNDED){
	value += (float)GAUnitGaussian();
	value = GAMax(child.alleleset(idx).lower(), value);
	value = GAMin(child.alleleset(idx).upper(), value);
      }
      child.gene(idx, value);
    }
  }
  return((int)nMut);
}


// Arithmetic crossover generates a new value that is the average of the parent
// values.  Note that this means both children in a sexual crossover will be
// identical.  If parents are not the same length, the extra elements are not
// set!  You might want to add some noise to this so that both children are not
// the same...
int
GARealArithmeticCrossover(const GAGenome& p1, const GAGenome& p2,
			  GAGenome* c1, GAGenome* c2) {
  const GA1DArrayGenome<float> &mom=
    DYN_CAST(const GA1DArrayGenome<float> &, p1);
  const GA1DArrayGenome<float> &dad=
    DYN_CAST(const GA1DArrayGenome<float> &, p2);

  int n=0;

  if(c1 && c2){
    GA1DArrayGenome<float> &sis=DYN_CAST(GA1DArrayGenome<float> &, *c1);
    GA1DArrayGenome<float> &bro=DYN_CAST(GA1DArrayGenome<float> &, *c2);

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      sis.gene(i, 0.5f * (mom.gene(i) + dad.gene(i)));
      bro.gene(i, 0.5f * (mom.gene(i) + dad.gene(i)));
    }
    n = 2;
  }
  else if(c1 || c2){
    GA1DArrayGenome<float> &sis = (c1 ?
				   DYN_CAST(GA1DArrayGenome<float> &, *c1) :
				   DYN_CAST(GA1DArrayGenome<float> &, *c2));

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      sis.gene(i, 0.5f * (mom.gene(i) + dad.gene(i)));
    }
    n = 1;
  }

  return n;
}


// Blend crossover generates a new value based on the interval between parents.
// We generate a uniform distribution based on the distance between parent
// values, then choose the child value based upon that distribution.
// modified malick.kane@epfl.ch, leni (Nov2000)
int
GARealBlendCrossover(const GAGenome& p1, const GAGenome& p2,
		     GAGenome* c1, GAGenome* c2) {
  /*const GA1DArrayGenome<float> &mom=
    DYN_CAST(const GA1DArrayGenome<float> &, p1);
  const GA1DArrayGenome<float> &dad=
    DYN_CAST(const GA1DArrayGenome<float> &, p2);*/

  GARealGenome &mom=(GARealGenome &)p1;
  GARealGenome &dad=(GARealGenome &)p2;
  int n=0;

  if(c1 && c2){
    //GA1DArrayGenome<float> &sis=DYN_CAST(GA1DArrayGenome<float> &, *c1);
    //GA1DArrayGenome<float> &bro=DYN_CAST(GA1DArrayGenome<float> &, *c2);

    // FRO : avant onavait :
    // GARealGenome &sis=(GARealGenome &)*c1;
    // Ce qui devait forcement planter puisque on caste un pointeur en objet pointe
    // ceci parait plus logique, la fonction n'est pas appellee par le programme de toute facon
    GARealGenome &sis=(GARealGenome &)*c1;
    GARealGenome &bro=(GARealGenome &)*c2;

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      float dist = 0;
      if(mom.gene(i) > dad.gene(i))
	dist = mom.gene(i) - dad.gene(i);
      else
	dist = dad.gene(i) - mom.gene(i);
      //float lo = (GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist;
      //float hi = (GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist;
      float lo = (float)GAMax(((GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist),mom.alleleset(i).lower());
      float hi = (float)GAMin(((GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist),mom.alleleset(i).upper());
      sis.gene(i, GARandomFloat(lo, hi));
      bro.gene(i, GARandomFloat(lo, hi));
    }
    n = 2;
  }
  else if(c1 || c2){
    /*GA1DArrayGenome<float> &sis = (c1 ?
				   DYN_CAST(GA1DArrayGenome<float> &, *c1) :
				   DYN_CAST(GA1DArrayGenome<float> &, *c2));*/
    GARealGenome &sis = (c1 ? (GARealGenome &)*c1: (GARealGenome &)*c2);

    int len = GAMax(mom.length(), dad.length());
    for(int i=0; i<len; i++) {
      float dist = 0;
      if(mom.gene(i) > dad.gene(i))
	dist = mom.gene(i) - dad.gene(i);
      else
	dist = dad.gene(i) - mom.gene(i);
      //float lo = (GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist;
      //float hi = (GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist;
      float lo = (float)GAMax(((GAMin(mom.gene(i), dad.gene(i))) - 0.5*dist),mom.alleleset(i).lower());
      float hi = (float)GAMin(((GAMax(mom.gene(i), dad.gene(i))) + 0.5*dist),mom.alleleset(i).upper());
      sis.gene(i, GARandomFloat(lo, hi));
    }
    n = 1;
  }

  return n;
}


//---------------------------------------------------------------
// simple comparator which calculated the euclidean distance
// between the two given individuals

float EucDisComparator(const GAGenome& g1, const GAGenome& g2) {
     GARealGenome & a = ( GARealGenome &) g1;
     GARealGenome & b = ( GARealGenome &) g2;
     float sum=0.0;
    // register int i;
    // register int len = a.length();    
     int i;
     int len = a.length();
     static float *range=NULL;

     if(!range) {
       range = new float[len];
       for(i=0;i<len;++i)
         range[i]=(a.alleleset(i).weight() )/(a.alleleset(i).upper() - a.alleleset(i).lower());
     }
     for(i=0;i<len;++i) {
       sum += SQR(range[i]*(a.gene(i)-b.gene(i)));
     }

     return sqrt(sum);
}


//---------------------------------------------------------------
// simple comparator which calculated the normalised sqr distance
// between the two given individuals
// It returns a number in the interval [0,1] where 0 means that
// the two genomes are identical (zero diversity) and 1 means they are
// completely different (maximum diversity).

float NormSqrDisComparator (const GAGenome& g1, const GAGenome& g2) {
     GARealGenome & a = ( GARealGenome &) g1;
     GARealGenome & b = ( GARealGenome &) g2;
     float sum=0.0;
  //   register int i;
  //   register int len = a.length();
  //   register int cpt = a.length();
     int i;
     int len = a.length();
     int cpt = a.length(); 
     static float *range=NULL;

     if(!range) {
       range = new float[len];
       for(i=0;i<len;++i){
         range[i]=(a.alleleset(i).weight() )/(a.alleleset(i).upper() - a.alleleset(i).lower());
         if (a.alleleset(i).weight()==0)
            cpt--;
       }
     }
     for(i=0;i<len;++i) {
       sum += SQR(range[i]*(a.gene(i)-b.gene(i)));
     }

     return sum/cpt;
}

// The comparator returns a number in the interval [0,1] where 0 means that
// the two genomes are identical (zero diversity) and 1 means they are
// completely different (maximum diversity).
/*
float ExpDisComparator (const GAGenome& g1, const GAGenome& g2) {
  GARealGenome& a = (GARealGenome &)g1;
  GARealGenome& b = (GARealGenome &)g2;
  return exp( -(a.gene(0) - b.gene(0)) * (a.gene(0) - b.gene(0)) / 100.0);
}
*/


// These must be included _after_ the instantiations because some compilers get
// all wigged out about the declaration/specialization order.  Note that some
// compilers require a syntax different than others when forcing the
// instantiation (i.e. GNU wants the 'template class', borland does not).
#ifndef USE_AUTO_INST
#include "GAAllele.cpp"
#include "GA1DArrayGenome.cpp"

GA_CPPUSE(ios)
GA_CPPUSE(istream)

#if defined(__BORLANDC__)
GAAlleleSet<float>;
GAAlleleSetCore<float>;
GAAlleleSetArray<float>;

GAArray<float>;
GA1DArrayGenome<float>;
GA1DArrayAlleleGenome<float>;
#else
template class GAAlleleSet<float>;
template class GAAlleleSetCore<float>;
template class GAAlleleSetArray<float>;

template class GAArray<float>;
template class GA1DArrayGenome<float>;
template class GA1DArrayAlleleGenome<float>;
#endif
#endif
