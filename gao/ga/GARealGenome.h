// $Header: /nfs/dsi/cvs/galib/ga/GARealGenome.h,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  real.h
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This header defines the specialization of the array genome of type float
for the real number genome.

Modified by malick.kane@epfl.ch
---------------------------------------------------------------------------- */
#ifndef _ga_real_h_
#define _ga_real_h_

#include "GAAllele.h"
#include "GA1DArrayGenome.h"

typedef GAAlleleSet<float> GARealAlleleSet;
typedef GAAlleleSetArray<float> GARealAlleleSetArray;
typedef GA1DArrayAlleleGenome<float> GARealGenome;

inline void GARealUniformInitializer(GAGenome& g){
  GA1DArrayAlleleGenome<float>::UniformInitializer(g);
}
inline void GARealOrderedInitializer(GAGenome& g){
  GA1DArrayAlleleGenome<float>::OrderedInitializer(g);
}

inline int GARealUniformMutator(GAGenome& g, float pmut){
  return GA1DArrayAlleleGenome<float>::FlipMutator(g, pmut);
}
inline int GARealSwapMutator(GAGenome& g, float pmut){
  return GA1DArrayGenome<float>::SwapMutator(g, pmut);
}
int GARealGaussianMutator(GAGenome &, float);

inline int GARealUniformCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::UniformCrossover(a,b,c,d);
}
inline int GARealEvenOddCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::EvenOddCrossover(a,b,c,d);
}
inline int GARealOnePointCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::OnePointCrossover(a,b,c,d);
}
inline int GARealTwoPointCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::TwoPointCrossover(a,b,c,d);
}
inline int GARealPartialMatchCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::PartialMatchCrossover(a,b,c,d);
}
inline int GARealOrderCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::OrderCrossover(a,b,c,d);
}
inline int GARealCycleCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d) {
  return GA1DArrayGenome<float>::CycleCrossover(a,b,c,d);
}

inline float GARealElementComparator(const GAGenome& a, const GAGenome& b) {
  return GA1DArrayGenome<float>::ElementComparator(a,b);
}

int GARealArithmeticCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d);
int GARealBlendCrossover(const GAGenome& a, const GAGenome& b,
				  GAGenome* c, GAGenome* d);
				  


//Rmq:
// The comparator is supposed to return a number that indicates how similar
// two genomes are. 
// The default comparator used here is the simple  called ElementComparator, 
// so here il just compare elements and return a number that
// indicates how many elements match.  



//---------------------------------------------------------------
// simple comparator which calculated the euclidean distance
// between the two givven individuals
float EucDisComparator(const GAGenome& g1, const GAGenome& g2);

//---------------------------------------------------------------
// simple comparator which calculated the normalized sqr distance
// between the two given individuals
// It returns a number in the interval [0,1] where 0 means that
// the two genomes are identical (zero diversity) and 1 means they are 
// completely different (maximum diversity).
float NormSqrDisComparator(const GAGenome& g1, const GAGenome& g2);


// The comparator returns also a number in the interval [0,1] where 0 means that
// the two genomes are identical (zero diversity) and 1 means they are 
// completely different (maximum diversity).
//float ExpDisComparator(const GAGenome& g1, const GAGenome& g2);



#endif
