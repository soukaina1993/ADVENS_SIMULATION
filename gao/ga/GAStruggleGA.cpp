/* ----------------------------------------------------------------------------
  StruggleHGA.C
  mbwall jan96
  Copyright (c) 1995-1996 Massachusetts Institute of Technology
                          all rights reserved

   Souce file for the struggle genetic algorithm object.  This algorithm was
developed by Thomas Grueninger while he studied at MIT.  It is very similar
to deterministic crowding.
---------------------------------------------------------------------------- */
#include "garandom.h"
#include "GAStruggleGA.h"
#include "ga.h"
#include <iostream>
#undef HOOPS

#define ONE_CHILD_ONLY

GAGenome * SelectFitnessProportional(GAPopulation *pop,GAGenome *mom,int minmax);
GAGenome * SelectUniform(GAPopulation *pop,GAGenome *mom,float pLocal,int minmax);

GAGenome * SelectUniform(GAPopulation *pop,GAGenome *mom,float pLocal,int minmax)
{
  GAGenome *oldInd =   ((DistanceScaling&)pop->scaling()).getInd();
  GAGenome *dad;

  if(!GAFlipCoin(pLocal)) {
    dad = &(pop->select());
  }
  else {
    ((DistanceScaling&)pop->scaling()).assignInd(mom);
    pop->scale(gaTrue);

  if (minmax == StruggleHGA::MINIMIZE) {
      int cc=1;
      while((pop->best(cc,GAPopulation::SCALED).fitness() == 0.0) && (cc < pop->size()))
	cc++;
      dad = &(pop->best(cc,GAPopulation::SCALED));
    }
    else {
      int cc=1;
      while((pop->worst(cc,GAPopulation::SCALED).fitness() == 0.0) && (cc < pop->size()))
	cc++;
      dad = &(pop->worst(1,GAPopulation::SCALED));
    }
    ((DistanceScaling&)pop->scaling()).assignInd(oldInd);
  }

  return(dad);
}



GAGenome * SelectFitnessProportional(GAPopulation *pop,GAGenome *mom,int minmax)
{
  GAGenome *oldInd =   ((DistanceScaling&)pop->scaling()).getInd();
  GAGenome *dad;
  static int n=-1;
  static  float *p=NULL;
  float psum=0.0;
  float cutoff;
  int i;

  if(pop->size() >n) {
    if(p)
	delete [] p;
    n = pop->size();
    p = new float[n];
  }

  ((DistanceScaling&)pop->scaling()).assignInd(mom);
  pop->scale(gaTrue);

  p[0] = 0;
  float xx;
  if (minmax == StruggleHGA::MINIMIZE) {
    for(i=1 ;i<n;i++) {
      xx = pop->best(i,GAPopulation::SCALED).fitness();
      if(xx > 0.0)
    	p[i] = (float)(1.0/(1 + pop->best(i,GAPopulation::SCALED).fitness()));
      else
        p[i] = 0.0;
      psum += p[i];
      p[i] += p[i-1];
    }
  } else {
    for(i=1 ;i<n;i++) {
      xx = pop->worst(i,GAPopulation::SCALED).fitness();
      if(xx > 0.0)
    	p[i] = (float)(1.0/(1 + pop->worst(i,GAPopulation::SCALED).fitness()));
      else
        p[i] = 0.0;
      psum += p[i];
      p[i] += p[i-1];
    }
  }
  for(i=1 ;i<n;i++) {
    p[i] /= psum;
  }

  //  for(i=0;i<n;++i)
  //    cerr << p[i] <<"  fitness = "<<pop->best(i,GAPopulation::SCALED).fitness()<<"\n";
  //  cerr << "Hit a key\n";
  //  getchar();

  cutoff = GARandomFloat();
  i=1;
  while(i < pop->size() && cutoff > p[i])
	i++;

  if (minmax == StruggleHGA::MINIMIZE) {
    dad = &( pop->best(i,GAPopulation::SCALED ));
  }
  else{
    dad = &( pop->worst(i,GAPopulation::SCALED ));
  }

  ((DistanceScaling&)pop->scaling()).assignInd(oldInd);
  return(dad);
}




void
DistanceScaling::evaluate(const GAPopulation & p) {
  for(int i=0; i<p.size(); i++) {
    p.individual(i).fitness( p.individual(i).compare(*genome) );
  }
}


GAParameterList&
StruggleHGA::registerDefaultParameters(GAParameterList& p) {
  GAGeneticAlgorithm::registerDefaultParameters(p);

  int ival = 1;
  p.add(gaNnReplacement, gaSNnReplacement, GAParameter::INT, &ival);
  p.add(gaNpReplacement, gaSNpReplacement, GAParameter::FLOAT, &gaDefPRepl);

  p.set(gaNscoreFrequency, gaDefScoreFrequency2);

  return p;
}

StruggleHGA::StruggleHGA(const GAGenome& g) : GAGeneticAlgorithm(g) {
  DistanceScaling sc;
  pop->scaling(sc);
  child = g.clone();
  child->geneticAlgorithm(*this);

  ((DistanceScaling&)pop->scaling()).assignInd(child);
  GAUniformSelector select(GAPopulation::RAW);
  pop->selector(select);
  _SelectMethod = 0;		// default to use MakeNewPop method
  pLocal=0.2;
  drawpop = DrawPopulationDefault;
}

StruggleHGA::StruggleHGA(const GAPopulation& p) : GAGeneticAlgorithm(p) {
  DistanceScaling sc;
  pop->scaling(sc);
  child = p.individual(1).clone();
  child->geneticAlgorithm(*this);
  ((DistanceScaling&)pop->scaling()).assignInd(child);
  GAUniformSelector select(GAPopulation::RAW);
  pop->selector(select);
  _SelectMethod = 0;
  drawpop = DrawPopulationDefault;
  pLocal=0.2;

}

void
StruggleHGA::initialize(unsigned int seed) {
  GARandomSeed(seed);
  pop->initialize();
  std::ofstream fa("popinit.dat",std::ios::out| std::ios::trunc);
  fa<<"Initial population\n";
  fa<<*pop<<"\n";
  fa.close();
  pop->evaluate(gaTrue);
  stats.reset(*pop);
  if(!scross) GAErr(GA_LOC, className(), "initialize", gaErrNoSexualMating);
}

GAPopulation & StruggleHGA::MakeChildPop(GAPopulation &cpop)
{
  GAGenome *mom, *dad;
  int ps;

  //ps = pcross * pop->size();
  //ps=5;// fix debug size of pop
  ps = pop->size();
  cpop.userData(pop->userData());  // need in evaluation of cpop
  cpop.evaluator(pop->evaluator());

  for (int i=0; i<ps; i++) {
    mom = &(pop->select());
    dad = &(pop->select());
    stats.numsel += 2;
    stats.numcro += (*scross)(*mom, *dad, child, 0);
    stats.nummut += child->mutate(pMutation());
    //child->setEvaluated(gaFalse);
    child->evaluate(gaFalse);
    cpop.add((*child));  // add a clone of child
  }
  // now have all the potential children
  cpop.evaluate(gaTrue);  //evaluate
  stats.numeval += cpop.size();
  // now have evaluated cpop
  return(cpop);
}


void
StruggleHGA::step() {

  GAPopulation cpop;
  int i;
  GAGenome *mom, *dad;

    for ( i=0; i<pop->size(); i++) {
      mom = &(pop->select());
	  switch(_SelectMethod) {
		case UNIFORM_WITH_CORRECTION:
		  dad = SelectUniform(pop,mom,pLocal,minmax);
		  break;
		case PROPORTIONAL:
		  dad = SelectFitnessProportional(pop,mom,minmax);
		  break;
		default:
		  dad = SelectUniform(pop,mom,pLocal,minmax);
		  break;
	  }

      stats.numsel += 2;
      stats.numcro += (*scross)(*mom, *dad, child, 0);
      stats.nummut += child->mutate(pMutation());
      child->evaluate(gaTrue);	//
      pop->scale(gaTrue);   // recalc the distance etc.
      if (minmax == MINIMIZE) {
	if (child->score() <= pop->best(0,GAPopulation::SCALED).score()) {
	  pop->best(0,GAPopulation::SCALED).copy(*child);
	  stats.numrep += 1;
	}

      }
      else {
	if (child->score() >= pop->worst(0,GAPopulation::SCALED).score()) {
	  pop->worst(0,GAPopulation::SCALED).copy(*child);
	  stats.numrep += 1;
	}
      }
  }

  pop->evaluate(gaTrue);
  stats.update(*pop);
  (*drawpop)(*pop);
  std::cout.flush();
}

int DrawPopulationDefault(GAPopulation &pop)
{
GAPopulation p=pop;
return(0);
}
