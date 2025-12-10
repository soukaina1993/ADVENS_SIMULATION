/****************************************************************************/
#include<fstream>
#include<string.h>
#include"gaiodata.h"
#define SQR(x) ((x)*(x))


/****************************************************************************/
GAParameterList * gaParams::registered (const char * ganame)
{
   if (Fields==NULL)
      return NULL;
   else{
      params->set(gaNpopulationSize, zpop  (ganame));	// population size
      params->set(gaNnGenerations  , ngen  (ganame)); // number of generations
      params->set(gaNpCrossover    , pcros (ganame));	// probability of crossover
      params->set(gaNpMutation     , pmut  (ganame));	// probability of mutation
      params->set(gaNnReplacement  , nrepl (ganame));	// how much of pop to replace each gen
      params->set(gaNnMigration    , nmig  (ganame));	// how much of pop to migrate each gen
      params->set(gaNnPopulations  , npop  (ganame));	// number of populations in parallele
      params->set(gaNscoreFrequency, sfreq (ganame));	// how often to record scores
      params->set(gaNflushFrequency, ffreq (ganame)); // how often to dump scores to file
      params->set(gaNselectScores, (int)GAStatistics::AllScores);
      params->set(gaNscoreFilename, "bog.dat");
  }
  return params;
}
/****************************************************************************/
gaData & operator<<(gaData &Data, TSequence &Seq){
   if (!Data.ofile->isopened())
         Data.ofile->open();
   if (!((strcmp(Data.ofilename(), "--") == 0) || Data.ofile==NULL)){
      Data.ofile->FileID() << Seq.noID <<"\t";
      Data.ofile->FileID() << Seq.Value<<"\t";
      Data.ofile->FileID() << Seq;
   }
   return Data;
}
/****************************************************************************/
gaData::gaData (TInputFile & inFile, TOutputFile& outFile):ioData(inFile, outFile) {
   ConAlles	= new GARealAlleleSetArray;
   cpnt	= NULL;

   iGnom 	= 0;
   nGnom		= 10;
   memo = new Array<TSequence*> (nGnom);
   for (unsigned int i=1; i<nGnom; i++)
      memo->Set(i, NULL);
}
//-------------------------------------------------------------
gaData::gaData (TInputFile & inFile):ioData(inFile){
   ConAlles	= new GARealAlleleSetArray;
   cpnt		= NULL;

   iGnom 	= 0;
   nGnom		= 10;
   memo 		= new Array<TSequence*> (nGnom);
   for (unsigned int i=0; i<nGnom; i++)
      memo->Set(i, NULL);
}
//-------------------------------------------------------------
gaData::gaData (const char * ifname, const char * ofname, unsigned int nTln, unsigned int nTcol):
   ioData(ifname, ofname, nTln, nTcol)
{
   ConAlles  = new GARealAlleleSetArray;
   cpnt	= NULL;

   iGnom 	= 0;
   nGnom		= 10;
   memo = new Array<TSequence*> (nGnom);

   for (unsigned int i=0; i<nGnom; i++)
      memo->Set(i, NULL);

}
//-------------------------------------------------------------
unsigned int gaData::nGenome (unsigned int gn)
{
   for (unsigned int i=0; i<nGnom; i++)
      if (!(memo->Get(i)==NULL)) delete memo->Get(i);
   delete memo;

   iGnom 	= 0;
   nGnom		= gn;
   memo 		= new Array<TSequence*> (nGnom);

   for (unsigned int i=0; i<nGnom; i++)
      memo->Set(i, NULL);

   return nGnom;
}
//-------------------------------------------------------------
TSequence * gaData::memorize (void)
{

   unsigned int Curr=iGnom;
   if (Curr<nGnom){
      memo->Set(Curr, genome().clone());
      iGnom++;

   }else return NULL;
   return memo->Get(Curr);
}
//-------------------------------------------------------------
TSequence::Evaluator gaData::evaluator(const TSequence::Evaluator f) {
   cpnt	= NULL;
   return eval = f;
}
//-------------------------------------------------------------
BasicES * gaData::evaluator(BasicES & thecpnt) 		 {
   eval	= (TSequence::Evaluator)0;
   return cpnt = &thecpnt;
}
//-------------------------------------------------------------
float gaData::evaluate (){
   if (!(eval==(TSequence::Evaluator)0)){
      genome().evaluator(eval);
      score = genome().evaluate();
   }else{
      if (!(cpnt==NULL))
         score	= cpnt->objective(this->genome());
      else return score = 0;
   }
   genome().Value = score;
   return score;
}
//-------------------------------------------------------------
gaData & gaData::store (GAGenome&  g){
      GARealGenome& chrom = (GARealGenome&)g;
      int j=0;
      for (unsigned int i=0; i<this->count(); i++){
          if (this->isopt(i)!=0){
              this->genome(i, chrom.gene(j));
              j++;
          }
      }
      //genome().Value=chrom.score();
      return *this;
}
//-------------------------------------------------------------
void gaData::title (void) {
      if (!ofile->isopened())
            ofile->open();
      ofile->FileID() <<"ngen"<<"\t";
      ofile->FileID() <<"score"<<"\t";
      ofile->titleline(*Records, "chrom");
}
//-------------------------------------------------------------
GARealAlleleSetArray * gaData::registered ()
{
   if (Fields==NULL)
       return NULL;
   else{
      int j=0;
      for (unsigned int i=0; i<count(); i++){
          if (isopt(i)!=0){
              ConAlles->add(lower(i),upper(i),
              		GAAllele::INCLUSIVE, GAAllele::INCLUSIVE);
              ConAlles->weight(j, gaData::weight(i));
              j++;
       	 }
      }
   }
   return ConAlles;
 }
//-------------------------------------------------------------
//historic of the component for the specified genome
bool
gaData::historic (TOutputFile & ofile, unsigned int noInd){
   unsigned int Size = memo->size();
   if (noInd<Size){

	if (!(cpnt==NULL) && !(memo->Get(noInd)==NULL)){
	    cpnt->variable(*memo->Get(noInd));
	    cpnt->performance();
	    cpnt->write(ofile);

	}
   }else return false;

   return true;
}
//-------------------------------------------------------------
//historic of the component for all the memorized genomes
bool
gaData::historic (TOutputFile & ofile){
   bool bb;
   unsigned int Size = memo->size();
   for (unsigned int i=0; i<Size; i++)
         bb= historic(ofile, i);

   return bb;
}
//-------------------------------------------------------------
//sensitivity analysis around of the specified gene
bool
gaData::sensitivity (TOutputFile & ofile, unsigned int noGene){
   unsigned int Size = genome().size();
   if (noGene < Size){
        float mini= lower (noGene);
        float maxi= upper (noGene);

        unsigned int Nb = 10;
        float dX = (maxi-mini) / Nb;
        float alleleX = mini;

        for (unsigned int k=0; k< Nb+1; k++){
	    alleleX = mini + k * dX;
	    genome().Set (noGene, alleleX);

    	    if (!(eval==(TSequence::Evaluator)0)){
		genome().evaluator(eval);
		genome().evaluate();
		ofile.FileID()<<genome().noID<<"\t";
		ofile.FileID()<<genome().Value<<"\t";
		ofile.FileID()<<genome();
	    }else{
		if (!(cpnt==NULL)){
		    cpnt->variable(this->genome());
    		    cpnt->performance();
    		    cpnt->write(ofile);
	    	}else return false;
        	    }
        }
   }else return false;
   return true;
}

/****************************************************************************/





/****************************************************************************/
gaMixedData::gaMixedData  (TInputFile & inFile, TOutputFile& outFile):
			gaData(inFile, outFile), MixedAlleleSet() {
   create();
}
//-------------------------------------------------
gaMixedData::gaMixedData  (TInputFile & inFile):gaData(inFile), MixedAlleleSet() {
   create();
}
//-------------------------------------------------
gaMixedData::gaMixedData  (const char * ifname, const char * ofname,
	    		unsigned int nTln, unsigned int nTcol):
gaData(ifname, ofname, nTln, nTcol), MixedAlleleSet() {
   create();
}
//-------------------------------------------------
gaData & gaMixedData::store (GAGenome&  g)
{
   MixedGenome& chrom = (MixedGenome&)g;
   if(nContinious())	gaMixedData::RealStore (chrom);
   if(nDiscretized())	gaMixedData::DisStore  (chrom);
   if(nInteger())		gaMixedData::IntStore  (chrom);
   if(nBinary())		gaMixedData::BinStore  (chrom);
   return *this;
}
//-------------------------------------------------
unsigned int gaMixedData::RealStore (GAGenome&  g){
      MixedGenome& chrom = (MixedGenome&)g;
      dynstring st="Continious";
      int j=0;
      for (unsigned int i=0; i<this->count(); i++){
          if (isopt(i) != 0 && st==*record(i+1).Title){
              this->genome(i, chrom.realgenome().gene(j));
              j++;
          }
      }
      return j;
}
//-------------------------------------------------
unsigned int gaMixedData::IntStore (GAGenome&  g){
      MixedGenome& chrom = (MixedGenome&)g;
      dynstring st="Integer";
      int j=0;
      for (unsigned int i=0; i<this->count(); i++){
          if (isopt(i) != 0 && st==*record(i).Title){
              this->genome(i, chrom.intgenome().gene(j));
              j++;
          }
      }
      return j;
}
//-------------------------------------------------
unsigned int gaMixedData::DisStore (GAGenome&  g){
      MixedGenome& chrom = (MixedGenome&)g;
      dynstring st="Discretized";
      int j=0;
      for (unsigned int i=0; i<this->count(); i++){
          if (isopt(i) != 0 && st==*record(i).Title){
              this->genome(i, chrom.disgenome().gene(j));
              j++;
          }
      }
      return j;
}
//-------------------------------------------------
unsigned int gaMixedData::BinStore (GAGenome&  g){
      MixedGenome& chrom = (MixedGenome&)g;
      dynstring st="Binary";
      int j=0;
      for (unsigned int i=0; i<this->count(); i++){
          if (isopt(i) != 0 && st==*record(i).Title){
              this->genome(i, chrom.bingenome().gene(j));
              j++;
          }
      }
      return j;
}
//-------------------------------------------------
GARealAlleleSetArray * gaMixedData::registered ()
{
   if(nContinious())	gaMixedData::RealRegister ();
   if(nDiscretized())	gaMixedData::DisRegister  ();
   if(nInteger())		gaMixedData::IntRegister  ();
   if(nBinary())		gaMixedData::BinRegister  ();
   return realals;
}
//-------------------------------------------------
unsigned int gaMixedData::RealRegister ()
{
   dynstring st="Continious";
   int j=0;
   for (unsigned int i=0; i<count(); i++){
       if (isopt(i) != 0 && st==*record(i+1).Title){
           realals->add(lower(i),upper(i),
              	GAAllele::INCLUSIVE, GAAllele::INCLUSIVE);
           realals->weight(j, gaData::weight(i));
           j++;
      }
   }
   return ncon = j;
}
//-------------------------------------------------
unsigned int gaMixedData::DisRegister ()
{
   dynstring st="Discretized";
   int j=0;
   for (unsigned int i=0; i<count(); i++){
       if (isopt(i) != 0 && st==*record(i+1).Title){
           disals->add(lower(i),upper(i), genome(i)-lower(i),
              	GAAllele::INCLUSIVE, GAAllele::EXCLUSIVE);
           disals->weight(j, gaData::weight(i));
           j++;
      }
   }
   return ndis = j;
 }
//-------------------------------------------------
unsigned int gaMixedData::IntRegister ()
{
   dynstring st="Integer";
   int j=0;
   for (unsigned int i=0; i<count(); i++){
       if (isopt(i) != 0 && st==*record(i+1).Title){
           intals->add(lower(i),upper(i), 1.0,
              	GAAllele::INCLUSIVE, GAAllele::INCLUSIVE);
           intals->weight(j, gaData::weight(i));
           j++;
      }
   }
   return nint = j;
 }
 //-------------------------------------------------
unsigned int gaMixedData::BinRegister ()
{
   dynstring st="Binary";
   GARealAlleleSet bin_alles;
   bin_alles.add(0);
   bin_alles.add(1);
   int j=0;
   for (unsigned int i=0; i<count(); i++){
       if (isopt(i) != 0 && st==*record(i+1).Title){
           binals->add(bin_alles);
           binals->weight(j, gaData::weight(i));
           j++;
      }
   }
   return nbin = j;
 }
//-------------------------------------------------
unsigned int gaMixedData::howoften (char * sname) {
   int cpt=0;
   dynstring st=sname;
    cout<<st<<endl;

   for (unsigned int i=0; i<3; i++){
      if (isopt(i)!=0 && st==*record(i+1).Title) // TODO changement ici en i+1 attention!
         cpt++;
   }
   return cpt;
}
//-------------------------------------------------
unsigned int gaMixedData::nUsed() {

   ncon	= howoften("Continious");
   nint	= howoften("Integer");
   ndis	= howoften("Discretized");
   nbin	= howoften("Binary");
   return ncon + ndis + nint + nbin;
}
//-------------------------------------------------
void gaMixedData::create (){

   nUsed();

   if(nContinious()){
      realals 	= ConAlles;
      scon	= new TSequence(nContinious());
      sequence (*scon, "Continious");
   }
   if(nInteger()){
      intals 	= new GARealAlleleSetArray;
      sint	= new TSequence(nInteger());
      sequence (*sint, "Integer");
   }
   if(nDiscretized()){
      disals = new GARealAlleleSetArray;
      sdis	   = new TSequence(nDiscretized());
      sequence (*sdis, "Discretized");
   }
   if(nBinary()){
      binals 	= new GARealAlleleSetArray;
      sbin	= new TSequence(nBinary());
      sequence (*sbin, "Binary");
   }
   //gaMixedData::registered ();

   ifile->firstcolum();
}
//-------------------------------------------------
void gaMixedData::distroy (){

   if(nContinious()) 	{delete scon;}
   if(nInteger()) 	{delete intals; delete sint;}
   if(nDiscretized()) 	{delete disals; delete sdis;}
   if(nBinary()) 		{delete binals; delete sbin;}
}
 //-------------------------------------------------
TSequence &  gaMixedData::sequence (TSequence &seq, char * sname)
{
   unsigned int cpt=0;
   dynstring st=sname;
   for (unsigned int i=0; i<count(); i++){
      if (isopt(i)!=0 && st==*record(i+1).Title){
         seq.Set(cpt, genome(i));
         cpt++;
      }
   }
   return seq;
}
//-------------------------------------------------------------
//sensitivity analysis around of the specified gene
bool
gaMixedData::sensitivity (TOutputFile & ofile, unsigned int noGene){

   dynstring st="Continious";
   unsigned int Size = genome().size();
   if (noGene < Size && st==*record(noGene).Title){
        float mini= lower (noGene);
        float maxi= upper (noGene);

        unsigned int Nb = 10;
        float dX = (maxi-mini) / Nb;
        float alleleX = mini;

        for (unsigned int k=0; k< Nb+1; k++){
	    alleleX = mini + k * dX;
	    genome().Set (noGene, alleleX);

    	    if (!(eval==(TSequence::Evaluator)0)){
		genome().evaluator(eval);
		genome().evaluate();
		ofile.FileID()<<genome().noID<<"\t";
		ofile.FileID()<<genome().Value<<"\t";
		ofile.FileID()<<genome();
	    }else{
		if (!(cpnt==NULL)){
		    cpnt->variable(this->genome());
    		    cpnt->performance();
    		    cpnt->write(ofile);
	    	}else return false;
	    }
        }
   }else return false;
   return true;
}

//-------------------------------------------------
void gaMixedData::title (void) {
   if (!ofile->isopened())
         ofile->open();
   ofile->FileID() <<"ngen"<<"\t";
   ofile->FileID() <<"score"<<"\t";
   //ofile->titleline(ifile->firstcolum());
   ofile->titleline(ifile->title());
}
//-------------------------------------------------
MixedAlleleSet * gaMixedData::mixedalleles (){
   MixedAlleleSet * mxals = new MixedAlleleSet;
   mxals->realals	= this->realals;
   mxals->intals	= this->intals;
   mxals->disals	= this->disals;
   mxals->binals	= this->binals;
   mxals->ncon	= this->ncon;
   mxals->ndis	= this->ndis;
   mxals->nint	= this->nint;
   mxals->nbin	= this->nbin;
   return mxals;
}




/****************************************************************************/

/****************************************************************************/

void gaMoData::create ()
{
   PoffAr = new Array<TSequence*> (nPoff);
   for (unsigned int i=0; i<nPoff; i++)
      PoffAr->Set(i, NULL);

   VdecAr = new Array<TSequence*> (nPoff);
   for (unsigned int i=0; i<nPoff; i++)
      VdecAr->Set(i, NULL);

   tgets = new Array<TTarget*> (nObj);
   for (unsigned int i=0; i<nObj; i++)
      tgets->Set(i, NULL);

   ifile->firstcolum();
}
//-------------------------------------------------
void gaMoData::distroy ()
{
   for (unsigned int i=0; i<nPoff; i++)
      if (!(VdecAr->Get(i)==NULL)) delete VdecAr->Get(i);
   delete VdecAr;
   for (unsigned int i=0; i<nPoff; i++)
      if (!(PoffAr->Get(i)==NULL)) delete PoffAr->Get(i);
   delete PoffAr;
   for (unsigned int i=0; i<nObj; i++)
      if (!(tgets->Get(i)==NULL)) delete tgets->Get(i);
   delete tgets;
}

/****************************************************************************/
gaMoData::gaMoData  (TInputFile & inFile, TOutputFile& outFile, unsigned int nfunc):
gaMixedData(inFile, outFile) {
   iPoff 	= 0;
   //nPoff	= nfunc;
   nPoff		= nfunc+1;
   iObj		= 0;
   nObj		= nfunc;
   create();
}
//-------------------------------------------------
gaMoData::gaMoData  (TInputFile & inFile, unsigned int nfunc):
gaMixedData(inFile) {
   iPoff 	= 0;
   //nPoff	= nfunc;
   nPoff		= nfunc+1;
   iObj		= 0;
   nObj		= nfunc;
   create();
}
//-------------------------------------------------
gaMoData::gaMoData  (const char * ifname, const char * ofname, unsigned int nfunc,
	    	   unsigned int nTln, unsigned int nTcol):
gaMixedData(ifname, ofname, nTln, nTcol) {
   iPoff 	= 0;
   //nPoff	= nfunc;
   nPoff		= nfunc+1;
   iObj		= 0;
   nObj		= nfunc;
   create();
}
//-------------------------------------------------------------
TSequence * gaMoData::payoff (void)
{
   unsigned int Curr=iPoff;
   if (Curr<nPoff){
      VdecAr->Set(Curr, genome().clone());
      PoffAr->Set(Curr, new TSequence(nObj));
      for (unsigned int i=0; i<iObj; i++){

         if (!(eval==(TSequence::Evaluator)0)){
	    vCurrent(Curr).evaluator(target(i).eval);
	    pCurrent(Curr).Set(i, vCurrent(Curr).evaluate());

	}else{
	    if (!(cpnt==NULL)){
		 cpnt->target (target(i).tgID);
    		 pCurrent(Curr).Set(i, cpnt->objective(vCurrent(Curr)));
	    }else return NULL;
        	}
      }
      iPoff++;
   }else return NULL;
   return &pCurrent(Curr);
}

//-------------------------------------------------------------
TSequence::Evaluator gaMoData::Set (const TSequence::Evaluator f, char*objname)
{
   cpnt	= NULL;
   unsigned int Curr=iObj;
   if (Curr<nObj){
      tgets->Set(Curr, new TTarget);
      tgets->Get(Curr)->eval = f;
      tgets->Get(Curr)->name = objname;
      iObj++;
   }else return (TSequence::Evaluator)0;
   return eval = f;
}
//-------------------------------------------------------------
BasicES * gaMoData::Set(BasicES & thecpnt, char*objname)
{
   eval	= (TSequence::Evaluator)0;
   thecpnt.target(objname);
   unsigned int Curr=iObj;
   if (Curr<nObj){
      tgets->Set(Curr, new TTarget);
      tgets->Get(Curr)->tgID = thecpnt.target();
      tgets->Get(Curr)->name = objname;
      iObj++;
   }else return NULL;
   return cpnt = &thecpnt;
}
/****************************************************************************/
Array <TTarget*> &
gaMoData::targets (void) {
   float mini, maxi;
   for (unsigned int i=0; i<iObj; i++) {
      mini = pCurrent(0).Get(i);
      maxi = pCurrent(0).Get(i);
      for (unsigned int k=0; k<iPoff; k++) {
          if (pCurrent(k).Get(i) > maxi)
             maxi = pCurrent(k).Get(i);
          else
             if (pCurrent(k).Get(i) < mini)
          	mini = pCurrent(k).Get(i);
      }
      target(i).min = mini;
      target(i).max = maxi;
      target(i).Wnf = maxi-mini;
      target(i).Wght= 1;
   }

   for (unsigned int i=0; i<iObj; i++)
      target(i).fOpt = pCurrent(i).Get(i);

   return *tgets;
}
/****************************************************************************/
void gaMoData::evaluator (unsigned int ii){

   if (!(eval==(TSequence::Evaluator)0)){
       gaData::evaluator(target(ii).eval);
   }else{
       if (!(cpnt==NULL))
	 cpnt->target (target(ii).tgID);
   }
}
/****************************************************************************/
float gaMoData::utility (unsigned int ii){

   if (!(eval==(TSequence::Evaluator)0)){
       genome().evaluator(target(ii).eval);
       score = genome().evaluate();

   }else{
       if (!(cpnt==NULL)){
	 cpnt->target (target(ii).tgID);
    	 score	= cpnt->objective(this->genome());
       }else return score = 0;
   }
   genome().Value = score;
   return score;
}
/****************************************************************************/
float
gaMoData::minimaxi (void) {
   targets ();
   float sm = 0.0;
   for (unsigned int i=0; i<iObj; i++) {
      sm = sm + target(i).Wght * (abs(target(i).fOpt-utility(i))) / target(i).Wnf;
   }
   return sm;
}
/****************************************************************************/
void gaMoData::motitle (void)	{
   if (!ofile->isopened())
       ofile->open();

   for (unsigned int i=0; i<iObj; i++){
       ofile->FileID() <<target(i).name<<"\t";
   }
   //ofile->FileID() <<ifile->title()<<"\n";
   ofile->titleline(ifile->title());

}


/****************************************************************************/
void gaMoData::mowrite   (void) {
   if (!ofile->isopened())
       ofile->open();
   if (!(strcmp(ofilename(), "--") == 0 || ofile==NULL)){
       for (unsigned int k=0; k<iPoff; k++){
	   for (unsigned int i=0; i<iObj; i++)
	       ofile->FileID() <<pCurrent(k).Get(i)<<"\t";
	   ofile->FileID() << vCurrent(k);
       }

   }
}
/****************************************************************************/
