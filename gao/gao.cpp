
/****************************************************************************/
#include<fstream>
#include<iostream>
#include<string.h>

#include"../BasicEnergySystem/BasicES.h"
#include"gao.h"

gaParams * gao::setting = new gaParams ("gao/gasettings.txt");
/****************************************************************************/
float objective (GAGenome& g){
      GARealGenome& chrom = (GARealGenome&)g;
      gaData * usdata = (gaData*)chrom.userData();
      usdata->store(chrom);
      return usdata->evaluate();
}
//-------------------------------------------------
float mxobjective (GAGenome& g){
      MixedGenome& chrom   = (MixedGenome&)g;
      gaMixedData * usdata = (gaMixedData*)chrom.userData();
      usdata->store(chrom);
      return usdata->evaluate();
}
//-------------------------------------------------
float moobjective (GAGenome& g){
      MixedGenome& chrom   = (MixedGenome&)g;
      gaMoData * usdata = (gaMoData*)chrom.userData();
      usdata->store(chrom);
      return usdata->evaluate();
}
//-------------------------------------------------
float moutility (GAGenome& g){
      MixedGenome& chrom   = (MixedGenome&)g;
      gaMoData * usdata = (gaMoData*)chrom.userData();
      usdata->store(chrom);
      return usdata->minimaxi();
}
/****************************************************************************/



/****************************************************************************/
gao::gao(const char * ifname , const char * ofname, type GA) {
     ThisGA 	= GA;
     nSen	= 1;
     InitGen	= NULL;
     if (!(ThisGA==gao::complex)){
        data	   = new gaData   (ifname, ofname);
        InitGen	   = new TSequence (data->genome());
        alleles      = data->registered();
        obj	   = objective;
        data->ofile->open();

     }
}
//-------------------------------------------------
gao::gao(TInputFile & inFile, TOutputFile& outFile, type GA) {
     ThisGA 	= GA;
     nSen	= 1;
     InitGen	= NULL;
     if (!(ThisGA==gao::complex)){
        data	   = new gaData   (inFile, outFile);
        InitGen	   = new TSequence (data->genome());
        alleles      = data->registered();
        obj	   = objective;
        data->ofile->open();
     }
}
//-------------------------------------------------
gao::gao(TInputFile & inFile, type GA) {
     ThisGA 	= GA;
     nSen	= 1;
     InitGen	= NULL;
     if (!(ThisGA==gao::complex)){
        data	   = new gaData   (inFile);
        InitGen	   = new TSequence (data->genome());
        alleles      = data->registered();
        obj	   = objective;
        data->ofile->open();
     }
}
//-------------------------------------------------
GAGeneticAlgorithm * gao::selectga (GAGenome& g, type GA){
     switch(GA){
     case gao::simple :{
          GASimpleGA::registerDefaultParameters (*setting->params);
          pm =  gao::setting->registered ("simple");
  	 ga = (GASimpleGA*) new GASimpleGA(g);
          break;
          }
     case gao::sstate :{
          GASteadyStateGA::registerDefaultParameters (*setting->params);
          pm =  gao::setting->registered ("sstate");
  	 ga = (GASteadyStateGA*) new GASteadyStateGA(g);
          break;
          }
     case gao::stggle :{
          StruggleHGA::registerDefaultParameters (*setting->params);
          pm =  gao::setting->registered ("stggle");
  	 ga = (StruggleHGA*) new StruggleHGA(g);
          break;
          }
     case gao::inc :{
          GAIncrementalGA::registerDefaultParameters (*setting->params);
          pm =  gao::setting->registered ("inc");
  	 ga = (GAIncrementalGA*) new GAIncrementalGA(g);
          break;
          }
     case gao::deme :{
          GADemeGA::registerDefaultParameters (*setting->params);
          pm =  gao::setting->registered ("deme");
  	 ga = (GADemeGA*) new GADemeGA(g);
          break;
          }
     }
     ga->parameters(*gao::setting->params);
     return ga;
}
//-------------------------------------------------
void gao::init  (){
     chrom   = new GARealGenome (*alleles, obj, (void*)data);
     selectga (*chrom, ThisGA);

     unsigned int nMemo = (unsigned int)(ga->nGenerations()/ndiplay + 1);
     data->nGenome (nMemo);
     data->genome().copy(initgenome());
     data->genome().noID	= 0;
     data->genome().Value	= data->evaluate ();
     data->memorize ();

}

//-------------------------------------------------
void gao::start (){
     data->title();
     switch(opt){
     case gao::max :{ga->maximize(); break;}
     case gao::min :{ga->minimize(); break;}
     }
     ga->initialize();

     while(!ga->done()) {
        ga->step();
        if(ga->generation() % ndiplay == 0) {
           write	(ga->population().best());
           memorize(ga->population().best());
        }
        display(ga->population().best());
     }
     //cout << "the ga generated:\n" << ga->statistics().bestIndividual() << endl;
     std::cout << std::endl<<"the ga generated:\n" << std::endl;
     display(ga->statistics().bestIndividual());
}
//-------------------------------------------------
void gao::stop  (){
     delete chrom;
     delete ga;
}
//-------------------------------------------------
void
gao::maximize (GAGenome::Evaluator fg, type GA)
{
      search (gao::max);
      ThisGA 	= GA;
      obj	= fg;
      init ();
      start();
      stop ();
}
//-------------------------------------------------
void
gao::minimize (GAGenome::Evaluator fg, type GA)
{
      search (gao::min);
      ThisGA 	= GA;
      obj	= fg;
      init ();
      start();
      stop ();
}
//-------------------------------------------------
void
gao::maximize (TSequence::Evaluator fs, type GA)
{
      search (gao::max);
      ThisGA 	= GA;
      data->evaluator(fs);
      init ();
      start();
      stop ();
}
//-------------------------------------------------
void
gao::minimize (TSequence::Evaluator fs, type GA)
{
      search (gao::min);
      ThisGA 	= GA;
      data->evaluator(fs);
      init ();
      start();
      stop ();
}

//-------------------------------------------------
void
gao::maximize (BasicES & thecpnt, type GA)
{
      search (gao::max);
      ThisGA 	= GA;
      data->evaluator(thecpnt);
      init ();
      start();
      stop ();
}
//-------------------------------------------------
void
gao::minimize (BasicES & thecpnt, type GA)
{
      search (gao::min);
      ThisGA 	= GA;
      data->evaluator(thecpnt);
      init ();
      start();
      stop ();
}
//-------------------------------------------------
void gao::write (const GAGenome& g) {
      GARealGenome& chrom = (GARealGenome&)g;
      gaData * usdata = (gaData*)chrom.userData();
      usdata->store(chrom);
      usdata->genome().noID = ga->generation();
      usdata->genome().Value=chrom.score();
      usdata->write ();
}
//-------------------------------------------------
void gao::display (const GAGenome& g) {
      GARealGenome& chrom = (GARealGenome&)g;
      gaData * usdata = (gaData*)chrom.userData();
      usdata->store(chrom);
      usdata->genome().noID = ga->generation();
      usdata->genome().Value=chrom.score();
      usdata->display ();
}
//-------------------------------------------------
void gao::write (const GAPopulation& p){
      for(int i=0; i<p.size(); i++)
         gao::write(p.individual(i));
}
//-------------------------------------------------
void gao::display (const GAPopulation& p){
      for(int i=0; i<p.size(); i++)
         gao::display(p.individual(i));
}
//-------------------------------------------------
TSequence * gao::memorize (const GAGenome& g) {
      GARealGenome& chrom = (GARealGenome&)g;
      gaData * usdata = (gaData*)chrom.userData();
      usdata->store(chrom);
      usdata->genome().noID = ga->generation();
      usdata->genome().Value=chrom.score();
      return usdata->memorize ();
}
//-------------------------------------------------
bool gao::historic (const char * ofname){
     if (!(data->cpnt==NULL)){
         resultfile (ofname);
         data->cpnt->titleline(data->ofileID());
         return data->historic (data->ofileID());
     }else return false;
     return true;
}
//-------------------------------------------------
bool gao::sensitivity (unsigned int noGene, const char * ofname){
     resultfile (ofname);
     if (!(data->eval==(TSequence::Evaluator)0))
	data->title();
     else
	data->cpnt->titleline(data->ofileID());
     return data->sensitivity (data->ofileID(), noGene);

}
//-------------------------------------------------
bool gao::sensitivities (unsigned int nbGene, const char * ofname){
     resultfile (ofname);
     nSen = nbGene;
     bool bb;
     if (!(data->eval==(TSequence::Evaluator)0))
	data->title();
     else
	data->cpnt->titleline(data->ofileID());

     for (unsigned int ii=0; ii<nSen; ii++)
     	bb = data->sensitivity (data->ofileID(), ii);
     return bb;

}

/****************************************************************************/




/****************************************************************************/
mxgao::mxgao(const char * ifname, const char * ofname, type GA ):
gao(ifname, ofname, gao::complex)
{
     ThisGA 	= GA;
     nSen	= 1;
     InitGen	= NULL;
     if (!(ThisGA==gao::mobjective)){
     	mxdata	= new gaMixedData   (ifname, ofname);
     	InitGen	= new TSequence (mxdata->genome());
     	mxdata->registered ();
     	mxalleles= mxdata->mixedalleles();
     	obj	   	= mxobjective;
     	mxdata->ofile->open();
     	data	= mxdata;
     }
}
mxgao::mxgao(TInputFile & inFile, TOutputFile& outFile, type GA):
gao(inFile, outFile, gao::complex)
{
     ThisGA 	= GA;
     nSen	= 1;
     InitGen	= NULL;
     if (!(ThisGA==gao::mobjective)){
     	mxdata	= new gaMixedData   (inFile, outFile);
     	InitGen	= new TSequence (mxdata->genome());
     	mxdata->registered ();
     	mxalleles    = mxdata->mixedalleles();
     	obj	   	= mxobjective;
     	mxdata->ofile->open();
     	data	= mxdata;
     }
}

mxgao::mxgao(TInputFile & inFile, type GA):
gao(inFile, gao::complex)
{
     ThisGA 	= GA;
     nSen	= 1;
     InitGen	= NULL;
     if (!(ThisGA==gao::mobjective)){
     	mxdata	= new gaMixedData   (inFile);
     	InitGen	= new TSequence (mxdata->genome());
     	mxdata->registered ();
     	mxalleles    = mxdata->mixedalleles();
     	obj	   	= mxobjective;
     	mxdata->ofile->open();
     	data	= mxdata;
     }
}
//-------------------------------------------------
void mxgao::init  (){
     chrom   = new MixedGenome (*mxalleles, obj, (void*)mxdata);
     gao::selectga (*chrom, ThisGA);

     unsigned int nMemo = (unsigned int)(ga->nGenerations()/ndiplay + 1);
     mxdata->nGenome (nMemo);
     mxdata->genome().copy(initgenome());
     mxdata->genome().noID	= 0;
     mxdata->genome().Value= mxdata->evaluate ();
     mxdata->memorize ();
}
//-------------------------------------------------
void mxgao::start (){
     mxdata->title();
     switch(opt){
     case gao::max :{ga->maximize(); break;}
     case gao::min :{ga->minimize(); break;}
     }
     ga->initialize();
     while(!ga->done()) {
        ga->step();
        if(ga->generation() % ndiplay == 0) {
           write	(ga->population().best());
           memorize(ga->population().best());
        }
        display(ga->population().best());
     }
     //cout << "the ga generated:\n" << ga->statistics().bestIndividual() << endl;
     std::cout << std::endl<<"the ga generated:\n" << std::endl;
     display(ga->statistics().bestIndividual());
}
//-------------------------------------------------

void mxgao::write (const GAGenome& g) {
      MixedGenome& chrom = (MixedGenome&)g;
      gaData * usdata = (gaData*)chrom.userData();
      usdata->store(chrom);
      usdata->genome().noID = ga->generation();
      usdata->genome().Value=chrom.score();
      usdata->write ();
}
//-------------------------------------------------
void mxgao::display (const GAGenome& g) {
      MixedGenome& chrom = (MixedGenome&)g;
      gaData * usdata = (gaData*)chrom.userData();
      usdata->store(chrom);
      usdata->genome().noID = ga->generation();
      usdata->genome().Value= chrom.score();
      usdata->display ();
}
//-------------------------------------------------
void mxgao::write (const GAPopulation& p){
      for(int i=0; i<p.size(); i++)
         mxgao::write(p.individual(i));
}
//-------------------------------------------------
void mxgao::display (const GAPopulation& p){
      for(int i=0; i<p.size(); i++)
         mxgao::display(p.individual(i));
}
//-------------------------------------------------
TSequence *  mxgao::memorize (const GAGenome& g) {
      MixedGenome& chrom = (MixedGenome&)g;
      gaData * usdata = (gaData*)chrom.userData();
      usdata->store(chrom);
      usdata->genome().noID = ga->generation();
      usdata->genome().Value=chrom.score();
      return usdata->memorize ();
}


/****************************************************************************/

/****************************************************************************/
mogao::mogao(const char * ifname, const char * ofname, unsigned int nfunc ) :
mxgao(ifname, ofname, gao::mobjective) {

     modata	= new gaMoData (ifname, ofname, nfunc);
     InitGen	= new TSequence (modata->genome());
     modata->registered ();
     mxalleles    = modata->mixedalleles();
     obj	   	= moobjective;
     modata->ofile->open();
     mxdata	= modata;
     data	= modata;

     nSen	= 1;
     //modata = new gaMoData   (ifname, ofname, nfunc);
     nMdl   = nfunc;
     iMdl   = 0;
     gamdl  = new Array<gaotype*> (nMdl);
     for (unsigned int i=0; i<nMdl; i++)
        gamdl->Set(i, NULL);

}
//-------------------------------------------------
mogao::mogao(TInputFile & inFile, TOutputFile& outFile, unsigned int nfunc) :
mxgao(inFile, outFile, gao::mobjective) {

     modata = new gaMoData   (inFile, outFile, nfunc);
     InitGen	= new TSequence (modata->genome());
     modata->registered ();
     mxalleles    = modata->mixedalleles();
     obj	   	= moobjective;
     modata->ofile->open();
     mxdata	= modata;
     data	= modata;

     nSen	= 1;
     //modata = new gaMoData   (inFile, outFile, nfunc);
     nMdl   = nfunc;
     iMdl   = 0;
     gamdl  = new Array<gaotype*> (nMdl);
     for (unsigned int i=0; i<nMdl; i++)
        gamdl->Set(i, NULL);
}
//-------------------------------------------------
mogao::mogao(TInputFile & inFile, unsigned int nfunc) :
mxgao(inFile, gao::mobjective) {

     modata = new gaMoData   (inFile, nfunc);
     InitGen	= new TSequence (modata->genome());
     modata->registered ();
     mxalleles    = modata->mixedalleles();
     obj	   	= moobjective;
     modata->ofile->open();
     mxdata	= modata;
     data	= modata;

     nSen	= 1;
     //modata = new gaMoData   (inFile, nfunc);
     nMdl   = nfunc;
     iMdl   = 0;
     gamdl  = new Array<gaotype*> (nMdl);
     for (unsigned int i=0; i<nMdl; i++)
        gamdl->Set(i, NULL);
}

/****************************************************************************/
void mogao::Set (TSequence::Evaluator fs, char*objname, type GA, extramum EX)
{
     modata->Set(fs, objname);
     unsigned int Curr=iMdl;
     if (Curr<nMdl){
	modele(Curr, new gaotype);
	modele(Curr).Name = objname;
	modele(Curr).GA   = GA;
	modele(Curr).EX   = EX;
	filename (Curr);
	iMdl++;
     };
}
/****************************************************************************/
void mogao::Set (BasicES & thecpnt, char*objname, type GA, extramum EX)
{
     //thecpnt.target(objname);
     modata->Set(thecpnt, objname);
     unsigned int Curr=iMdl;
     if (Curr<nMdl){
	modele(Curr, new gaotype);
	modele(Curr).Name = objname;
	modele(Curr).GA   = GA;
	modele(Curr).EX   = EX;
	filename (Curr);
	iMdl++;
     };
}
/****************************************************************************/
void mogao::filename (unsigned int i)
{
     dynstring corps= modele(i).Name;
     dynstring ext1 = ".rslt";
     dynstring ext2 = ".htry";
     dynstring ext3 = ".svty";

     dynstring fname   = corps + ext1;
     modele(i).rfname  = new dynstring (fname.liste());
     fname.empty();

     fname  = corps + ext2;
     modele(i).hfname  = new dynstring (fname.liste());
     fname.empty();

     fname  = corps + ext3;
     modele(i).sfname  = new dynstring (fname.liste());
     fname.empty();
}
/****************************************************************************/
void mogao::init  (){
     chrom   = new MixedGenome (*mxalleles, obj, (void*)modata);
     gao::selectga (*chrom, ThisGA);

     unsigned int nMemo = (unsigned int)(ga->nGenerations()/ndiplay + 1);
     modata->nGenome (nMemo);
     mxdata->genome().copy(initgenome());
     modata->genome().noID	= 0;
     modata->genome().Value= modata->evaluate ();
     modata->memorize ();
}
//-------------------------------------------------
void mogao::start (){
     modata->title();
     switch(opt){
     case gao::max :{ga->maximize(); break;}
     case gao::min :{ga->minimize(); break;}
     }
     ga->initialize();
     while(!ga->done()) {
        ga->step();
        if(ga->generation() % ndiplay == 0) {
           write	(ga->population().best());
           memorize(ga->population().best());
        }
        display(ga->population().best());
     }
     std::cout << std::endl<<"the ga generated:\n" << std::endl;
     display(ga->statistics().bestIndividual());
}
/****************************************************************************/
void mogao::optimize (unsigned int i)
{
     resultfile (modele(i).rfname->liste());		//Set the file name of the ieme target
     modata->evaluator (i);			//Select the ieme objectif
     ThisGA = modele(i).GA;			//Select the type of the GA
     search  (modele(i).EX); 			//Define the direction of the research (min or max)
     init  (); 					//Initialization
     start (); 					//Optimization and write whole data
     stop  ();					//End of the optimization
     modata->payoff ();				//evaluate all objective with the same genome
     modata->evaluator (i); 			//Select again the current objectif
     historic (modele(i).hfname->liste());		//Give the history of the evaluation (VI)
     sensitivities (nSen, modele(i).sfname->liste());	//First variable sensitivity
}
/****************************************************************************/
void mogao::optimize (void)
{
//     char * Oldfile = resultfile();
     for (unsigned int i=0; i<iMdl; i++)
     	optimize (i);

     resultfile ("utility.rslt");
     search (gao::min);
     ThisGA = gao::sstate;
     obj	= moutility;
     init ();
     start();
     stop ();
     modata->payoff ();
     modata->evaluator (0); 			//Select the first objectif
     historic ("utility.htry");			//Give the history of the evaluation (VI)
     //sensitivity (0, "utility.svty");		//First variable sensitivity
     sensitivities (nSen, "utility.svty");		//n-First variables sensitivity

     resultfile ("payofftable.txt");
     modata->motitle ();
     modata->mowrite ();
}
/****************************************************************************/
void mogao::write (const GAGenome& g) {
      mxgao::write (g);
}
/****************************************************************************/
void mogao::display (const GAGenome& g) {
      mxgao::display (g);
}
/****************************************************************************/
void mogao::write (const GAPopulation& p){
      mxgao::write (p);
}
/****************************************************************************/
void mogao::display (const GAPopulation& p){
      mxgao::display (p);
}
/****************************************************************************/
TSequence *  mogao::memorize (const GAGenome& g) {
      return mxgao::memorize (g);
}

/****************************************************************************/
