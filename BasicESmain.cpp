#include <iostream>
#include <string>
#include <fstream>
#include "iofile/iofile.h"
#include "BasicEnergySystem/exchanger/TTubular.h"
#include "BasicEnergySystem/exchanger/TPlate.h"
#include "fluids/PengRobinson.h"
#include "fluids/MWater.h"
#include "gao/ga/ga.h"
#include "gao/gao.h"



using namespace std;



float fObjective(GAGenome &);
float fObjective2(TSequence& vect);
float mxObjective1(TSequence& vect);




int main (){

    // optimizer

    cout << "mk ga main test"<<endl;
        mxgao mxga ("mxdata.txt", "mxxrslts.txt"); //11-05-2020 le fichier joue le rôle de génome (avec des allèles)
    mxga.minimize(fObjective2, gao::simple);
    cout<<"That's it!"<<endl;


    ////plate1
    const char* file_inFLUID = "fluids/LKDataFluidsNEW5.txt";
    //DECLARATION DES FLUIDES
    FluidCore *fluidMain=new MWater;
    FluidCore *fluidPlateHot=new MWater;

    //DECLARATION DE L ECHANGEUR
    TPlate plateOne("P1",TExchanger::COUNTER, TExchanger::DEFAULT);
    plateOne.init();

    //On initialise les fluides chauds et froids (récuperer les données des états)
    plateOne.Hot.init(fluidPlateHot);
    plateOne.Cold.init(fluidMain);

    plateOne.performance();
    cout <<"/****************************************************************************/"<<endl;
    cout<<"---------------------Plates Heat Exchanger n°1-------------------------------"<<endl;
    plateOne.display();


        /////Evaporator
   FluidCore *fluidEvapHot=new  MWater("H2O");
    TTubular evap("Cond",TExchanger::COUNTER, TExchanger::EVAP, TTubular::HOT);


    evap.Cold.pStateIn.SetT(plateOne.Cold.pStateOut.T+273.15);
    evap.Cold.pStateIn.SetP(1E5);
    evap.Cold.Massflow=plateOne.Cold.Massflow;

    evap.Hot.pStateIn.SetT(423.15);
    evap.Hot.pStateOut.SetT(323.15);
    evap.Hot.pStateIn.SetP(1E5);



    double Ta;
    Ta=((evap.Cold.pStateIn.GetT()+evap.Cold.pStateOut.GetT())/2+(evap.Hot.pStateIn.GetT()+evap.Hot.pStateOut.GetT())/2)/2;
    evap.Cold.pStateWall->SetT(Ta);
    evap.Hot.pStateWall->SetT(Ta);

    evap.Hot.init(fluidEvapHot);
    evap.Cold.init(fluidMain);

    evap.performance();


    cout <<"/****************************************************************************/"<<endl;
    cout<<"-----------------------------Evaporator--------------------------------------"<<endl;

    evap.display();

    evap.Q=1;

    evap.performance();
    cout<<evap.Q<<endl;
    evap.display();
     return 0;


}


//take in input the chromose. To quantify how good it is
float
fObjective2(TSequence& vect)
{
    //fixe
    float HotTin =  353.0;
    float ColdTin = 300.0;

    //inconnue
    float HotTout =  vect.Get(0);
    float ColdTout  = vect.Get(1);

    //fichier param
    float k = 4066.84;
    //fixer mais est cencé être calculé
    float HotSpecHeat = 4000; //cp


    float Q= 10000.0;
    float DeltaTlog = 0.0;
    float LeftDeltaT = abs(HotTin - ColdTin); //co-courant
    float RightDeltaT = abs(HotTout - ColdTout); //co-courant
    DeltaTlog=(LeftDeltaT-RightDeltaT)/ln(LeftDeltaT/RightDeltaT);


    float Pinch = 0.0;

    if (LeftDeltaT<RightDeltaT)
    {
        Pinch	=LeftDeltaT;
    }

    else
    {
        Pinch	=RightDeltaT;
    }
  //  cout << "Pinch" << Pinch << endl;
  //return Q / k / DeltaTlog + abs(Pinch);

    if(HotTout < ColdTout)
  {
      if (Pinch >=0 && Pinch <= 100){
         // cout << "Pinch" << Pinch << endl;
          cout << "area" << Q / k / DeltaTlog << endl;
          return Q / k / DeltaTlog ;

      }
      else
      {
          return 1000;
      }
  }
  else
  {
      return 1000;
  }









        //flow HotSpecHeat=(Hot.tStateIn.GetH()-Hot.tStateOut.GetH())/(Hot.tStateIn.GetT()-Hot.tStateOut.GetT());

    //Q = HotMassflow * HotSpecHeat * abs(HotTin - HotTout);
    Q= 10000.0;

    //cout << "Q value" << Q  << endl;
    //cout << "surface" <<  Q/k/DeltaTlog << endl;


    //eturn (HotMassflow * HotSpecHeat * abs(HotTin - HotTout))/k/DeltaTlog;

}
