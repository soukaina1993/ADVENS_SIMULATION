//
// Created by robin.lemaire on 12.07.2021.
//
#include <fstream>
#include "iofile/iofile.h"
#include "BasicEnergySystem/exchanger/TExchanger.h"
#include "BasicEnergySystem/BasicES.h"


using namespace std;

const char* file_inEXCH = "BasicEnergySystem/exchanger/exchparam.txt";


int main(){


    //Take in input the name corresponding to the one in the exchparamfile

    TExchanger ex(TExchanger::P1);
    cout<<"Creation"<<endl;

    //exchParams exchParams(inputFileEXCH);
    //exchParams.RegisterParams(&exchangerTest, "P1");



    cout<<"Init"<<endl;
ex.Hot.SpecHeat=4180;
ex.Cold.SpecHeat=4180;
 //   cout<<ex.flw<<endl;
        ex.performance();

    ex.optimize(TExchanger::COST);
    ex.display();


    return 0;
}