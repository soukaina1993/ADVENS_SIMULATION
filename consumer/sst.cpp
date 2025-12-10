//
// Created by lucile.schulthe on 19.01.2022.
//

#include "sst.h"

void sst::calcul_flow(Tnetwork &network,double & power){

    if(dT>0 && network.flow->StateIn->T<connec.Hot.StateIn->T)
            cout<<"Attention le reseau n'est pas assez chaud pour combler les besoins"<<endl;
    if(dT<0 && network.flow->StateIn->T>connec.Cold.StateIn->T)
            cout<<"Attention le reseau n'est pas assez froid pour combler les besoins"<<endl;

    flow->StateIn->T = network.flow->StateIn->T;
    flow->StateIn->P=network.flow->StateIn->P;
    flow->StateIn->P=network.flow->StateIn->P-( network.alt_plant-alt)*9.81;
    flow->StateOut->P=flow->StateIn->P-dP;
    flow->StateCalculation(flow->StateIn,network.fluid);
    flow->StateCalculation(flow->StateOut,network.fluid);
    flow->MassflowCalculation(power);
    connec.dP=connec.getloss(); // encore à voir
    flow->setStateOut_Ph(flow->StateOut->P-connec.dP,flow->StateIn->h-flow->Power/ *flow->Massflow);

}


void sst::init_connec(){

    if(sst_type=="HP"){
     //   THeatpump HP("R134a",0.8,Tpinch/2,Tpinch/2);
     //   HP.init();
     //   HP.typ=THeatpump::WATER;
        // connec=HP;
        cout<<"Heat pump"<<endl;
    }
    else if(sst_type=="HX"){
        cout<<"Heat echanger"<<endl;
    }
    else{
        cout<<"Prends les valeurs de reference"<<endl;
        // eventuellement recalibre le pinch et tout
    }

    if(dT>0){
        connec.sor=HEAT;
        connec.Cold.StateOut->T=Tmin;
        connec.Hot.StateIn->T=Tmin+Tpinch;
        connec.Hot.StateOut->T=connec.Hot.StateIn->T-dT;
        flow=&connec.Cold;
        cout<<"cold"<<endl;
    }
    else if(dT<0){
        connec.sor=COOL;
        connec.Hot.StateOut->T=Tmin;
        connec.Cold.StateIn->T=Tmin-abs(Tpinch);
        connec.Cold.StateOut->T=connec.Cold.StateIn->T-dT;
        flow=&connec.Hot;
    }

};

void sst::Display(void){
    cout << endl;
    cout <<"/****************************************************************************/" << endl;
    cout<<"                         Building                                "<<endl<<endl;
    cout<<"Energy total                                      [kwh] :       :"<<energy<<endl;
    cout<<"Power max total                                   [kwh] :       :"<<power_max<<endl;
    cout<<"Minimal temperature distribution network          [degC] :      :"<<Tmin-Tpinch<<endl;
    cout<<"Difference temperature                            [degC] :      :"<<dT<<endl;
    cout<<endl;
    cout <<"/****************************************************************************/" << endl;
    cout << endl;
}





void sst::input_file(vector <string> &data_input){
    const char* file_in=const_cast<char*>(data_input[1].c_str());
    string type=data_input[1];
    int Time;
    sst_type=data_input[4];
    dT=stod(data_input[5]); // est-ce que cest cote cold ou hot?
    Tmin=stod(data_input[6]);
    Tpinch=stod(data_input[7]);
    L_racc=stod(data_input[10]);
    DN_racc=stod(data_input[11]);
    alt=stod(data_input[12]);
    if(type=="Power")
        Getdata(file_in,power,Time);
    else if(type=="Flow"){
        double a =1/dT/4186;
        Getdata(file_in,flow_vec,Time);
        power= multi(flow_vec,a);}
    else if(type=="Power-Flow"){
        Getdata(file_in,power,Time);
        Getdata(file_in,flow_vec,Time,2);}
    else  cout<<"Unknown parameter"<<endl;


};

