//
// Created by cornelia.blanke on 08.03.2023.
//

#include "network/QNetwork.h"
#include "fluids/PengRobinson.h"
#include "fluids/ExtCoolProp.h"
#include "BasicEnergySystem/heatpump/THeatpump.h"
#include "network/TElectronHP.h"


using namespace std;

int main(){
    cout << "Create QNetwork..." << endl;
//    QNetwork mynetwork;
//    // test access of attributes
//    mynetwork.Z = 2;
//    mynetwork.N = 1;
//    // test access of heritated attributes
//    mynetwork.m = mynetwork.Z;
//    // test display function
//    mynetwork.display();

    QNetwork mynetwork;
/*** Test Tour-de-Peilz ***/
//    mynetwork.read_networktable((char*)"C:\\Users\\cornelia.blanke\\Documents\\adv\\TourDe2.csv");
//    double Mdot_of_all_electrons[] = {38.38,2.82,2.95,1.75,51.18,
//            7.14,10.92,4,1.41,4.38,
//            1.47,2.48,10.88,6.18,8.49,
//            8.18,6.62,4.99,3.22,5.54};
//    double mu_of_all_electrons[] = {4.2,0.31,0.32,0.19,5.59,
//                                    0.78,1.19,0.44,0.15,0.48,
//                                    0.16,0.27,1.19,0.68,0.93,
//                                    0.89,0.72,0.55,0.35,0.61};

/*** Test Estavayer ***/
    mynetwork.read_networktable((char*)"C:\\Users\\cornelia.blanke\\Documents\\adv\\Esta_NEW.csv");
    double Mdot_of_all_electrons[] = {0.242539943, 0.38877595, 0.455692426, 0.914188189, 1.074035312,
                                         0.195264975, 0.093932055, 0.126004202, 2.381312899, 0.794739696,
                                         0.752609309, 0.42928054, 1.930920464, 0.341876066, 0.264407073,
                                         1.135406852, 0.350180986, 1.251091285, 0.421840272, 0.529806894,
                                         0.941381134, 0.095498568, 0.110387111, 0.089253866, 0.061256787,
                                         0.35885952, 0.257132262, 0.23571347, 0.269135967, 0.256353584,
                                         0.231743547, 0.220914674, 0.284830704, 0.290762183, 0.280819524,
                                         0.005337352, 0.124357629, 0.124357629, 0.124357629, 0.098735672,
                                         0.071653949, 0, 0, 0.079670651, 0.007885937,
                                         0.043141815, 2.67E-06, 2.67E-06, 0, 0, 0, 0, 0};
    mynetwork.compute_ANZ();
    mynetwork.display();
    mynetwork.compute_m();
    cout << "m: " << mynetwork.m <<endl;
    mynetwork.compute_mtable();
    mynetwork.display_mtable();

    cout << endl << "Create TBranch..." << endl;
    TBranch mybranch;
    mybranch.QNet = &mynetwork;
    mybranch.l = -2;
    mybranch.compute_m();
    cout << "m: " << mybranch.m << endl;
    mybranch.compute_m_from_table();
    cout << "m: " << mybranch.m << endl;
    mybranch.compute_sstnumber();
    cout << "m: " << mybranch.sstnumber << endl;
    mybranch.compute_sstnumber_from_table();
    cout << "m: " << mybranch.sstnumber << endl;

    cout << endl << "Create TSegment..." << endl;
    TSegment mysegment;
    mysegment.QNet = &mynetwork;
    mysegment.l = -2;
    mysegment.k = 1;
    mysegment.compute_m();
    cout << "m: " << mysegment.m <<endl;
    mysegment.compute_m_from_table();
    cout << "m: " << mysegment.m <<endl;

    cout << endl << "Create TProton..." << endl;
    TProton myproton;
    myproton.QNet = &mynetwork;
    myproton.l = -11;
    myproton.k = 1;
    myproton.compute_m();
    cout << "m: " << myproton.m <<endl;
    myproton.compute_m_from_table();
    cout << "m: " << myproton.m <<endl;

    cout << endl << "Create TElectrons..." << endl;
    TElectron myelectron;
    TElectron myelectron2(TElectron::NEUTRAL);
    TElectron myelectron3(TElectron::POS);
    myelectron.QNet = &mynetwork;
    myelectron.ID = 10;
    myelectron.find_kl_from_ID();
    cout << "ID, k, l: " << myelectron.ID << ", " << myelectron.k << ", " << myelectron.l << endl;


    cout << endl << "Create TNeutron..." << endl;
    TNeutron myneutron;
    myneutron.QNet = &mynetwork;
    myneutron.n = 3;
    myneutron.compute_m();
    cout << "m: " << myneutron.m <<endl;
    myneutron.m = 0;
    myneutron.compute_m_from_table();
    cout << "m: " << myneutron.m <<endl;

    cout << endl << "Create TCentralPlant..." << endl;
    TCentralPlant myplant;

    cout << endl << "Test muTable..." << endl;
    //vector<TElectron> elec(mynetwork.Z);
    //mynetwork.define_electrons(elec);
    mynetwork.create_electrons();
    //mynetwork.create_protons_from_electrons();
    mynetwork.create_neutrons();
    //mynetwork.define_centralplant(myplant);
    //myplant.max_massflow = 182.98;
    mynetwork.create_centralplant();
    mynetwork.CentralPlant->max_massflow = 182.98;          // TODO
    //mynetwork.init_electrons(mu_of_all_electrons, "mu");
    mynetwork.init_electrons(Mdot_of_all_electrons, "Mdot");
    mynetwork.compute_muTable();
    mynetwork.display_muTable();

    mynetwork.compute_mu_from_table();
    mynetwork.compute_massflow_from_mu();
    cout << "mu Network: " << mynetwork.mu <<
        ", Mdot Network: " << mynetwork.flow->massflow() << endl;
    myneutron.compute_mu_from_table();
    myneutron.compute_massflow_from_mu();
    cout << "mu Neutron " << myneutron.n << ": " << myneutron.mu <<
        ", Mdot Neutron: " << myneutron.flow->massflow() << endl;
    mysegment.compute_mu_from_table();
    mysegment.compute_mubar_from_mu();
    mysegment.compute_massflow_from_mu();
    cout << "mu Segment: " << mysegment.mu << "\t mu_bar: " << mysegment.mu_bar <<
        ", Mdot Segment: " << mysegment.flow->massflow() <<endl<<endl;



/*** How to use fluid library ***/
    mynetwork.flow->PIn (100000);
    mynetwork.flow->TIn (353.15);
    mynetwork.flow->POut (100000);
    mynetwork.flow->TOut (333.15);

    mynetwork.init_fluid("Water");
    mynetwork.flow->StateIn->Display();
    cout << endl;
    mynetwork.flow->pStateIn->Display();

    cout << "Cp: " << mynetwork.flow->Cp << endl;
    mynetwork.compute_k0();
    cout << "k0: " << mynetwork.k0 << endl;
    mynetwork.Tambient = & mynetwork.Toutside;
    mynetwork.compute_nubar();
    cout << "nubar: " << mynetwork.nu_bar << endl << endl;

    cout << "deltabar: " << mynetwork.delta_bar << endl << endl;

/*** Couple segment with HPipe ***/
    mysegment.flow->PIn (500000);
    mysegment.flow->TIn(350);
    mysegment.flow->POut (500000);
    mysegment.flow->TOut(350);
    mysegment.init_fluid("Water");

    mysegment.Pipe.Diam = 0.05;
    mysegment.Pipe.Length = 1.0;
    mysegment.Pipe.RugAbs = 0.0001;
    mysegment.Pipe.LossCalculation(293);
    cout << "Pressure Loss: " << mysegment.Pipe.DeltaPr << endl;
    cout << "Temperature Loss: " << mysegment.Pipe.DeltaT << endl;
    mysegment.Pipe.U_pipe = 100;
    mysegment.Pipe.LossCalculation(293);
    cout << "Temperature Loss: " << mysegment.Pipe.DeltaT << endl;
    //mysegment.Pipe.Display();


/*** add information about all segments to mynetwork ***/
    mynetwork.read_diamtable((char*)"C:\\Users\\cornelia.blanke\\Documents\\adv\\Esta_NEW_diam.csv");
    mynetwork.read_lengthtable((char*)"C:\\Users\\cornelia.blanke\\Documents\\adv\\Esta_NEW_length.csv");
    //mynetwork.read_pipetable((char*)"C:\\Users\\cornelia.blanke\\Documents\\adv\\Esta_NEW_pipes.csv");
    mynetwork.read_pipetable((char*)"C:\\Users\\cornelia.blanke\\Documents\\adv\\Esta_NEW_pipes2.csv");

    mynetwork.create_segments();

    for (int n=0; n<mynetwork.nbBranches; n++) {      // loop over table
        for (int m = 0; m < mynetwork.lengthtable.size(); m++) {
            if (mynetwork.Table_of_Segments[m][n] != nullptr) {
                mynetwork.Table_of_Segments[m][n]->init_DN();
                mynetwork.Table_of_Segments[m][n]->init_length();
                mynetwork.Table_of_Segments[m][n]->init_pipeprop();
                mynetwork.Table_of_Segments[m][n]->flow->PIn(500000);
                mynetwork.Table_of_Segments[m][n]->flow->POut(500000);
                mynetwork.Table_of_Segments[m][n]->flow->TIn(350);
                mynetwork.Table_of_Segments[m][n]->flow->TOut(350);
                mynetwork.Table_of_Segments[m][n]->init_fluid("Water");
                mynetwork.Table_of_Segments[m][n]->compute_mu_from_table();
                mynetwork.Table_of_Segments[m][n]->compute_massflow_from_mu();
                mynetwork.Table_of_Segments[m][n]->compute_dP();
                mynetwork.Table_of_Segments[m][n]->compute_tau();

                cout << "Segment " <<
                     ": l=" << mynetwork.Table_of_Segments[m][n]->l << ", k=" << mynetwork.Table_of_Segments[m][n]->k <<
                     ", DN=" << mynetwork.Table_of_Segments[m][n]->DN << ", L="
                     << mynetwork.Table_of_Segments[m][n]->length <<
                     ", dP=" << mynetwork.Table_of_Segments[m][n]->dP << ", tau="
                      << mynetwork.Table_of_Segments[m][n]->tau <<
                      endl;
            }
        }
    }


/*** Buildings, Districts, Environment ***/
    Building B;
    B.define_heatmode("Cool");
    char* file = (char*) "C:\\Users\\cornelia.blanke\\OneDrive - HESSO\\Lucile\\ADVENS\\adv - Copie_old\\dataBase\\Fribourg_B.txt";
    B.init_norm(file, "9032961");
    //B.Display();

    char* file_meteo = (char*) "C:\\Users\\cornelia.blanke\\OneDrive - HESSO\\Lucile\\ADVENS\\adv - Copie_old\\dataBase\\Meteo_Fribourg_6mois.txt";
    environment fribourg(file_meteo,fribourg.T,"2018");   // starts in August
    //fribourg.Display();

    B.calcul_power(fribourg);

    District D(file);
    //D.Display();
    D.calcul_power(fribourg);

    vector<double> energy = {20000.0, 15000.0, 22000.0};
    vector<double> Tref = {330.0, 350.0, 345.0};
    vector<int> type = {1, 4, 6};
    vector<const char*> heat_type = {"Cool", "Heat", "ECS"};
    District D2(energy, Tref, type, heat_type);
    D2.calcul_power(fribourg);
    District D3(12000, 355.0, 1, "HeatECS");
    D3.calcul_power(fribourg);

    myelectron.D = D;
    //myelectron.compute_massflow_from_power();

    cout << "start " << endl;
    vector<vector<string>> content = mynetwork.read_csv(file_meteo);
    vector<double> temperature(content.size()-1);
    for (int i=0; i<content[0].size(); i++)
        if (content[0][i] == "2019"){
            for (int j=1; j<content.size(); j++)
                temperature[j-1] = stod(content[j][i]);
        }
    cout << temperature[0] << " " << temperature[content.size()-2] << endl;

    cout << "end" << endl;

    mynetwork.read_settings((char*)"C:\\Users\\cornelia.blanke\\Documents\\adv\\Cases\\Esta_NEW\\Esta_NEW_settings.txt");
//    for (const auto& setting : mynetwork.settings)
//        cout << setting.first << " " << setting.second << endl;

    mynetwork.CentralPlant->init_boilers();
    for (auto& boiler : mynetwork.CentralPlant->boilers)
        cout << boiler.type << " " << boiler.max_power << " " << boiler.res_cost_per_kWh << endl;
    cout << endl;
    unordered_map<Boiler::Type, double> new_prices {{Boiler::WOOD, 0.13}};
    mynetwork.CentralPlant->modify_pricelist(new_prices);
    mynetwork.CentralPlant->init_boilers();
    for (auto& boiler : mynetwork.CentralPlant->boilers)
        cout << boiler.type << " " << boiler.max_power << " " << boiler.res_cost_per_kWh << endl;

    /*** fluid libraries ***/
    PRfluid prfluid((char*)"R134a");
    double h=271230, P=418900;
    double v=prfluid.v_Ph_a(P, h);
    cout << v << endl;
    double T=prfluid.T_vP( v, P);
    double s,p,u,htemp;
    prfluid.Prop_vT(v,T,p,u,htemp,s);
    cout << s << endl;

    double Ts = prfluid.T_P_s(P);
    double vl,vg;
    prfluid.v_PT(P,Ts,vl,vg);
    double x = Mfraction(vl,vg,v);
    cout << x << endl;
    prfluid.Prop_Tx_s(Ts,x,v,p,u,htemp,s);
    cout << s << endl;
//
//    MWater water((char*)"H2O");
//    double h, u, s= 648.5, dpv;
//    P=4000000;
//    v = water.v_Ps_a(P, s);
//    cout << v << endl;
//
//    T = water.T_vP(v, P);
//    cout << T << endl;
//    water.Prop_vT(v,T,P,u,h,s);
//    cout << h << endl;

    /*** Heat Pumps ***/
    THeatpump myHP("R134a");
    double Thot=60+273.15, Tcold=10+273.15, Text=20+273.15;
    cout << "COP_max = " << myHP.COP_max(Thot, Tcold) << endl;
    myHP.COP_real(Thot, Tcold);     // from cycle
    cout << myHP.COP << endl;
    myHP.eta_comp = 0.91;
    myHP.COP_real(Thot, Tcold, Text);   // from eta_ex
    cout << myHP.COP << endl;



    /*** Integrate CoolProp ***/
    ExtCoolProp coolfluid((char*)"R134a");
    double cvo, cpo;
    P = 100000;
    T = 293.15;
    coolfluid.CvCp_T_o(T, cvo, cpo);
    cout << coolfluid.Mw << endl;
    prfluid.CvCp_T_o(T, cvo, cpo);
    cout << prfluid.Mw << endl;


    cout << endl << "finished" << endl;


    return 0;
}