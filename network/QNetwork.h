//
// Created by cornelia.blanke on 06.03.2023.
//

#ifndef FLUIDS_QNETWORK_H
#define FLUIDS_QNETWORK_H

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include "../db_connect/DBConnection.h"
#include "Quantum.h"
#include "TBranch.h"
#include "TSegment.h"
#include "TProton.h"
#include "TElectron.h"
//#include "TElectronHP.h"     //circular include!
#include "TNeutron.h"
#include "TCentralPlant.h"


class TCentralPlant;
class TElectron; class TProton; class TNeutron;
class TBranch; class TSegment;


class QNetwork: public Quantum {
    // Fonction non implementée mais pour avoir l'héritage pour l'instant
//    BasicES & performance();
//    float objective(TSequence& seq);

public:
    QNetwork();
    explicit QNetwork(const char* inputfile);           // specify settings file
    QNetwork(const char* case_dir, const char* case_name);        // specify case dir + name

    QNetwork(const QNetwork&) = delete;                     // no copy constructor
    QNetwork&  operator=(const  QNetwork &)  = delete;      // no assignment
    ~QNetwork() override;

    void display() override;

    bool stringtobool(const char* yesno);
    bool stringtobool(char yesno);

    void set_setting(const char* key, const char* value);
    string get_setting(const char* key) const;
    vector<string> parse_setting(const char* key, unsigned int until=10000) const;

    vector<vector<string>> read_csv(const char * inputfile, char sep='\t');
    vector<vector<string>> read_csv_if_exists(const char * inputfile, char sep='\t');
    vector<vector<string>> read_gpkg(const char * inputfile, const char* layerName);        // read gpkg layer and return string table

    void read_settings(const char * inputfile);
    void set_filepath();
    void connect_to_database();
    void get_geopackage_from_database(const char* keyword) const;
    void get_geopackages_from_database() const;
    void read_meteo(const char * inputfile, const char* label_air, int nTln=1, int nTcol=1, char sep='\t');   // compatibility old version (air=sol)
    void read_meteo(const char * inputfile, const char* label_air, const char* label_sol, int nTln=1, int nTcol=1, char sep='\t');
    void read_meteo(const char * tablefile, const char * meteoLayerName,    // read from gpkg
                    const char* label_air, const char* label_sol);
    void read_meteo();

    void read_networktable(const char * inputfile, unsigned int nTln=0, unsigned int nTcol=2, char sep = '\t');
    void read_lengthtable(const char * inputfile, unsigned int nTln=1, unsigned int nTcol=1, char sep='\t');
    void read_diamtable(const char * inputfile, unsigned int nTln=1, unsigned int nTcol=1, char sep='\t');
    // read information from gpkg file
    void get_networktable(const char* inputfile, const char* nodeLayerName="Nodes");
    void get_diam_length_table(const char* inputfile, const char* edgeLayerName="Edges");
    void read_tables(const char* inputfile, const char* nodeLayerName="Nodes", const char* edgeLayerName="Edges");
    void read_tables();

    void read_pipetable(const char * inputfile, unsigned int nTln=2, unsigned int nTcol=1, char sep='\t');
    void read_pipetable(const char * tablefile, const char * pipeLayerName);    // gpkg
    void read_pipetable();

    void read_SST(const char * inputfile, unsigned int nTln=1, char sep='\t');

    void create_electrons_from_gpkg(const char * networkfile, const char* tablefile, const char* sstLayerName="SST",
                                    const char* hpLayerName="HP Catalogue", const char* demandLayerName="Demand Data", const char* dTLayerName="dT Data");
    void create_electrons_from_SST(const char * inputfile, const char * heatpumpfile="", const char * demandfile="",
                                   const char * dTfile="", unsigned int nTln=1, char sep='\t');
    void create_electrons_from_SST();

    template<typename Q> void for_each_element(vector<Q*> List_of_Elements, void(Q::*ptrfunc)());
    template<typename Q, typename T> void for_each_element(vector<Q*> List_of_Elements, void(Q::*ptrfunc)(T), T t);
    void for_each_electron(void(TElectron::*ptrfunc)());
    //void for_each_electron(void(TElectron::*ptrfunc)(bool), bool b);
    template<typename T> void for_each_electron(void(TElectron::*ptrfunc)(T), T t);

    void for_each_proton(void(TProton::*ptrfunc)());
    template<typename T> void for_each_proton(void(TProton::*ptrfunc)(T), T t);
    void for_each_neutron(void(TNeutron::*ptrfunc)());
    template<typename T> void for_each_neutron(void(TNeutron::*ptrfunc)(T), T t);
    void for_each_branch(void(TBranch::*ptrfunc)());
    template<typename T> void for_each_branch(void(TBranch::*ptrfunc)(T), T t);
    void for_each_segment(void(TSegment::*ptrfunc)());
    template<typename T> void for_each_segment(void(TSegment::*ptrfunc)(T), T t);

    void define_electrons(vector<TElectron> &elec);         // link to already existing electrons
    void create_electrons();                        // create new electrons (of type TElectron)
    void define_centralplant(TCentralPlant& central);  // link to already existing central plant
    void create_centralplant();                        // create new central plant
    void define_protons(vector<TProton> &prot);
    void create_protons();
    void create_protons_from_electrons();           // create protons belonging to electrons
    void define_neutrons(vector<TNeutron> &neut);
    void create_neutrons();
    void define_branches(vector<TBranch> &branch);
    void create_branches();
    void define_segments(vector<TSegment> &seg);
    void create_segments();

    void init_electrons(double all_electrons[], string property);   // init electrons from a list of data

    template<typename T> T** allocate_table(unsigned int rows, unsigned int cols, T **M);
    template<typename T> void deallocate_table(unsigned int rows, T **M);
    template<typename T> void deallocate_list(unsigned int len, T **L, bool elements_created);
    template<typename T> void deallocate_ptrvector(vector<T*> &L, bool elements_created);


    template<typename T> T** read_table(unsigned int &rows, unsigned int &cols, T **M,
                                        const char* inputfile, unsigned int nTln, unsigned int nTcol, char sep='\t');
    template<typename T> void read_table(vector<vector<T>> &M,
                                        const char* inputfile, unsigned int nTln, unsigned int nTcol, char sep='\t');
    template<typename T> void display_table(unsigned int rows, unsigned int cols, T **M,
                                            const char*  name_of_table="Table", unsigned int prec=3);
    template<typename T> void display_table(const vector<vector<T>> &M,
                                            const char* name_of_table="Table", unsigned int prec=3);


    void compute_mtable();
    void display_mtable();


    void compute_mu();
    void compute_mu(Quantum* pQuantum);
    void compute_muTable();
    void display_muTable();

    void compute_demand();
    void compute_demand(Quantum* pQuantum);

    void get_list_of_ID();
    void compute_ANZ();
    void compute_m_from_table() override;
    void compute_m() override;
    void compute_mu_from_table() override;

    void compute_massflow_from_mu() override;
    void compute_mu_from_massflow() override;
    void compute_mubar_from_mu() override;

    void compute_k0();

    void get_Toutside(environment &ext);
    void get_Tground(environment &ext);
    void get_Toutside();
    void get_Tground();
    void init_Tambient();               // inits the Tambient pointer
    void init_Tambient(bool below) override;     // should not be used for QNetwork
    void init_Tref () override;         // inits the Tref pointer
    void init_Tsource();                // inits the Tsource pointer

    void connect_network() override;
    void connect_flow() override;

    void compute_deltabar();

    void compute_dP() override;
    void compute_tau() override;

    void compute_COP();

    //estimate_powerloss() estimates dPower without computing the temperatures
    void estimate_powerloss() override;     //sum of all dPower
    void compute_nubar() override;          //loop over all elements
    void compute_epsilon() override;        //loop over all elements
    void compute_power_from_massflow() override;    //loop over all elements except TElectron

    void compute_material_cost() override;                          //sum of all elements
    void compute_engineering_cost(const char* region) override;     //sum of all elements
    void compute_engineering_cost() {compute_engineering_cost((char*) "");}  // call default
    void compute_total_invest() override;                           //loop over all elements
    void compute_operating_cost() override;                         //sum of all elements
    void compute_resources_cost() override;                         //sum of all elements

    // solvers
    void init() override;
    void extract_quantum_network();
    void init_network();
    void init_timestep();
    void solve_timestep();      // use settings from settings file
    void solve_timestep(bool MF, bool T, bool Tsolver, bool COP, bool C);
    void solve_massflow();
    void solve_temperature();
    void solve_temperature(Quantum* pQuantum);
    void next_timestep();

    void save_results_to_file();    // use settings from settings file
    void save_results_to_file(const vector<vector<bool>>& quantities, unsigned int frequency=1e8);
    // 7 element types: Network, CentralPlant, Electrons, Protons, Neutrons, Segments, Branches
    // use function pointers? https://www.tutorialspoint.com/function-pointer-to-member-function-in-cplusplus

    vector<Quantum*> print_headline(ostream &results, const vector<bool>& elements);  //returns print order
    void print_timestep(ostream &results, vector<Quantum*>& printElements, int property) const;
    double get_property_by_idx(Quantum*& q, int& prop_id) const;

    void save_results_to_gpkg();
    void save_results_to_gpkg(const char* gpkgfile, const vector<vector<bool>>& quantities, unsigned int frequency=1e8);

    void save_results();

    void move_geopackage_to_database(const char* keyword, bool delete_local=true) const;
    void clear_tmp(const char* keyword) const;
    void finalize(const bool delete_local=true) const;






public:
    unsigned int nbBranches=0, nbRecords=0;    // table dimensions
    vector<vector<int>> networktable;
    vector<int> List_of_ID;

    vector<vector<double>> lengthtable;
    vector<vector<int>> diamtable;       // nominal diameters DN are always int-values
    vector<vector<double>> pipetable;

    vector<vector<int>> mtable;
    vector<vector<double>> muTable;     // not needed any more



    unsigned int A=0, N=0, Z=0;  // N: number of neutrons, Z: number of protons, A = Z + N
    environment* meteo= nullptr;              // meteo includes T and Tsol

    unsigned int timestep=0, max_timestep=0;
    double k0;          // Kane factor
    double delta_bar;   // network working performance
    double Toutside=273.15, Tground=273.15, T0=273.15;
    double res_cost_per_kWh=0.0;


    TCentralPlant* CentralPlant=nullptr;

    vector<TElectron*> List_of_Electrons;      // length of list is Z
    vector<TProton*> List_of_Protons;          // length of list is Z
    /* should always be created with create_protons_from_electrons() */

    vector<TNeutron*> List_of_Neutrons;        // length of list is N, ascending order
    vector<TBranch*> List_of_Branches;         //length of list is 2N+1, ascending order
    vector<vector<TSegment*>> Table_of_Segments; // same shape as lengthtable/diamtable

private:
    unordered_map<string, string> settings;
    string filepath;
    DBConnection* connection=nullptr;
    GPKGTable* GPKGtable=nullptr;

    bool use_database=false;
    bool use_gpkg=false;
    bool use_postgis=false;

    // output properties
    const vector<const char*> properties =
        {"massflow", "speed", "mu", "mubar", "dP", "dP_m",
         "demand", "power", "dPower", "TIn", "TOut", "dT", "delta", "nubar", "epsilon",
         "material_cost", "engineering_cost", "operating_cost", "resources_cost"};

    bool centralplant_created=false, electrons_created=false, protons_created=false,
        neutrons_created=false, branches_created=false, segments_created=false,
        meteo_created=false;
};








/*** General functions ***/
// allocate 2d array
template<typename T> T** QNetwork::allocate_table(unsigned int rows, unsigned int cols, T **M)
{
    M = new T*[rows];
    for (int i = 0; i < rows; i++){
        M[i] = new T[cols];
    }
    return M;
}

// deallocate 2d array
template<typename T> void QNetwork::deallocate_table(unsigned int rows, T **M)
{
    if (M != NULL) {
        for (int i = 0; i < rows; ++i)
            delete[] M[i];
        delete[] M;
    }
    M = NULL;
}

// deallocate list of pointers and its elements
template<typename T> void QNetwork::deallocate_list(unsigned int len, T **L, bool elements_created)
{
    if (L != NULL) {
        if (elements_created) {
            for (int i = 0; i < len; i++)
                delete L[i];
        }
        delete[] L;
    }
    L = NULL;
}

// deallocate elements of a vector of pointers
template<typename T> void QNetwork::deallocate_ptrvector(vector<T*> &L, bool elements_created)
{
    if (elements_created) {
        for (int i = 0; i < L.size(); i++)
            delete L[i];
    }
    L.clear();
}

// display 2d table
template<typename T> void QNetwork::display_table(unsigned int rows, unsigned int cols, T **M, const char* name_of_table, unsigned int prec)
{
    if (M == NULL){
        cerr << "table not defined" << endl;
        exit(1);
    }
    else {
        // save current precision
        streamsize ss = cout.precision();
        cout << fixed << setprecision((int)prec);
        cout << "Display " << name_of_table << "..." << endl;
        for (int i=0; i<rows; i++) {
            for (int j=0; j<cols; j++) {
                cout << M[i][j] << '\t';
            }
            cout << endl;
        }
        // reset precision to default
        cout << defaultfloat << setprecision((int)ss) << endl;
    }
}

template<typename T> void QNetwork::display_table(const vector<vector<T>> &M, const char* name_of_table, unsigned int prec)
{
    if (M.size() == 0){
        cerr << "table not defined" << endl;
        exit(1);
    }
    else {
        // save current precision
        streamsize ss = cout.precision();
        cout << fixed << setprecision(prec);
        cout << "Display " << name_of_table << "..." << endl;
        for (int i=0; i<M.size(); i++) {
            for (int j=0; j<M[0].size(); j++) {
                cout << M[i][j] << '\t';
            }
            cout << endl;
        }
        // reset precision to default
        cout << defaultfloat << setprecision(ss) << endl;
    }
}



// read table into T **M, find its dimensions (rows x cols) and return M:
template<typename T> T** QNetwork::read_table(unsigned int &rows, unsigned int &cols, T **M,
                                              const char* inputfile, unsigned int nTln, unsigned int nTcol, char sep) {
    // can handle files with nTln=0 or nTcol=0
    unsigned int ln_offset=0, col_offset=0;
    if (nTln*nTcol==0 && !is_arithmetic<T>::value){
        cerr << "Not implemented" << endl;
        exit(1);
    }
    if (nTln == 0)
        ln_offset = 1;
    if (nTcol == 0)
        col_offset = 1;

    TInputFile InputFile((char *) inputfile, nTln+ln_offset, nTcol+col_offset, sep);
    InputFile.open();
    InputFile.config(Tiofile::onrcds);
    TSeqsArray Records = *InputFile.GetRecords();
    TSeqsArray Fields = *InputFile.GetFields();
    InputFile.close();

    // create table
    rows = InputFile.nRecords() + ln_offset;
    cols = InputFile.nFields() + col_offset;
    M = allocate_table(rows, cols, M);

    // read file into table
    if (nTln==0 && nTcol==0) {
        dynstring v = InputFile.firstline();
        int begin, end = 0;
        for (begin = 0; begin < v.size(); begin++)
            if (isdigit(v[begin]) || v[begin] == '+' || v[begin] == '-' || v[begin] == '.') {
                for (end = begin+1; end < v.size(); end++)
                    if (v[end] == sep) {
                        end -= 1;
                        break;
                    }
                break;
            }
        string word;
        for (int i=begin; i<=end; i++)
            word += v[i];
        M[0][0]=stod(word);
    }
    if (ln_offset > 0){
        for (int j = 0; j < InputFile.nFields(); j++)
            M[0][j+col_offset] = stod((Fields.Get(j)->Title)->liste());
    }
    if (col_offset > 0){
        for (int i = 0; i < InputFile.nRecords(); i++)
            M[i+ln_offset][0] = stod((Records.Get(i)->Title)->liste());
    }
    for (unsigned int i = 0; i < InputFile.nRecords(); i++) {
        TSequence *PtrSeq = Get(Records, i+1);
        for (int j = 0; j < InputFile.nFields(); j++) {
            M[i+ln_offset][j+col_offset] = PtrSeq->Get(j);
        }
    }

    return M;
}

template<typename T> void QNetwork::read_table(vector<vector<T>> &M,
                                              const char* inputfile, unsigned int nTln, unsigned int nTcol, char sep) {
    // can handle files with nTln=0 or nTcol=0
    unsigned int ln_offset=0, col_offset=0;
    if (nTln*nTcol==0 && !is_arithmetic<T>::value){
        cerr << "Not implemented" << endl;
        exit(1);
    }
    if (nTln == 0)
        ln_offset = 1;
    if (nTcol == 0)
        col_offset = 1;

    TInputFile InputFile((char *) inputfile, nTln+ln_offset, nTcol+col_offset, sep);
    InputFile.open();
    InputFile.config(Tiofile::onrcds);
    TSeqsArray Records = *InputFile.GetRecords();
    TSeqsArray Fields = *InputFile.GetFields();
    InputFile.close();

    // create table
    M.clear();
    M.resize(InputFile.nRecords()+ln_offset, vector<T>(InputFile.nFields()+col_offset, 0));

    // read file into table
    if (nTln==0 && nTcol==0) {
        dynstring v = InputFile.firstline();
        int begin, end = 0;
        for (begin = 0; begin < v.size(); begin++)
            if (isdigit(v[begin]) || v[begin] == '+' || v[begin] == '-' || v[begin] == '.') {
                for (end = begin+1; end < v.size(); end++)
                    if (v[end] == sep) {
                        end -= 1;
                        break;
                    }
                break;
            }
        string word;
        for (int i=begin; i<=end; i++)
            word += v[i];
        M[0][0] = stod(word);
    }
    if (ln_offset > 0){
        for (int j = 0; j < InputFile.nFields(); j++) {
            string s;
            s = Fields.Get(j)->Title->liste();
            M[0][j + col_offset] = stod(s);
        }
    }
    if (col_offset > 0){
        for (int i = 0; i < InputFile.nRecords(); i++) {
            string s;
            s = Records.Get(i)->Title->liste();
            M[i + ln_offset][0] = stod(s);
        }
    }
    for (unsigned int i = 0; i < InputFile.nRecords(); i++) {
        TSequence *PtrSeq = Get(Records, i+1);
        for (int j = 0; j < InputFile.nFields(); j++) {
            float f;
            f = PtrSeq->Get(j);
            M[i+ln_offset][j+col_offset] = f;
        }
    }
}


#endif //FLUIDS_QNETWORK_H
