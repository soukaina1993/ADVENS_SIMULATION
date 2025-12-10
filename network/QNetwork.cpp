//
// Created by cornelia.blanke on 06.03.2023.
//

#include <regex>
#include <functional>
#include <filesystem>
#include "QNetwork.h"
#include "TElectronHP.h"
#include "../geodata/GeoDataset.h"


using namespace std;

QNetwork::QNetwork() {}

QNetwork::QNetwork(const char * inputfile) {
    read_settings(inputfile);
}

QNetwork::QNetwork(const char* case_dir, const char* case_name) {
    string file;
    file = case_dir;
    file += "/";
    file += case_name;
    file += "_settings.txt";
    read_settings((char *) file.c_str());
}

QNetwork::~QNetwork() {
    deallocate_ptrvector(List_of_Electrons, electrons_created);
    deallocate_ptrvector(List_of_Protons, protons_created);
    deallocate_ptrvector(List_of_Neutrons, neutrons_created);
    deallocate_ptrvector(List_of_Branches, branches_created);
    if (segments_created)
        for (int m = 0; m < lengthtable.size(); m++)      // loop over table
            for (int n = 0; n < nbBranches; n++) {
                delete Table_of_Segments[m][n];
            }
    if (centralplant_created) delete CentralPlant;
    if (meteo_created) delete meteo;
    delete connection;
    delete GPKGtable;
}

void QNetwork::display() {
    if (networktable.empty()){
        cerr << "Network Table not defined" << endl;
        exit(1);
    }
    else {
        display_table(networktable, "Network Table");
        cout << "Z = " << Z << endl;
        cout << "N = " << N << endl;
        cout << "A = " << A << endl;
    }
}

/*** read functions ***/
bool QNetwork::stringtobool(const char* yesno){
    if (strcmp(yesno, "Y")==0 || strcmp(yesno, "y")==0)
        return true;
    else if (strcmp(yesno, "N")==0 || strcmp(yesno, "n")==0)
        return false;
    else {
        cerr << "Specify Y or N" << endl;
        exit(1);
    }
}

bool QNetwork::stringtobool(const char yesno){
    if (yesno == 'Y' || yesno == 'y')
        return true;
    else if (yesno == 'N' || yesno == 'n')
        return false;
    else {
        cerr << "Specify Y or N" << endl;
        exit(1);
    }
}

void QNetwork::set_setting(const char* key, const char* value){
    string s1, s2;
    s1 = key;
    s2 = value;
    settings[s1] = s2;
}

string QNetwork::get_setting(const char* key) const {
    string s;
    s = key;
    return settings.at(s);      // throws std::out_of_range exception if key is not available
}

vector<string> QNetwork::parse_setting(const char* key, unsigned int until) const {
    stringstream ss(get_setting(key));
    string word;
    vector<string> result;

    if (getline(ss, word, '[')) {          // get first part
        result.push_back(word);
        until--;
        while (until >= 0 && getline(ss, word, '/')) {    // get second part
            result.push_back(word);
            until--;
        }
        if (result.back().back() == ']')    // remove ']' from last entry
            result.back().pop_back();
    }
    else {
        cerr << "Missing key " << key << " in settings file" << endl;
        exit(1);
    }
    return result;
}

// read_csv function found on internet
vector<vector<string>> QNetwork::read_csv(const char * inputfile, char sep) {
    vector<vector<string>> content = read_csv_if_exists(inputfile, sep);
    if (content.empty()) {
        cerr << "Could not open the file " << inputfile << "\n";
        exit(1);
    }
    return content;
}

vector<vector<string>> QNetwork::read_csv_if_exists(const char * inputfile, char sep) {
    vector<vector<string>> content;
    fstream file(inputfile, ios::in);
    if (file.is_open()) {
        vector<string> row;
        string line, word;
        while (getline(file, line)) {
            row.clear();
            stringstream str(line);
            while (getline(str, word, sep))
                row.push_back(word);
            content.push_back(row);
        }
        file.close();
    }
    return content;
}

vector<vector<string>> QNetwork::read_gpkg(const char * inputfile, const char* layerName){
    vector<vector<string>> content;
    GeoDataset readData (inputfile, GeoDataset::Accessmode::READ);
    content = readData.get_layer_as_table(layerName);

    return content;
}

void QNetwork::read_settings(const char * inputfile){
    vector<const char*> keys =
            {"Case Directory", "Case Name", "Connection", "GPKG Table",
             "Network", "Q-Network", "Data", "Results",
             "Meteo File", "Meteo Year", "Air Temp", "Ground Temp",
             "Fluid", "Boiler Type", "Boiler Power", "Boiler Price/kWh",
             "Supply Temperature", "Estimated Return Temperature",
             "Solve Mass Flow", "Solve Temperature", "Solver", "Update COP", "Solve Cost",
             "Save Results", "Save Frequency"};
    // append properties to keys
    keys.insert(keys.end(), properties.begin(), properties.end());
    string line;
    fstream file(inputfile, ios::in);
    if (file.is_open()) {
        //vector<const char*> patterns = keys;
        set<const char*> patterns(keys.begin(), keys.end());
        while (getline(file, line)) {
            for (auto & pattern : patterns) {
                // check if beginning of line matches key plus a whitespace (space or tab):
                int i=0;
                bool found=true;
                while (pattern[i] != '\0'){
                    if (pattern[i] != line[i]) {
                        found = false;
                        break;
                    }
                    i++;
                }
                if (found && (line[i] == ' ' || line[i] == '\t') ){
                    string substring;
                    bool details=false;     // details inside [...]
                    for (auto c : line.substr(i+1, line.size())) {
                        if (details || (c != '\t' && c !=' '))   // ignore whitespaces
                            substring.push_back(c);
                        if (c == '[')
                            details = true;
                    }
                    set_setting(pattern, substring.c_str());
                    //patterns.erase(std::remove(patterns.begin(), patterns.end(), pattern), patterns.end());
                    patterns.erase(pattern);
                    break;
                }
            }
        }

// check if all necessary keys are set
        if (patterns.find("Network") == patterns.end() && patterns.find("Q-Network") == patterns.end() &&
                patterns.find("Data") == patterns.end() && patterns.find("Results") == patterns.end()) {
            // new version using .gpkg files
            use_gpkg = true;
            patterns.erase("Meteo File");        // erase unnecessary key

            if (patterns.find("Connection") == patterns.end())
                use_database = true;
            else
                for (auto &key : {"Case Name", "Connection", "GPKG Table"})   // keys not used
                    patterns.erase(key);
        }
        else if (patterns.find("Case Name") == patterns.end() && patterns.find("Meteo File") == patterns.end()) {
            // old version using .csv files
            for (auto &key : {"Connection", "GPKG Table", "Network", "Q-Network", "Data", "Results"})   // keys not used
                patterns.erase(key);
        }
        else {
            cerr << R"(Specify either "Case Name" and "Meteo File" or all GPKG file names)" << endl;
            exit(1);
        }

        if (patterns.find("Meteo Year") == patterns.end()) {    // for compatibility with old file style
            for (auto &key : {"Air Temp", "Ground Temp"}) {
                set_setting(key, get_setting("Meteo Year").c_str());
                patterns.erase(key);
            }
        }
        else {
            patterns.erase("Meteo Year");
        }

        // now if everything is there, patterns should be empty
        if (!patterns.empty()) {
            for (const auto& pat : patterns)
                cerr << "\"" << pat << "\" ";
            cerr << "not specified in settings file" << endl;
            exit(1);
        }
    }
    else {
        cerr << "Could not open the file " << inputfile << "\n";
        exit(1);
    }
    file.close();
}

void QNetwork::set_filepath() {
    filepath = get_setting("Case Directory");
    filesystem::create_directory(filepath);     // create dir if it does not exist
    filepath += "/";
    if (!use_gpkg)
        filepath += get_setting("Case Name");
}

void QNetwork::connect_to_database(){
    // get infos from settings
    string dbname, dbuser, dbpasswd;
    stringstream ss(get_setting("Connection"));


    if (getline(ss, dbname, '['))               // get first part
        if (getline(ss, dbuser, '/'))           // get first part of second part
            if(getline(ss, dbpasswd, ']'))      // get second part of second part
                ;

    if (dbname.empty() || dbuser.empty() || dbpasswd.empty()) {
        cerr << "Define database connection in settings file" << endl;
        exit(1);
    }

    auto setting = parse_setting("GPKG Table");
    string id="id", name="name", data="data", time="time";

    if (setting.size() > 1){
        if (setting[1] != "-")
            id = setting[1];
        if (setting.size() > 2){
            if (setting[2] != "-")
                name = setting[2];
            if (setting.size() > 3){
                if (setting[3] != "-")
                    data = setting[3];
                if (setting.size() > 4 && setting[4] != "-")
                    time = setting[4];
            }
        }
    }

    // set members
    connection = new DBConnection (dbname.c_str(), dbuser.c_str(), dbpasswd.c_str());

    GPKGtable = new GPKGTable (connection, get_setting("Case Name").c_str(),
                        setting[0].c_str(), id.c_str(), name.c_str(), data.c_str(),
                        false, true, time.c_str());
}

void QNetwork::get_geopackage_from_database(const char* keyword) const {
    string file;
    stringstream ss1(get_setting(keyword));
    getline(ss1, file, '[');
    file += ".gpkg";

    if (! GPKGtable->download_gpkg(file.c_str(), (filepath + file).c_str())) {
        cerr << "Geopackage " << file << " could not be downloaded from database" << endl;
        exit(1);
    }
}

void QNetwork::get_geopackages_from_database() const {
    get_geopackage_from_database("Network");
    get_geopackage_from_database("Data");
}

void QNetwork::read_meteo(const char * inputfile, const char* label_air, const int nTln, const int nTcol, const char sep) {
    read_meteo(inputfile, label_air, label_air, nTln, nTcol, sep);
}

void QNetwork::read_meteo(const char * inputfile, const char* label_air, const char* label_sol, const int nTln, const int nTcol, const char sep) {
    meteo = new environment;
    meteo_created = true;
    //TODO Getdata() is very slow!    meteo->init(inputfile, meteo->T, "2019", sep, nTln, nTcol);

    meteo->init(inputfile, sep, nTln, nTcol);   // inits private data
    vector<vector<string>> content = read_csv(inputfile);
    meteo->Time = (int)content.size()-nTln;
    meteo->T.clear();
    meteo->T.resize(content.size()-nTln);
    meteo->Tsol.clear();
    meteo->Tsol.resize(content.size()-nTln);

    // lambda function working for T and Tsol
    auto fill_T = [](vector<double>& T, const int j, const int nTln, const vector<vector<string>>& content)
    {
        bool celsius = false, kelvin = false;
        for (int i = nTln; i < content.size(); i++)
            if (celsius)
                T[i - nTln] = stod(content[i][j]) + 273.15;
            else if (kelvin)
                T[i - nTln] = stod(content[i][j]);
            else if (stod(content[i][j]) < 100) {   // in first time-step, set to C or K
                celsius = true;
                T[i - nTln] = stod(content[i][j]) + 273.15;
            } else {
                kelvin = true;
                T[i - nTln] = stod(content[i][j]);
            }
    };

    for (int j=nTcol; j<content[0].size(); j++) {
        if (content[0][j] == label_air) {
            fill_T(meteo->T, j, nTln, content);
            if (strcmp(label_air,label_sol)==0){      // don't read twice
                meteo->Tsol = meteo->T;
                break;
            }
        }
        else if (strcmp(content[0][j].c_str(), label_sol)==0) {
            fill_T(meteo->Tsol, j, nTln, content);
            if (label_air==label_sol){      // don't read twice
                meteo->T = meteo->Tsol;
                break;
            }
        }
    }
    meteo->Tmean = sum(meteo->T)/meteo->Time;
}

void QNetwork::read_meteo(const char * tablefile, const char* meteoLayerName, const char* label_air, const char* label_sol) {
    meteo = new environment;
    meteo->init(tablefile);   // inits private data
    meteo_created = true;
    GeoDataset tabledata(tablefile, GeoDataset::Accessmode::READ);
    OGRLayer* meteoLayer = tabledata.get_Layer(meteoLayerName);

    meteo->Time = (int) meteoLayer->GetFeatureCount();
    meteo->T.clear();
    meteo->T.reserve(meteo->Time);
    meteo->Tsol.clear();
    meteo->Tsol.reserve(meteo->Time);

    // lambda function working for T and Tsol
    auto fill_T = [](vector<double>& T, const char* label, OGRLayer*& meteoLayer)
    {
        bool celsius = false, kelvin = false;
        for (auto &data : meteoLayer)
            if (celsius)
                T.push_back(data->GetFieldAsDouble(label) + 273.15);
            else if (kelvin)
                T.push_back(data->GetFieldAsDouble(label));
            else if (data->GetFieldAsDouble(label) < 100) {   // in first time-step, set to C or K
                celsius = true;
                T.push_back(data->GetFieldAsDouble(label) + 273.15);
            } else {
                kelvin = true;
                T.push_back(data->GetFieldAsDouble(label));
            }
    };

    fill_T(meteo->T, label_air, meteoLayer);
    if (strcmp(label_air,label_sol)==0)      // don't read twice
        meteo->Tsol = meteo->T;
    else
        fill_T(meteo->Tsol, label_sol, meteoLayer);

    meteo->Tmean = sum(meteo->T)/meteo->Time;
}

void QNetwork::read_meteo() {
    if (use_gpkg){
        string layer;
        auto setting = parse_setting("Data", 1);

        if (setting.size() > 1 && setting[1] != "-")
            layer = setting[1];
        else
            layer = "Meteo Data";                      // set default

        read_meteo((filepath + setting[0] + ".gpkg").c_str(), layer.c_str(),
                   get_setting("Air Temp").c_str(),
                   get_setting("Ground Temp").c_str());
    }
    else
        read_meteo(get_setting("Meteo File").c_str(),
                   get_setting("Air Temp").c_str(),
                   get_setting("Ground Temp").c_str());
}

void QNetwork::read_networktable(const char * inputfile, unsigned int nTln, unsigned int nTcol, char sep){
    read_table(networktable, inputfile, nTln, nTcol, sep);
    nbRecords = networktable.size();
    nbBranches = networktable[0].size();
}

void QNetwork::read_lengthtable(const char * inputfile, unsigned int nTln, unsigned int nTcol, char sep){
    read_table(lengthtable, inputfile, nTln, nTcol, sep);
    if (!(lengthtable.size() == nbRecords - 2 || lengthtable.size() == nbRecords - 1) || lengthtable[0].size() != nbBranches){
        cerr << "Length table has wrong dimensions" << endl;
        exit(1);
    }
}

void QNetwork::read_diamtable(const char * inputfile, unsigned int nTln, unsigned int nTcol, char sep){
    read_table(diamtable, inputfile, nTln, nTcol, sep);
    if (!(diamtable.size() == nbRecords - 2 || diamtable.size() == nbRecords - 1) || diamtable[0].size() != nbBranches){
        cerr << "Diameter table has wrong dimensions" << endl;
        exit(1);
    }
}

/*** networktable functions ***/
void QNetwork::get_list_of_ID(){
    for (int m=1; m<nbRecords-1; m++)
        for (int n=0; n<nbBranches; n++)
            if (networktable[m][n] > 0)
                List_of_ID.push_back(networktable[m][n]);
    sort(List_of_ID.begin(), List_of_ID.end());
}

void QNetwork::compute_ANZ() {
    if (networktable.empty()) {
        cout << "Error: networktable not defined" << endl;
        exit(1);
    }
    N = (nbBranches-1)/2;
    if (List_of_ID.empty())
        get_list_of_ID();
    Z = List_of_ID.size();
    A = N + Z;
}

/*** pipetable functions ***/
void QNetwork::read_pipetable(const char* inputfile, unsigned int nTln, unsigned int nTcol, char sep){
    read_table(pipetable, inputfile, nTln, nTcol, sep);
    // convert units to SI
    if (pipetable[0].size()==8) {    //two-layer pipe
        for (int row = 0; row < pipetable.size(); row++)
            for (int col = 1; col < pipetable[0].size(); col++) {
                pipetable[row][col] /= 1000.;
                if (col == 3)
                    col = 6;  //jump over columns
            }
    }
    else if (pipetable[0].size()==10){    //three-layer pipe
        for (int row = 0; row < pipetable.size(); row++)
            for (int col = 1; col < pipetable[0].size(); col++) {
                pipetable[row][col] /= 1000.;
                if (col == 4)
                    col = 8;  //jump over columns
            }
    }
    else{
        cerr << "pipetable has wrong dimensions" << endl;
        exit(1);
    }
}

void QNetwork::read_pipetable(const char * tablefile, const char * pipeLayerName){
    GeoDataset tabledata(tablefile, GeoDataset::Accessmode::READ);
    OGRLayer* pipeLayer = tabledata.get_Layer(pipeLayerName);
    vector<const char*> fields;
    int fieldcount=0;
    pipetable.clear();

    if (pipeLayer != nullptr){
        if ((fieldcount=pipeLayer-> GetLayerDefn()->GetFieldCount()) == 8)
            fields = {"DN", "D_int", "D_ext", "D_ext_inso", "Conduct_steel", "Conduct_inso", "U", "Roughness"};
        else if (fieldcount == 10){
            fields = {"DN", "D_int", "D_ext", "D_ext_inso", "D_ext_mant", "Conduct_steel", "Conduct_inso", "Conduct_mant", "U", "Roughness"};
        }
        else{
            cerr << "Pipe info incomplete" << endl;
            exit(1);
        }

        for (auto &feature : pipeLayer) {
            vector<double> pipeinfo;
            for (auto &field: fields) {
                double value = feature->GetFieldAsDouble(field);
                // scale by 0.001
                if (strcmp(field, "DN")!=0 && strncmp(field, "Conduct_", 8)!=0 && strcmp(field, "U")!=0)
                    value /= 1000.0;
                pipeinfo.push_back(value);
            }
            pipetable.push_back(pipeinfo);
        }
    }
    else {
        cerr << "pipetable cannot be read" << endl;
        exit(1);
    }
}

void QNetwork::read_pipetable(){
    if (use_gpkg) {
        string layer;
        auto setting = parse_setting("Data", 2);

        if (setting.size() > 2 && setting[2] != "-")
            layer = setting[2];
        else
            layer = "Pipe Catalogue";       // set default

        read_pipetable((filepath + setting[0] + ".gpkg").c_str(), layer.c_str());
    }
    else
        read_pipetable((filepath + "_pipes.csv").c_str());
}

/*** read list of substations ***/
void QNetwork::read_SST(const char * inputfile, unsigned int nTln, char sep){
    vector<vector<string>> content = read_csv(inputfile, sep);
    if (List_of_Electrons.empty() || content.size()-nTln < Z) {
        cerr << "Size of SST file does not match" << endl;
        exit(1);
    }
    for(unsigned int i=0; i<Z; i++){      // loop over existing List_of_Electrons
        for(unsigned int j=nTln; j<content.size(); j++) {
            if (List_of_Electrons[i]->ID == stoi(content[j][0])){   // find corresponding entry on file
                List_of_Electrons[i]->D =
                        District(1000*stod(content[j][1]),
                                 stod(content[j][8]),
                                 content[j][4].c_str(),
                                 content[j][5].c_str());
                List_of_Electrons[i]->altitude = stod(content[j][6]);
                List_of_Electrons[i]->dT = stod(content[j][9]);
                List_of_Electrons[i]->length = stod(content[j][10]);
                List_of_Electrons[i]->DN = stoi(content[j][11]);
                break;
            }
        }
    }
}

///*** get information from gpkg ***/
void QNetwork::get_networktable(const char* inputfile, const char* nodeLayerName){
    GeoWriteDataset dataset (inputfile, GeoDataset::Accessmode::READ, nodeLayerName, "", "");
    dataset.extract_networktable(networktable);
}

void QNetwork::get_diam_length_table(const char* inputfile, const char* edgeLayerName){
    GeoWriteDataset dataset (inputfile, GeoDataset::Accessmode::READ, "", edgeLayerName, "");
    dataset.extract_diam_length_table(diamtable, lengthtable);
}

void QNetwork::read_tables(const char* inputfile, const char* nodeLayerName, const char* edgeLayerName){
    GeoWriteDataset dataset (inputfile, GeoDataset::Accessmode::READ, nodeLayerName, edgeLayerName, "");
    dataset.extract_networktable(networktable);
    dataset.extract_diam_length_table(diamtable, lengthtable);
    nbRecords = networktable.size();
    nbBranches = networktable[0].size();
}

void QNetwork::read_tables() {
    if (use_gpkg) {
        string nodelayer="Nodes", edgelayer="Edges";
        auto setting = parse_setting("Q-Network", 2);
        string table = setting[0] + ".gpkg";

        if (setting.size() > 1) {    // get first part of second part (if exists)
            if (setting[1] != "-")
                nodelayer = setting[1];
            if (setting.size() > 2 && setting[2] != "-")    // get second part of second part (if exists)
                edgelayer = setting[2];
        }
        read_tables((filepath+table).c_str(), nodelayer.c_str(), edgelayer.c_str());
    }
    else {
        read_networktable((filepath + ".csv").c_str());
        read_diamtable((filepath + "_diam.csv").c_str());
        read_lengthtable((filepath + "_length.csv").c_str());
    }
}




void QNetwork::create_electrons_from_gpkg(const char* networkfile, const char* tablefile, const char* sstLayerName,
                                          const char* hpLayerName, const char* demandLayerName, const char* dTLayerName){

    OGRLayer* hpLayer=nullptr;
    OGRLayer* demandLayer=nullptr;
    OGRLayer* dTLayer=nullptr;

    GeoDataset tabledata(tablefile, GeoDataset::Accessmode::READ);
    GeoWriteDataset netdata(networkfile, GeoDataset::Accessmode::READ, "", "", sstLayerName);
    netdata.set_Layers();

    if (netdata.sstLayer->GetFeatureCount() < Z) {
        cerr << "Size of SST layer does not match" << endl;
        exit(1);
    }
    electrons_created=true;
    List_of_Electrons.clear();
    List_of_Electrons.reserve(Z);

    OGRFeature* feature;
    for(auto &id : List_of_ID) {      // loop over existing List_of_ID
        if ((feature=netdata.sstLayer->GetFeature(id)) == nullptr){
            cerr << "No information about Substation " << id << endl;
            exit(1);
        }

        string s;
        s = to_string(id);
        const char* sstid = s.c_str();

        // create electrons
        TElectron* pElectron;
        if (strcmp(feature->GetFieldAsString("Substation Type"), "HP") == 0) {
            auto* pelectronHP = new TElectronHP(id);

            if (hpLayer == nullptr)     // set only once
                hpLayer = tabledata.get_Layer(hpLayerName);

            // search Heatpump on file
            OGRFeature* hp;
            while ((hp = hpLayer->GetNextFeature()) != nullptr)
                if (strcmp(feature->GetFieldAsString("Substation Name"), hp->GetFieldAsString("Name")) == 0)
                    break;      // heatpump found on file

            if (hp!=nullptr && strcmp(hp->GetFieldAsString("COP"), "Default")==0) {   // heatpump found, compute COP
                pelectronHP->create_heatpump(hp->GetFieldAsString("Fluid"),
                                             hp->GetFieldAsDouble("eta_comp"),
                                             hp->GetFieldAsDouble("pinchCond"),
                                             hp->GetFieldAsDouble("pinchEvap"));
            }
            else if (hp!=nullptr && strcmp(hp->GetFieldAsString("COP"), "Data")==0) {
                cerr << "Option Data not implemented" << endl;
                exit(1);
            }
            else if (hp != nullptr) {   // heatpump found, read COP from file
                pelectronHP->create_heatpump();
                pelectronHP->set_COP(hp->GetFieldAsDouble("COP"));   // set constant COP
            }
            else { // no file or heatpump not found on file
                pelectronHP->create_heatpump();
            }
            pElectron = pelectronHP;
        }
        else    // "Substation Type"=="HX"
            pElectron = new TElectron (id);
        pElectron->QNet = this;
        List_of_Electrons.push_back(pElectron);

        // fill data
        if (strcmp(feature->GetFieldAsString("Energy [kWh/an]"), "Data") == 0) {
            // read energy data from file
            if (demandLayer == nullptr)     // set only once
                demandLayer = tabledata.get_Layer(demandLayerName);

            pElectron->D =
                    District(-1,       // undefined energy / not necessary here
                             feature->GetFieldAsDouble("T_ref (C)"),
                             feature->GetFieldAsDouble("T_need (C)"),
                             Building::get_type(feature->GetFieldAsString("Building Type")),
                             feature->GetFieldAsString("Heat Type"));
            pElectron->D.power.clear();
            for (auto &data : demandLayer)         // fill data into D.power
                pElectron->D.power.push_back(data->GetFieldAsDouble(sstid));
            pElectron->D.Q_total = sum(pElectron->D.power);
            pElectron->D.B[0].energy = pElectron->D.Q_total;  // update building
            pElectron->D.power_max = max(pElectron->D.power);
        }
        else {                                  // use SIA norm
            pElectron->D =
                    District(1000 * feature->GetFieldAsDouble("Energy [kWh/an]"),
                             feature->GetFieldAsDouble("T_ref (C)"),
                             feature->GetFieldAsDouble("T_need (C)"),
                             Building::get_type(feature->GetFieldAsString("Building Type")),
                             feature->GetFieldAsString("Heat Type"));
        }
        pElectron->altitude = feature->GetFieldAsDouble("Altitude [m]");
        if (strcmp(feature->GetFieldAsString("dT"), "Data") == 0) {
            if (dTLayer == nullptr)
                dTLayer = tabledata.get_Layer(dTLayerName);

            pElectron->dTvar.clear();
            for (auto &data : dTLayer)         // fill data into dTvar
                pElectron->dTvar.push_back(data->GetFieldAsDouble(sstid));
            pElectron->dT = pElectron->dTvar[0];
        }
        else
            pElectron->dT = feature->GetFieldAsDouble("dT");
        pElectron->length = feature->GetFieldAsDouble("L_raccord (m)");
        pElectron->DN = feature->GetFieldAsInteger("DN_raccord");
    }
}

void QNetwork::create_electrons_from_SST(const char * inputfile,const char * heatpumpfile,
                                         const char * demandfile, const char * dTfile,
                                         const unsigned int nTln, const char sep){
    vector<vector<string>> content_SST = read_csv(inputfile, sep);
    vector<vector<string>> content_demand;      // for time-dependent demand data
    vector<vector<string>> content_dT;          // for time-dependent dT data
    vector<vector<string>> content_HP;          // heatpump file

    unsigned int i,j,k,m,n;     // iterators


    if (content_SST.size() - nTln < Z) {
        cerr << "Size of SST file does not match" << endl;
        exit(1);
    }
    electrons_created=true;
    List_of_Electrons.clear();
    List_of_Electrons.reserve(Z);

    for(i=0; i<Z; i++){      // loop over existing List_of_ID
        for(j=nTln; j < content_SST.size(); j++) {
            if (List_of_ID[i] == stoi(content_SST[j][0])){   // find corresponding entry on file
                // create electrons
                TElectron* pElectron;

                if (content_SST[j][2] == "HP") {
                    auto* pelectronHP = new TElectronHP(List_of_ID[i]);

                    if (content_HP.empty())     // read only once
                        content_HP = read_csv_if_exists(heatpumpfile);

                    // search Heatpump on file
                    k=0;
                    for (k=1; k<content_HP.size(); k++)
                        if (content_HP[k][0] == content_SST[j][3])
                            break;      // heatpump found on file

                    if (k>0 && k<content_HP.size() && content_HP[k][2] == "Default") {   // heatpump found, compute COP
                        pelectronHP->create_heatpump(content_HP[k][1],
                                                     stod(content_HP[k][3]),
                                                     stod(content_HP[k][4]),
                                                     stod(content_HP[k][5]));
                    }
                    else if (k>0 && k<content_HP.size() && content_HP[k][2] == "Data") {
                        cerr << "Option Data not implemented" << endl;
                        exit(1);
                    }
                    else if (k>0 && k<content_HP.size()) {   // heatpump found, read COP from file
                        pelectronHP->create_heatpump();
                        pelectronHP->set_COP(stod(content_HP[k][2]));   // set constant COP
                    }
                    else { // no file or heatpump not found on file
                        pelectronHP->create_heatpump();
                    }
                    pElectron = pelectronHP;
                }

                else    // content_SST[j][2]=="HX"
                    pElectron = new TElectron (List_of_ID[i]);
                pElectron->QNet = this;
                List_of_Electrons.push_back(pElectron);

                // fill data
                if (content_SST[j][1] == "Data") {      // read energy data from file
                    if (content_demand.empty())     // read only once
                        content_demand = read_csv_if_exists(demandfile);

                    for (n=1; n<content_demand[0].size(); n++)      // identify column
                        if (content_demand[0][n] == content_SST[j][0])
                            break;
                    List_of_Electrons[i]->D =
                            District(-1,       // undefined energy / not necessary here
                                     stod(content_SST[j][8]),
                                     stod(content_SST[j][7]),
                                     content_SST[j][4].c_str(),
                                     content_SST[j][5].c_str());
                    List_of_Electrons[i]->D.power.clear();
                    for (m=1; m<content_demand.size(); m++)         // fill data into D.power
                        List_of_Electrons[i]->D.power.push_back(stod(content_demand[m][n]));
                    List_of_Electrons[i]->D.Q_total = sum(List_of_Electrons[i]->D.power);
                    List_of_Electrons[i]->D.B[0].energy = List_of_Electrons[i]->D.Q_total;  // update building
                    List_of_Electrons[i]->D.power_max = max(List_of_Electrons[i]->D.power);
                }
                else {                                  // use SIA norm
                    List_of_Electrons[i]->D =
                            District(1000 * stod(content_SST[j][1]),
                                     stod(content_SST[j][8]),
                                     stod(content_SST[j][7]),
                                     content_SST[j][4].c_str(),
                                     content_SST[j][5].c_str());
                }
                List_of_Electrons[i]->altitude = stod(content_SST[j][6]);
                if (content_SST[j][9] == "Data"){
                    if (content_dT.empty())
                        content_dT = read_csv_if_exists(dTfile);
                    for (n=1; n<content_dT[0].size(); n++)      // identify column
                        if (content_dT[0][n] == content_SST[j][0])
                            break;
                    List_of_Electrons[i]->dTvar.clear();
                    for (m=1; m<content_dT.size(); m++)         // fill data into dTvar
                        List_of_Electrons[i]->dTvar.push_back(stod(content_dT[m][n]));
                    List_of_Electrons[i]->dT = List_of_Electrons[i]->dTvar[0];
                }
                else
                    List_of_Electrons[i]->dT = stod(content_SST[j][9]);
                List_of_Electrons[i]->length = stod(content_SST[j][10]);
                List_of_Electrons[i]->DN = stoi(content_SST[j][11]);
                break;
            }
            else if (j == content_SST.size() - 1){
                cerr << "No information about Substation " << List_of_ID[i] << endl;
                exit(1);
            }
        }
    }
}

void QNetwork::create_electrons_from_SST(){
    if (use_gpkg) {
        string qnettable, sstlayer;
        // find SST layer
        auto setting = parse_setting("Q-Network");
        qnettable = setting[0] + ".gpkg";
        if (setting.size() > 3 && setting[3] != "-")
                sstlayer = setting[3];
        else
            sstlayer = "SST";       // set default

        // find other layers
        string datatable, hplayer = "HP Catalogue", demlayer = "Demand Data", dtlayer = "dT Data";
        setting.clear();
        setting = parse_setting("Data");
        datatable = setting[0] + ".gpkg";

        if (setting.size() > 3) {    // modify defaults
            if (setting[3] != "-")
                hplayer = setting[3];
            if (setting.size() > 4) {
                if (setting[4] != "-")
                    demlayer = setting[4];
                if (setting.size() > 5 && setting[5] != "-")
                    dtlayer = setting[5];
            }
        }
        create_electrons_from_gpkg((filepath+qnettable).c_str(), (filepath+datatable).c_str(),sstlayer.c_str(),
                                   hplayer.c_str(), demlayer.c_str(), dtlayer.c_str());
    }
    else {
        create_electrons_from_SST((filepath + "_SST.csv").c_str(),
                                  (filepath + "_HP.csv").c_str(),
                                  (filepath + "_SST_demand.csv").c_str(),
                                  (filepath + "_SST_dT.csv").c_str());
    }
}

/*** mtable functions ***/
void QNetwork::compute_mtable() {
    // create mtable
        // first line: number of sst on a branch
        // last line: number of sst below neutrons
        // between is number of sst below protons (connected electron + all below)
    mtable.clear();
    mtable.resize(nbRecords, vector<int>(nbBranches));
    int branch0 = (int)(nbBranches/2);

    for (int j=0; j<nbBranches; j++) {      // initialize first+last line = 0
        mtable[0][j] = 0;
        mtable[nbRecords - 1][j] = 0;
    }
    for (int j=0; j<nbBranches; j++)      // compute first line: sst on branch
        for (unsigned int i = nbRecords - 2; i > 0; i--) {
            if (networktable[i][j] > 0)
                (mtable[0][j])++;
        }

    for (int j=branch0; j>=0; j--) {
        int pos = branch0 + j;
        int neg = branch0 - j;
        // protons
        for (unsigned int i = nbRecords - 2; i > 0; i--) {
            if (networktable[i][pos] > 0)
                if (mtable[i+1][pos] > 0)
                    mtable[i][pos] = mtable[i+1][pos] + 1;
                else
                    mtable[i][pos] = mtable[nbRecords-1][pos] + 1;
            else
                mtable[i][pos] = 0;
            if (networktable[i][neg] > 0)
                if (mtable[i+1][neg] > 0)
                    mtable[i][neg] = mtable[i+1][neg] + 1;
                else
                    mtable[i][neg] = mtable[nbRecords-1][neg] + 1;
            else
                mtable[i][neg] = 0;
        }
        // neutron
        for (int k=0; k<nbBranches; k++){
            if( networktable[nbRecords-1][k] == abs(j)  && j!=0){
                mtable[nbRecords - 1][k] = max(mtable[1][pos], mtable[nbRecords-1][pos]) +
                                            max(mtable[1][neg], mtable[nbRecords-1][neg]);
                break;
            }
        }

    }
}
void QNetwork::display_mtable(){
    display_table(mtable, "mtable");
}

void QNetwork::compute_m_from_table() {
    if (mtable.empty()){
        cerr << "mtable not defined" << endl;
        exit(1);
    }
    else
        m = Z;
}

void QNetwork::compute_m() {
    m = Z;
}

/*** mu functions ***/
void QNetwork::compute_mu(){
    compute_mu(CentralPlant);
    mu = CentralPlant->mu;
}

void QNetwork::compute_mu(Quantum* pQuantum){
    if (!pQuantum->next.empty()) {     // does nothing for TElectron
        pQuantum->mu = 0.0;
        for (int i = 0; i < pQuantum->next.size(); i++) {
            compute_mu(pQuantum->next[i]);
            pQuantum->mu += pQuantum->next[i]->mu;
        }
    }
}


void QNetwork::compute_muTable(){
    muTable.clear();
    muTable.resize(nbRecords, vector<double>(nbBranches));
    int branch0 = (int)(nbBranches/2);

    for (int j=0; j<nbBranches; j++) {      // initialize first+last line = 0
        muTable[0][j] = 0;
        muTable[nbRecords - 1][j] = 0;
    }
    // electrons
    for (int i=0; i<Z; i++){
        unsigned int k_ = List_of_Electrons[i]->k;
        int l_ = List_of_Electrons[i]->l;
        muTable[k_][l_+branch0] = List_of_Electrons[i]->mu;
    }
    // branches and neutrons
    for (int n=N; n>0; n--){
        for (int i=1; i<nbRecords; i++){
            muTable[0][branch0 + n] += muTable[i][branch0 + n];     // positive branch
            muTable[0][branch0 - n] += muTable[i][branch0 - n];     // negative branch

            for (int j=0; j<nbBranches; j++) {
                if (networktable[nbRecords - 1][j] == n) {
                    muTable[nbRecords - 1][j] = muTable[0][branch0 - n] + muTable[0][branch0 + n];
                    break;
                }
            }
        }
    }
    for (int i=1; i<nbRecords; i++)      // branch zero
        muTable[0][branch0] += muTable[i][branch0];
}

void QNetwork::display_muTable(){
    display_table(muTable, "muTable");
}

void QNetwork::compute_mu_from_table(){
    if (muTable.empty()){
        cerr << "muTable not defined" << endl;
        exit(1);
    }
    else {
        int branch0 = (int) (nbBranches / 2);
        mu = muTable[0][branch0];
    }
    // loop over all elements
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_mu_from_table();
        for_each_proton(&TProton::compute_mu_from_table);
        for_each_neutron(&TNeutron::compute_mu_from_table);
        for_each_segment(&TSegment::compute_mu_from_table);
        for_each_branch(&TBranch::compute_mu_from_table);
    }
}

void QNetwork::compute_massflow_from_mu(){
    double Mdot = massflow_from_mu(mu, CentralPlant->max_massflow, Z);
    flow->massflow(Mdot);
    // loop over all elements
    if(Table_of_Segments.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_massflow_from_mu();
        for_each_proton(&TProton::compute_massflow_from_mu);
        for_each_neutron(&TNeutron::compute_massflow_from_mu);
        for_each_segment(&TSegment::compute_massflow_from_mu);
        for_each_branch(&TBranch::compute_massflow_from_mu);
    }
}



void QNetwork::compute_mu_from_massflow(){
    mu = mu_from_massflow(flow->massflow(), CentralPlant->max_massflow, Z);
    // loop over all elements
    if(Table_of_Segments.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_mu_from_massflow();
        for_each_proton(&TProton::compute_mu_from_massflow);
        for_each_neutron(&TNeutron::compute_mu_from_massflow);
        for_each_segment(&TSegment::compute_mu_from_massflow);
        for_each_branch(&TBranch::compute_mu_from_massflow);
    }
}

void QNetwork::compute_mubar_from_mu(){
    Quantum::compute_mubar_from_mu();
    // loop over all elements
    if(Table_of_Segments.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_mubar_from_mu();
        for_each_proton(&TProton::compute_mubar_from_mu);
        for_each_neutron(&TNeutron::compute_mubar_from_mu);
        for_each_segment(&TSegment::compute_mubar_from_mu);
        for_each_branch(&TBranch::get_mu_from_segment);
        for_each_branch(&TBranch::compute_mubar_from_mu);
    }
}

/*** compute demand functions ***/
void QNetwork::compute_demand(){
    compute_demand(CentralPlant);
    demand = CentralPlant->demand;
}

void QNetwork::compute_demand(Quantum* pQuantum){
    if (!pQuantum->next.empty()) {     // does nothing for TElectron
        pQuantum->demand = 0.0;
        for (int i = 0; i < pQuantum->next.size(); i++) {
            compute_demand(pQuantum->next[i]);
            pQuantum->demand += pQuantum->next[i]->demand;
        }
    }
}

/*** loop over list/table of elements and apply member function ***/
template<typename Q> void QNetwork::for_each_element(vector<Q*> List_of_Elements, void(Q::*ptrfunc)()){
    for (auto element : List_of_Elements)
    (element->*ptrfunc)();
}

template<typename Q, typename T>
void QNetwork::for_each_element(vector<Q*> List_of_Elements, void(Q::*ptrfunc)(T), T t){
    for (auto element : List_of_Elements)
        (element->*ptrfunc)(t);
}

void QNetwork::for_each_electron(void(TElectron::*ptrfunc)()){
    for_each_element(List_of_Electrons, ptrfunc);
}

template<typename T> void QNetwork::for_each_electron(void(TElectron::*ptrfunc)(T), T t){
    for_each_element(List_of_Electrons, ptrfunc, t);
}

void QNetwork::for_each_proton(void(TProton::*ptrfunc)()){
    for_each_element(List_of_Protons, ptrfunc);
}

template<typename T> void QNetwork::for_each_proton(void(TProton::*ptrfunc)(T), T t){
    for_each_element(List_of_Protons, ptrfunc, t);
}

void QNetwork::for_each_neutron(void(TNeutron::*ptrfunc)()){
    for_each_element(List_of_Neutrons, ptrfunc);
}

template<typename T> void QNetwork::for_each_neutron(void(TNeutron::*ptrfunc)(T), T t){
    for_each_element(List_of_Neutrons, ptrfunc, t);
}

void QNetwork::for_each_branch(void(TBranch::*ptrfunc)()){
    for_each_element(List_of_Branches, ptrfunc);
}

template<typename T> void QNetwork::for_each_branch(void(TBranch::*ptrfunc)(T), T t){
    for_each_element(List_of_Branches, ptrfunc, t);
}

void QNetwork::for_each_segment(void(TSegment::*ptrfunc)()){
    for (int j = 0; j < Table_of_Segments[0].size(); j++)
        for (auto & Line_of_Table : Table_of_Segments)
            if (Line_of_Table[j] == nullptr)
                break;      // goto next branch
            else
                (Line_of_Table[j]->*ptrfunc)();
}

template<typename T> void QNetwork::for_each_segment(void(TSegment::*ptrfunc)(T), T t){
    for (int j = 0; j < Table_of_Segments[0].size(); j++)
        for (auto & Line_of_Table : Table_of_Segments)
            if (Line_of_Table[j] == nullptr)
                break;      // goto next branch
            else
                (Line_of_Table[j]->*ptrfunc)(t);
}

/*** functions to define/create elements of network ***/
void QNetwork::define_electrons(vector<TElectron> &elec){
    if (elec.size() != Z){
        cerr << "Size do not match!" << endl;
        exit(1);
    }
    List_of_Electrons.clear();
    List_of_Electrons.reserve(Z);
    for (int i=0; i<Z; i++) {
        TElectron* pElectron = &elec[i];
        pElectron->QNet = this;
        if (pElectron->ID == 0)
            pElectron->ID = List_of_ID[i];
        List_of_Electrons.push_back(pElectron);
    }
}

void QNetwork::create_electrons(){
    electrons_created=true;
    List_of_Electrons.clear();
    List_of_Electrons.reserve(Z);
    for (int i = 0; i < Z; i++) {
        auto* pElectron = new TElectron (List_of_ID[i]);
        pElectron->QNet = this;
        List_of_Electrons.push_back(pElectron);
    }
}

void QNetwork::define_centralplant(TCentralPlant& central){
    CentralPlant = &central;
    CentralPlant->QNet = this;
}

void QNetwork::create_centralplant(){
    centralplant_created=true;
    CentralPlant = new TCentralPlant;
    CentralPlant->QNet = this;
}

void QNetwork::define_protons(vector<TProton> &prot){
    if(prot.size() != Z){
        cerr << "Size do not match!" << endl;
        exit(1);
    }
    List_of_Protons.clear();
    List_of_Protons.reserve(Z);
    for (int i=0; i<Z; i++) {
        TProton* pProton = &prot[i];
        pProton->QNet = this;
        if (pProton->ID == 0)
            pProton->ID = i + 1;
        List_of_Protons.push_back(pProton);
    }
}

void QNetwork::create_protons(){
    protons_created=true;
    List_of_Protons.clear();
    List_of_Protons.reserve(Z);
    for (int i = 0; i < Z; i++){
        auto* pProton = new TProton (i + 1);
        pProton->QNet = this;
        List_of_Protons.push_back(pProton);
    }
}

void QNetwork::create_protons_from_electrons(){
    protons_created=true;
    List_of_Protons.clear();
    List_of_Protons.reserve(Z);
    for (int i = 0; i < Z; i++){
        auto* pProton = new TProton (List_of_Electrons[i]->ID, List_of_Electrons[i]->l, List_of_Electrons[i]->k);
        pProton->QNet = this;
        List_of_Protons.push_back(pProton);
    }
}

void QNetwork::define_neutrons(vector<TNeutron> &neut){
    if (N==0){
        cerr << "Compute N first!" << endl;
        exit(1);
    }
    if(neut.size() != N){
        cerr << "Size do not match!" << endl;
        exit(1);
    }
    List_of_Neutrons.clear();
    List_of_Neutrons.resize(N);
    for (int i=0; i<N; i++) {
        List_of_Neutrons[neut[i].n -1] = &neut[i];
        List_of_Neutrons[neut[i].n -1]->QNet = this;
    }
}

void QNetwork::create_neutrons(){
    neutrons_created=true;
    List_of_Neutrons.clear();
    List_of_Neutrons.reserve(N);
    for (int i = 0; i < N; i++){
        auto* pNeutron = new TNeutron (i + 1);
        pNeutron->QNet = this;
        List_of_Neutrons.push_back(pNeutron);
    }
}

void QNetwork::define_branches(vector<TBranch> &branch){
    if (N==0){
        cerr << "Compute N first!" << endl;
        exit(1);
    }
    if(branch.size() != 2*N+1){
        cerr << "Size do not match!" << endl;
        exit(1);
    }
    List_of_Branches.clear();
    List_of_Branches.resize(N);
    for (int i=0; i<2*N+1; i++) {
        List_of_Branches[branch[i].l + N] = &branch[i];
        List_of_Branches[branch[i].l + N]->QNet = this;
    }
}

void QNetwork::create_branches(){
    branches_created=true;
    List_of_Branches.clear();
    List_of_Branches.reserve(N);
    for (int i = 0; i < 2*N+1; i++){
        auto* pBranch = new TBranch (i - N);
        pBranch->QNet = this;
        List_of_Branches.push_back(pBranch);
    }
}

void QNetwork::define_segments(vector<TSegment> &seg){
    if(seg.size() != A){
        cerr << "Size do not match!" << endl;
        exit(1);
    }
    if (lengthtable.empty() || diamtable.empty()){
        cerr << "Read diameters and lengths first!" << endl;
        exit(1);
    }
    Table_of_Segments.clear();
    Table_of_Segments.resize(lengthtable.size(), vector<TSegment*>(nbBranches, nullptr));
    int branch0 = (int)(nbBranches/2);
    for (int i=0; i<A; i++) {
        Table_of_Segments[seg[i].k - 1][seg[i].l + branch0] = &seg[i];
        Table_of_Segments[seg[i].k - 1][seg[i].l + branch0]->QNet = this;
    }
}

void QNetwork::create_segments(){
    if (lengthtable.empty() || diamtable.empty()){
        cerr << "Read diameters and lengths first!" << endl;
        exit(1);
    }
    segments_created=true;
    Table_of_Segments.clear();
    Table_of_Segments.resize(lengthtable.size(), vector<TSegment*>(nbBranches, nullptr));
    int branch0 = (int)(nbBranches/2);
    bool found;

    for (int n=0; n<nbBranches; n++){      // loop over table
        found = false;
        for (int m=0; m < lengthtable.size(); m++){
            if ((m < nbRecords-2) && networktable[m+1][n] > 0){    // segment followed by proton
                Table_of_Segments[m][n] = new TSegment (n-branch0, m+1);
                Table_of_Segments[m][n]->QNet = this;
            }
            else if (networktable[nbRecords-1][n] > 0 && (!found)){     // segment followed by neutron
                Table_of_Segments[m][n] = new TSegment (n-branch0, m+1);
                Table_of_Segments[m][n]->QNet = this;
                found = true;
            }
        }
    }
}

/*** function to init electrons from a list of data ***/
void QNetwork::init_electrons(double all_electrons[], string property){
    for_each_electron(&TElectron::find_kl_from_ID);

    if (property=="mu")
    {
        for (int i=0; i<List_of_Electrons.size(); i++)
            List_of_Electrons[i]->mu = all_electrons[i];
    }
    else if (property=="Mdot")
    {
        for (int i=0; i<List_of_Electrons.size(); i++) {
            List_of_Electrons[i]->flow->massflow(all_electrons[i]);
            List_of_Electrons[i]->compute_mu_from_massflow();
        }
    }
    else if (property=="delta")
    {
        for (int i=0; i<List_of_Electrons.size(); i++)
            List_of_Electrons[i]->delta = all_electrons[i];
    }
    else if (property=="dT")
    {
        for (int i=0; i<List_of_Electrons.size(); i++)
            List_of_Electrons[i]->dT = all_electrons[i];
    }
    else {
        cerr << "property not implemented" << endl;
        exit(1);
    }
}

/*** function k0 ***/
void QNetwork::compute_k0() {
    if (networktable.empty() || CentralPlant==nullptr){
        cerr << "Network needs to be defined" << endl;
        exit(1);
    }
    else
        k0 = CentralPlant->max_massflow * flow->Cp/Z;      // max_massflow
        //k0 = flow->massflow() * flow->Cp/Z;                  // time-dependent massflow
}

/*** functions ambient temperature ***/
void QNetwork::get_Toutside(environment &ext){
    Toutside = ext.T[timestep];
}

void QNetwork::get_Tground(environment &ext){
    Tground = ext.Tsol[timestep];
}

void QNetwork::get_Toutside(){
    get_Toutside(*meteo);
}

void QNetwork::get_Tground(){
    get_Tground(*meteo);
}

void QNetwork::init_Tambient() {
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cerr << "Init network first" << endl;
        exit(1);
    }
    Tambient = &Tground;        // TODO implement switch between above/below ground for each element
    for_each_electron(&TElectron::init_Tambient, true);
    for_each_proton(&TProton::init_Tambient, true);
    for_each_neutron(&TNeutron::init_Tambient, true);
    for_each_segment(&TSegment::init_Tambient, true);
    for_each_branch(&TBranch::init_Tambient, true);       //TODO does not make sense for a branch
    CentralPlant->init_Tambient(false);
}

void QNetwork::init_Tambient(const bool below) {
    cerr << "Do not use!" << endl;
    exit(1);
}

void QNetwork::init_Tref() {
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cerr << "Init network first" << endl;
        exit(1);
    }
    Tref = &T0;
    for_each_electron(&TElectron::init_Tref);
    for_each_proton(&TProton::init_Tref);
    for_each_neutron(&TNeutron::init_Tref);
    for_each_segment(&TSegment::init_Tref);
    for_each_branch(&TBranch::init_Tref);
    CentralPlant->init_Tref();
}

void QNetwork::init_Tsource() {
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cerr << "Init network first" << endl;
        exit(1);
    }
    Quantum::init_Tsource(CentralPlant);
    for (int j=0; j<Table_of_Segments[0].size(); j++)
        for (int i=0; i<Table_of_Segments.size(); i++)
            if (Table_of_Segments[i][j] == nullptr)
                break;  //go to next branch
            else
                Table_of_Segments[i][j]->init_Tsource(CentralPlant);

    for (auto & pElectron : List_of_Electrons)
        pElectron->init_Tsource(CentralPlant);
    for (auto & pProton : List_of_Protons)
        pProton->init_Tsource(CentralPlant);
    for (auto & pNeutron : List_of_Neutrons)
        pNeutron->init_Tsource(CentralPlant);
    for (auto & pBranch : List_of_Branches)
        pBranch->init_Tsource(CentralPlant);
    CentralPlant->init_Tsource(CentralPlant);
}


/*** functions connect network ***/
void QNetwork::connect_network(){
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cerr << "Init network first" << endl;
        exit(1);
    }
    for_each_segment(&TSegment::connect_network);
    for_each_electron(&TElectron::connect_network);
    for_each_proton(&TProton::connect_network);
    for_each_neutron(&TNeutron::connect_network);
    for_each_branch(&TBranch::connect_network);
    CentralPlant->connect_network();
}

void QNetwork::connect_flow(){
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cerr << "Init network first" << endl;
        exit(1);
    }
    for_each_segment(&TSegment::connect_flow);
    for_each_electron(&TElectron::connect_flow);
    for_each_proton(&TProton::connect_flow);
    for_each_neutron(&TNeutron::connect_flow);
    for_each_branch(&TBranch::connect_flow);
    CentralPlant->connect_flow();
}

/*** functions delta ***/
void QNetwork::compute_deltabar(){
    delta_bar=0.0;
    for (int i=0; i<Z; i++)
        delta_bar += List_of_Electrons[i]->delta;
    delta_bar /= Z;
}

/*** function pressure loss ***/
void QNetwork::compute_dP(){
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
    List_of_Branches.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    dP = 0.0;
    for_each_segment(&TSegment::compute_dP);
    for_each_electron(&TElectron::compute_dP);
    for_each_proton(&TProton::compute_dP);
    for_each_neutron(&TNeutron::compute_dP);
    for_each_branch(&TBranch::compute_dP);

    //find most unfavourable way with highest dP
    vector<double> dP_maxbelow (N+1, 0);  //stores highest dP below neutron n, with dP_maxbelow[0]=0
    int branch0 = (int)(nbBranches/2);
    for (int n = N; n>0; n--){
            dP_maxbelow[n] = max(List_of_Branches[-n+branch0]->dP + dP_maxbelow[networktable[nbRecords-1][-n+branch0]],
                                 List_of_Branches[n+branch0]->dP + dP_maxbelow[networktable[nbRecords-1][n+branch0]]);
    }
    dP = List_of_Branches[branch0]->dP + dP_maxbelow[networktable[nbRecords-1][branch0]];
    CentralPlant->dP = dP;
}


/*** function tau ***/
void QNetwork::compute_tau(){
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    tau = 0.0;
    for_each_segment(&TSegment::compute_tau);
    for_each_electron(&TElectron::compute_tau);
    for_each_proton(&TProton::compute_tau);
    for_each_neutron(&TNeutron::compute_tau);
    for_each_branch(&TBranch::compute_tau);
}

/*** function COP ***/
void QNetwork::compute_COP(){
    if(List_of_Electrons.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    for (auto &pElectron: List_of_Electrons)
        if (auto *pElectronHP = dynamic_cast<TElectronHP *>(pElectron))  // only valid for Heatpumps
            pElectronHP->compute_COP();
}

/*** function approximative power loss ***/
void QNetwork::estimate_powerloss() {
    if(Table_of_Segments.empty() || List_of_Electrons.empty() || List_of_Protons.empty() ||
       List_of_Branches.empty()){
        cerr << "Init network first" << endl;
        exit(1);
    }
    dPower = 0.0;
    for (int j=0; j<Table_of_Segments[0].size(); j++)
        for (int i=0; i<Table_of_Segments.size(); i++)
            if (Table_of_Segments[i][j] == nullptr)
                break;
            else {
                Table_of_Segments[i][j]->estimate_powerloss();
                dPower += Table_of_Segments[i][j]->dPower;
            }
    for (auto & pElectron : List_of_Electrons) {
        if (pElectron->DN !=0) {
            pElectron->estimate_powerloss();
            dPower += pElectron->dPower;
        }
    }
    for (auto & pProton : List_of_Protons) {
        pProton->estimate_powerloss();
        dPower += pProton->dPower;
    }
    for (auto & pNeutron : List_of_Neutrons) {
        pNeutron->estimate_powerloss();
        dPower += pNeutron->dPower;
    }
    for (auto & pBranch : List_of_Branches) {
        pBranch->estimate_powerloss();
    }
    CentralPlant->dPower = dPower;
}

void QNetwork::compute_nubar(){
    Quantum::compute_nubar();
    // loop over all elements
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_nubar();
        for_each_electron(&TElectron::compute_nubar);
        for_each_proton(&TProton::compute_nubar);
        for_each_neutron(&TNeutron::compute_nubar);
        for_each_segment(&TSegment::compute_nubar);
        for_each_branch(&TBranch::compute_nubar);
    }
}

void QNetwork::compute_epsilon(){
    Quantum::compute_epsilon();
    // loop over all elements
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_epsilon();
        for_each_electron(&TElectron::compute_epsilon);
        for_each_proton(&TProton::compute_epsilon);
        for_each_neutron(&TNeutron::compute_epsilon);
        for_each_segment(&TSegment::compute_epsilon);
        for_each_branch(&TBranch::compute_epsilon);
    }
}

void QNetwork::compute_power_from_massflow(){
    Quantum::compute_power_from_massflow();
    // loop over all elements except TElectron
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_power_from_massflow();
        for_each_proton(&TProton::compute_power_from_massflow);
        for_each_neutron(&TNeutron::compute_power_from_massflow);
        for_each_segment(&TSegment::compute_power_from_massflow);
        for_each_branch(&TBranch::compute_power_from_massflow);
    }
}

void QNetwork::compute_material_cost() {
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        material_cost = 0.0;
        for (int j=0; j<Table_of_Segments[0].size(); j++)
            for (int i=0; i<Table_of_Segments.size(); i++)
                if (Table_of_Segments[i][j] == nullptr)
                    break;
                else {
                    Table_of_Segments[i][j]->compute_material_cost();
                    material_cost += Table_of_Segments[i][j]->material_cost;
                }
        for (auto & pElectron : List_of_Electrons) {
            if (pElectron->DN !=0) {
                pElectron->compute_material_cost();
                material_cost += pElectron->material_cost;
            }
        }
        for (auto & pProton : List_of_Protons) {
            pProton->compute_material_cost();
            material_cost += pProton->material_cost;
        }
        for (auto & pNeutron : List_of_Neutrons) {
            pNeutron->compute_material_cost();
            material_cost += pNeutron->material_cost;
        }
        for (auto & pBranch : List_of_Branches) {
            pBranch->compute_material_cost();
        }
        CentralPlant->compute_material_cost();
        material_cost += CentralPlant->material_cost;
    }
}

void QNetwork::compute_engineering_cost(const char* region) {
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        engineering_cost = 0.0;
        for (int j=0; j<Table_of_Segments[0].size(); j++)
            for (int i=0; i<Table_of_Segments.size(); i++)
                if (Table_of_Segments[i][j] == nullptr)
                    break;
                else {
                    Table_of_Segments[i][j]->compute_engineering_cost(region);
                    engineering_cost += Table_of_Segments[i][j]->engineering_cost;
                }
        for (auto & pElectron : List_of_Electrons) {
            if (pElectron->DN !=0) {
                pElectron->compute_engineering_cost(region);
                engineering_cost += pElectron->engineering_cost;
            }
        }
        for (auto & pProton : List_of_Protons) {
            pProton->compute_engineering_cost(region);
            engineering_cost += pProton->engineering_cost;
        }
        for (auto & pNeutron : List_of_Neutrons) {
            pNeutron->compute_engineering_cost(region);
            engineering_cost += pNeutron->engineering_cost;
        }
        for (auto & pBranch : List_of_Branches) {
            pBranch->compute_engineering_cost(region);
        }
        CentralPlant->compute_engineering_cost(region);
        engineering_cost += CentralPlant->engineering_cost;
    }
}

void QNetwork::compute_total_invest() {
    Quantum::compute_total_invest();
    // loop over all elements
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        CentralPlant->compute_total_invest();
        for_each_electron(&TElectron::compute_total_invest);
        for_each_proton(&TProton::compute_total_invest);
        for_each_neutron(&TNeutron::compute_total_invest);
        for_each_segment(&TSegment::compute_total_invest);
        for_each_branch(&TBranch::compute_total_invest);
    }
}

void QNetwork::compute_operating_cost() {
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        operating_cost = 0.0;
        for (int j=0; j<Table_of_Segments[0].size(); j++)
            for (int i=0; i<Table_of_Segments.size(); i++)
                if (Table_of_Segments[i][j] == nullptr)
                    break;
                else {
                    Table_of_Segments[i][j]->compute_operating_cost();
                    operating_cost += Table_of_Segments[i][j]->operating_cost;
                }
        for (auto & pElectron : List_of_Electrons) {
            if (pElectron->DN !=0) {
                pElectron->compute_operating_cost();
                operating_cost += pElectron->operating_cost;
            }
        }
        for (auto & pProton : List_of_Protons) {
            pProton->compute_operating_cost();
            operating_cost += pProton->operating_cost;
        }
        for (auto & pNeutron : List_of_Neutrons) {
            pNeutron->compute_operating_cost();
            operating_cost += pNeutron->operating_cost;
        }
        for (auto & pBranch : List_of_Branches) {
            pBranch->compute_operating_cost();
        }
        CentralPlant->compute_operating_cost();
        operating_cost += CentralPlant->operating_cost;
    }
}

void QNetwork::compute_resources_cost() {
    if(Table_of_Segments.empty() || List_of_Electrons.empty()|| List_of_Protons.empty() ||
       List_of_Branches.empty() || CentralPlant==nullptr){
        cout << "Warning: Network not initialized" << endl;
    }
    else {
        resources_cost = 0.0;   // resources cost in this timestep
        if (timestep == 0)
            total_resources_cost = 0.0;     // cumulated resources cost

        for (int j=0; j<Table_of_Segments[0].size(); j++)
            for (int i=0; i<Table_of_Segments.size(); i++)
                if (Table_of_Segments[i][j] == nullptr)
                    break;
                else {
                    Table_of_Segments[i][j]->compute_resources_cost();
                }
        for (auto & pElectron : List_of_Electrons) {
            if (pElectron->DN !=0) {
                pElectron->compute_resources_cost();
            }
        }
        for (auto & pProton : List_of_Protons) {
            pProton->compute_resources_cost();
        }
        for (auto & pNeutron : List_of_Neutrons) {
            pNeutron->compute_resources_cost();
        }
        for (auto & pBranch : List_of_Branches) {
            pBranch->compute_resources_cost();
        }
        CentralPlant->compute_resources_cost();
        resources_cost = CentralPlant->resources_cost;
        total_resources_cost += resources_cost;
    }
}



/*** solvers ***/
void QNetwork::init() {
    if (settings.empty()){
        cerr << "Read settings file first" << endl;
        exit(1);
    }

    // set filepath variable
    set_filepath();

    // get files from database
    if (use_database) {
        connect_to_database();
        get_geopackages_from_database();    // download of Network and Data
    }
}

void QNetwork::extract_quantum_network() {
    if (filepath.empty() || settings.empty() || (use_database && (connection==nullptr || GPKGtable==nullptr)))
        init();

    if (settings.find("Network") != settings.end() && settings.find("Q-Network") != settings.end()) {   // both keys found
        auto settings_prep = parse_setting("Network");
        auto settings_qnet = parse_setting("Q-Network");

        string input_network = filepath + settings_prep[0] + ".gpkg";
        string output_network = filepath + settings_qnet[0] + ".gpkg";

        string sstfile = filepath + "test2_SST.csv";  // TODO version gpkg
        const double tol = 0.1;
        const char *distrib_filter = "Genre_cond = 'Retour' AND Fonction = 'Conduite de distribution'";
        const char *raccord_filter = "Genre_cond = 'Retour' AND Fonction = 'Conduite de raccordement'";

        // function from "GeoDataset.h"
        ::extract_quantum_network(input_network.c_str(), output_network.c_str(), sstfile.c_str(),
                                  settings_prep[1].c_str(), settings_prep[2].c_str(),
                                  settings_qnet[1].c_str(), settings_qnet[2].c_str(), settings_qnet[3].c_str(),
                                  tol, distrib_filter, raccord_filter);
    }


}

void QNetwork::init_network(){
    if (filepath.empty() || settings.empty() || (use_database && (connection==nullptr || GPKGtable==nullptr)))
        init();

    // read files
    read_meteo();
    read_tables();      // networktable, diamtable, lengthtable

    // compute some quantum numbers
    compute_ANZ();
    compute_mtable();
    compute_m_from_table();

    read_pipetable();  // pipetable

    create_centralplant();
    CentralPlant->compute_m_from_table();
    try {
        CentralPlant->max_massflow = stod(get_setting("Maximum Mass Flow"));
    } catch(...) {}    // just ignore exception
    try {
        CentralPlant->max_power = stod(get_setting("Maximum Power"));
    } catch(...) {}
    CentralPlant->init_pressure(500000, 500000);      //TODO
    CentralPlant->init_temperature(stod(get_setting("Supply Temperature")),
                                   stod(get_setting("Estimated Return Temperature")));    // estimation, Celsius or Kelvin
    CentralPlant->init_fluid(get_setting("Fluid"));
    CentralPlant->init_boilers();
    CentralPlant->compute_maxvalues_from_boilers();

    create_electrons_from_SST();
    for (auto & pElectron : List_of_Electrons){
        pElectron->compute_m_from_table();
        pElectron->find_kl_from_ID();
        pElectron->init_DN();
        pElectron->init_pipeprop();
        if (pElectron->D.power.empty())     // not yet read from file, use SIA norm
            pElectron->D.calcul_power(*meteo);            // sum of all buildings, only one building here
        else if (pElectron->D.power.size() != meteo->Time){
            cerr << "Lengths of meteo file and SST demand file are not the same" << endl;
            exit(1);
        }
        if ((!pElectron->dTvar.empty()) && (pElectron->dTvar.size() != meteo->Time)){
            cerr << "Lengths of meteo file and dT file are not the same" << endl;
            exit(1);
        }
        pElectron->D.B.clear();                              // free memory
    }
    create_protons_from_electrons();
    for_each_proton(&TProton::compute_m_from_table);

    create_neutrons();
    for_each_neutron(&TNeutron::compute_m_from_table);

    create_segments(); // init Pipe
    for (int j=0; j<Table_of_Segments[0].size(); j++){
        for (int i=0; i<Table_of_Segments.size(); i++){
            TSegment* pSegment=Table_of_Segments[i][j];
            if (pSegment == nullptr)
                break;
            else {
                pSegment->compute_m_from_table();
                pSegment->init_length();
                pSegment->init_DN();
                pSegment->init_pipeprop();
            }
        }
    }

    create_branches();
    for_each_branch(&TBranch::compute_m_from_table);
    for_each_branch(&TBranch::get_segments);
    for_each_branch(&TBranch::get_protons);


    // free memory
    List_of_ID.clear();
    diamtable.clear();
    lengthtable.clear();
    pipetable.clear();

    // init T-pointers
    init_Tambient();
    init_Tref();
    init_Tsource();

    // init fluid properties of all elements
    flow->copy(CentralPlant->flow);
    dT = CentralPlant->dT;
    for (auto & pElectron : List_of_Electrons) {
        pElectron->flow->copy(flow);
        pElectron->compute_delta_from_dT();
    }
    compute_deltabar();
    for (auto & pProton : List_of_Protons) {
        pProton->flow->copy(flow);
        pProton->dT = dT;
    }
    for (auto & pNeutron : List_of_Neutrons) {
        pNeutron->flow->copy(flow);
        pNeutron->dT = dT;
    }
    for (int j=0; j<Table_of_Segments[0].size(); j++)
        for (int i=0; i<Table_of_Segments.size(); i++) {
            TSegment *pSegment = Table_of_Segments[i][j];
            if (pSegment == nullptr)
                break;
            else {
                pSegment->flow->copy(flow);
                pSegment->dT = dT;
            }
        }
    for (auto & pBranch : List_of_Branches) {
        pBranch->flow->copy(flow);
        pBranch->dT = dT;
    }
    compute_k0();
    connect_network();
    //connect_flow();

    if (stringtobool(get_setting("Solve Cost").c_str())) {
        compute_material_cost();
        compute_engineering_cost();
        compute_total_invest();
        compute_operating_cost();
    }

    //solver settings
    max_timestep=meteo->Time;

    cout << "init() completed\n";
}

void QNetwork::init_timestep(){
    // init time-dependent "physics" of network
    get_Toutside();
    get_Tground();
    //T0 = Toutside;      //TODO does this make sense?
    compute_COP();
    for_each_electron(&TElectron::update_dT);
    for_each_electron(&TElectron::compute_power_from_district);
    for_each_electron(&TElectron::compute_mu_from_power);
    for_each_electron(&TElectron::compute_mubar_from_mu);
}

void QNetwork::solve_timestep(){
    solve_timestep(stringtobool(get_setting("Solve Mass Flow").c_str()),
                   stringtobool(get_setting("Solve Temperature").c_str()),
                   get_setting("Solver")=="Full",
                   stringtobool(get_setting("Update COP").c_str()),
                   stringtobool(get_setting("Solve Cost").c_str()));
}

void QNetwork::solve_timestep(const bool MF, const bool T, const bool Tsolver, const bool COP, const bool C){
    cout << "Timestep " << timestep;
    if (MF){    // solve mass flow
        solve_massflow();
        compute_dP();
        compute_demand();
    }
    if (T) {    // solve temperature
        compute_tau();

        if (Tsolver)  // detailed solver
            solve_temperature();
        else {        // approximate solver
            estimate_powerloss();
            compute_nubar();
        }
        compute_power_from_massflow();
        compute_epsilon();
        if (COP)            // update COP for next timestep
            compute_COP();
    }
    if (C) {
        compute_resources_cost();
    }
    cout << "\n";
}

void QNetwork::solve_massflow(){
//    compute_muTable();
//    compute_mu_from_table();
    compute_mu();
    compute_mubar_from_mu();
    compute_massflow_from_mu();
}

void QNetwork::solve_temperature(){
    if (CentralPlant== nullptr){
        cerr << "Init network first" << endl;
        exit(1);
    }
    dPower = 0.0;
    solve_temperature(CentralPlant);

    // dPower, dT, delta for all Branches
    for_each_branch(&TBranch::estimate_powerloss);
    for_each_branch(&TBranch::get_TIn_TOut_dT);
    for_each_branch(&TBranch::compute_delta_from_dT);

    // TOut, dT, delta for *this
    flow->pStateOut->T = CentralPlant->flow->TOut();
    dT = CentralPlant->dT;
    compute_delta_from_dT();
}

void QNetwork::solve_temperature(Quantum* pQuantum){
    pQuantum->compute_T_supply();
    for (auto & next : pQuantum->next) {
        solve_temperature(next);
    }
    pQuantum->compute_T_return();
    dPower += pQuantum->dPower;
    pQuantum->compute_nubar_from_T(k0);
}

void QNetwork::next_timestep(){
    timestep++;
}

/*** save result files ***/
void QNetwork::save_results_to_file(){
    if(stringtobool(get_setting("Save Results").c_str())) {
//        static vector<bool[7]> quantities (properties.size());
        static vector<vector<bool>> quantities;
        quantities.resize(properties.size(), vector<bool>(get_setting("massflow").length()));
        static unsigned int frequency=1e8;
        // initialise in first time step
        if (timestep == 0) {
            for (int i = 0; i < quantities.size(); i++) {
                string setting;
                setting = get_setting(properties[i]);
                for (int j = 0; j < quantities[0].size(); j++)
                    quantities[i][j] = stringtobool(setting[j]);
            }
            try {
                frequency = stoul(get_setting("Save Frequency"));
            } catch (const std::invalid_argument &) {}    // just ignore exception
        }
        save_results_to_file(quantities, frequency);
    }
}


void QNetwork::save_results_to_file(const vector<vector<bool>>& quantities, unsigned int frequency){
    static vector<vector<Quantum*>> printElements (properties.size());   //elements to be printed
    static vector<ofstream> results (properties.size());
    static vector<int> prop;

    if (timestep==0) {
        // find properties in question
        for (int i = 0; i < properties.size(); i++)
            for (int j = 0; j < quantities[0].size(); j++)
                if (quantities[i][j]) {
                    prop.push_back(i);
                    break;
                }

        //char buf[256 * 1024];
        string file, outfile;
        file = get_setting("Case Directory");
        file += "/";
        file += get_setting("Case Name");
        file += "_results_";
        for (auto &p: prop) {
            //TODO not working: results[p].rdbuf()->pubsetbuf(buf, sizeof(buf));

            // create and open all necessary files
            outfile = file;
            outfile += properties[p];
            outfile += ".csv";
            results[p].open(outfile);

            // first timestep -> write title line
            printElements[p] = print_headline(results[p], quantities[p]);
        }
    }

    else if (timestep%frequency == 0)          // all periodical timesteps
        for (auto &p: prop)
            results[p].flush();

    // all timesteps
    for (auto &p: prop)
        if (p < 15 || (timestep == max_timestep-1))     // print costs only for last time-step
            print_timestep(results[p], printElements[p], p);

}

vector<Quantum*> QNetwork::print_headline(ostream &results, const vector<bool>& elements) {
    vector<Quantum*> printElements;
    printElements.clear();

    results << "Timestep";

    if (elements[0]) {
        results << "\tTotal";
        printElements.push_back(this);
    }
    if (elements[1]) {
        results << "\tCentral Plant";
        printElements.push_back(CentralPlant);
    }
    if(elements[2])
        for (auto &pElectron: List_of_Electrons) {
            results << "\tSubstation: " << pElectron->ID;
            printElements.push_back(pElectron);
        }
    if(elements[3])
        for (auto &pProton: List_of_Protons) {
            results << "\tTee: " << pProton->ID;
            printElements.push_back(pProton);
        }
    if(elements[4])
        for (auto &pNeutron: List_of_Neutrons) {
            results << "\tBifurcation: " << pNeutron->n;
            printElements.push_back(pNeutron);
        }
    if(elements[5])
        for (int j=0; j<Table_of_Segments[0].size(); j++)
            for (int i=0; i<Table_of_Segments.size(); i++) {
                TSegment *pSegment = Table_of_Segments[i][j];
                if (pSegment == nullptr)
                    break;
                else {
                    results << "\tSegment: " << pSegment->l;
                    results << "," << pSegment->k;
                    printElements.push_back(pSegment);
                }
            }
    if(elements[6])
        for (auto &pBranch: List_of_Branches) {
            results << "\tBranch: " << pBranch->l;
            printElements.push_back(pBranch);
        }

    results << "\n";

    return printElements;
}

void QNetwork::print_timestep(ostream &results, vector<Quantum*>& printElements, int property) const {
    results << timestep + 1;
    for (auto & pQuantum : printElements)
        results << "\t"<< get_property_by_idx(pQuantum, property);
    results << "\n";
}

double QNetwork::get_property_by_idx(Quantum*& q, int& prop_id) const {
    switch (prop_id) {
        case 0:     return q->flow->massflow();
        case 1:     return q->speed;
        case 2:     return q->mu;
        case 3:     return q->mu_bar;
        case 4:     return q->dP;
        case 5:     return q->dP_m;
        case 6:     return q->demand;
        case 7:     return q->flow->Power;
        case 8:     return q->dPower;
        case 9:     return q->flow->TIn();
        case 10:    return q->flow->TOut();
        case 11:    return q->dT;
        case 12:    return q->delta;
        case 13:    return q->nu_bar;
        case 14:    return q->epsilon;
        case 15:    return q->material_cost;
        case 16:    return q->engineering_cost;
        case 17:    return q->operating_cost;
        case 18:    return q->total_resources_cost;
        default:
            cerr << "Wrong property index" << endl;
            exit(1);
    }
}

void QNetwork::save_results_to_gpkg(){
    string gpkgfile;
    gpkgfile = filepath;
    gpkgfile += get_setting("Results");
    gpkgfile += ".gpkg";

    if(stringtobool(get_setting("Save Results").c_str())) {
//        static vector<bool[7]> quantities (properties.size());
        static vector<vector<bool>> quantities;
        quantities.resize(properties.size(), vector<bool>(get_setting("massflow").length()));
        static unsigned int frequency=1e8;
        // initialise in first time step
        if (timestep == 0) {
            for (int i = 0; i < quantities.size(); i++) {
                string setting;
                setting = get_setting(properties[i]);
                for (int j = 0; j < quantities[0].size(); j++)
                    quantities[i][j] = stringtobool(setting[j]);
            }
            try {
                frequency = stoul(get_setting("Save Frequency"));
            } catch (const std::invalid_argument &) {}    // just ignore exception
        }
        save_results_to_gpkg(gpkgfile.c_str(), quantities, frequency);
    }
}

void QNetwork::save_results_to_gpkg(const char* gpkgfile, const vector<vector<bool>>& quantities, unsigned int frequency){
    static vector<OGRLayer*> poLayers;
    static vector<vector<Quantum*>> printElements (properties.size());   //elements to be printed
    static vector<vector<int>> idxElements (properties.size());
    static vector<int> prop;    // indices of properties to be printed

    static GeoDataset outdata (gpkgfile, GeoDataset::Accessmode::WRITE);

    if (timestep==0) {
        // start transaction (if supported)
        outdata.stage_dataset();

        // find properties in question
        for (int i = 0; i < properties.size(); i++)
            for (int j = 0; j < quantities[0].size(); j++)
                if (quantities[i][j]) {
                    prop.push_back(i);
                    break;
                }

        for (auto &p: prop) {
            // first timestep -> get elements
            ostringstream headline;
            string word;
            vector<string> elementNames;
            printElements[p] = print_headline(headline, quantities[p]);
            stringstream ss(headline.rdbuf()->str());
            bool skip=true;
            while (getline(ss, word, '\t')) {
                if (skip)   // skip the first entry "Timestep"
                    skip = false;
                else
                    elementNames.push_back(word);
            }
            OGRLayer* poLayer = outdata.create_Layer(properties[p], elementNames, GeoDataset::Type::DOUBLE, "Timestep");
            //outdata.get_Layer(properties[p])->SetMetadataItem("OGR_GPKG_SYNCHRONOUS", "OFF");
            poLayers.push_back(poLayer);
            for (auto & element : elementNames)
                idxElements[p].push_back(poLayer->GetLayerDefn()->GetFieldIndex(element.c_str()));
        }
    }
    else if (timestep%frequency == 0) {         // all periodical timesteps: start new transaction
        outdata.commit_dataset();
        outdata.stage_dataset();
    }

    // all timesteps
    for (auto &p: prop) {
        if (p < 15 || (timestep == max_timestep - 1)) {    // print costs only for last time-step
            // Create a feature
            OGRLayer *&poLayer = poLayers[p];
            OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
            poFeature->SetFID(timestep + 1);

            for (int q = 0; q < printElements[p].size(); q++)
                poFeature->SetField(idxElements[p][q], get_property_by_idx(printElements[p][q], p));

            if (poLayer->CreateFeature(poFeature) != OGRERR_NONE) {
                cerr << "Failed to create feature." << endl;
                exit(1);
            }
            OGRFeature::DestroyFeature(poFeature);
        }
    }

    // last timestep: commit output to gpkg
    if (timestep == max_timestep - 1) {
        outdata.commit_dataset();
        outdata.vacuum_dataset();
    }
}

void QNetwork::save_results() {
    if (use_gpkg)
        save_results_to_gpkg();
    else
        save_results_to_file();
}

void QNetwork::move_geopackage_to_database(const char* keyword, bool delete_local) const {  // use for Q-Network and Results
    string file;
    stringstream ss(get_setting(keyword));
    getline(ss, file, '[');
    file += ".gpkg";

    if (! GPKGtable->upload_gpkg((filepath + file).c_str(), file.c_str(), delete_local)) {
        cerr << "Geopackage " << file << " could not be uploaded to database" << endl;
    }
    else if (filesystem::is_empty(filepath))
        filesystem::remove(filepath);
}

void QNetwork::clear_tmp(const char* keyword) const {       // use for Preprocessed Network and Input Data (unmodified files)
    string file;
    stringstream ss(get_setting(keyword));
    getline(ss, file, '[');
    std::remove((filepath + file).c_str());

    if (filesystem::is_empty(filepath))
        filesystem::remove(filepath);
}

void QNetwork::finalize(const bool delete_local) const {
    if (settings.find("Q-Network") != settings.end() && settings.find("Results") != settings.end()) {   // both keys found
        move_geopackage_to_database("Q-Network", delete_local);
        move_geopackage_to_database("Results", delete_local);
        if (delete_local) {
            clear_tmp("Preprocessed Network");
            clear_tmp("Input Data");
        }
    }
}