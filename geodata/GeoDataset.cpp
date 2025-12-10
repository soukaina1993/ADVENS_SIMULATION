//
// Created by jamil.maj on 19.09.2024.
//

#include "GeoDataset.h"

#include <cstdlib>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;


/*** ------------------------- Base Class GeoDataset ------------------------------ ***/

// Constructor: opens or creates file
GeoDataset::GeoDataset(const char* ifilename, Accessmode rw, bool verbose){
    filename = ifilename;

    // Initialize GDAL
    GDALAllRegister();

    // Options
    char **papszOptions = nullptr;
    unsigned int nOpenFlags=GDAL_OF_VECTOR;

    if (rw == READ)
        nOpenFlags |= GDAL_OF_READONLY;
    else
        nOpenFlags |= GDAL_OF_UPDATE;

    // Try to open file (if it exists)
    poDS = GDALDataset::FromHandle(GDALOpenEx(ifilename, nOpenFlags,
                                              nullptr, papszOptions, nullptr));

    if (poDS != nullptr) {  // file exists and can be accessed
        // Get the driver used to open the dataset
        if ((poDriver=poDS->GetDriver())==nullptr) {
            cerr << "Driver cannot be determined" << endl;
        }
        else if (verbose) {
            cout << "Network File " << filename << " opened successfully.\n";

            if (poDriver != nullptr)
                cout << "Driver used: " << poDriver->GetDescription() << " - " << poDriver->
                        GetMetadataItem(GDAL_DMD_LONGNAME) << endl;

            cout << "Available layers in the dataset:\n";  // List layers in the dataset
            vector<string> layers = get_Layers();
            for (auto &layer : layers)
                cout << "'" << layer << "'" << "\n";
        }
    }
    else if (rw == WRITE){    // file does not yet exist --> create file if WRITE mode
        // Define default driver
        const char *pszFormat = "GPKG";
        // get driver from filename suffix
        if (filename.rfind(".shp") == (filename.size() - 4))
            pszFormat = "ESRI Shapefile";
        else if (filename.rfind(".gpkg") != (filename.size() - 5)) {
            cerr << "File Type not supported. GPKG used." << endl;
        }
        poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
        // driver info
        if (verbose) {
            if (poDriver == nullptr) {
                cerr << "Driver error" << endl;
                exit(1);
            } else {
                char **papszMetadata = poDriver->GetMetadata();
                if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATE, FALSE))
                    printf("Driver %s supports Create() method.\n", pszFormat);
                if (CSLFetchBoolean(papszMetadata, GDAL_DCAP_CREATECOPY, FALSE))
                    printf("Driver %s supports CreateCopy() method.\n", pszFormat);
            }
        }
        poDS = poDriver->Create(ifilename, 0, 0, 0, GDT_Unknown, papszOptions);
    }

    // result of operation
    if (poDS == nullptr) {
        cerr << "Output file " << filename << " cannot be accessed." << endl;
        exit( 1 );
    }
    else if (verbose)
        cout << "Successfully opened file " << filename << endl;
}

// Destructor: if a GDAL dataset is open, close it
GeoDataset::~GeoDataset() {
    if (poDS != nullptr) {
        GDALClose(poDS);
        cout << "Network File " << filename << " closed successfully!" << endl;
    }
}

OGRLayer* GeoDataset::create_Layer(const char* name,
                                   const vector<const char*>& fieldNames, const vector<Type>& fieldTypes, const char* fid,
                                   const GeoCoord::CRS crs) const {
    if (!(poDS->TestCapability(ODsCCreateLayer))) {
        cerr << "Layer creation for " << filename << " is not supported" << endl;
        exit(1);
    }

    if (crs != GeoCoord::CRS::WGS84)
        cerr << "Not yet implemented. WGS84 is used." << endl;

    OGRSpatialReference oSRS;
    //Swiss CRS = EPSG:2056 not supported
    if (oSRS.SetWellKnownGeogCS( "CRS84" ) == OGRERR_FAILURE)
        cerr << "SRS failed" << endl;

    char **papszOptions = nullptr;
    papszOptions = CSLSetNameValue(papszOptions, "OVERWRITE", "YES");
    papszOptions = CSLSetNameValue(papszOptions, "FID", fid);
    OGRLayer* poLayer = poDS->CreateLayer(name, &oSRS, wkbUnknown, papszOptions);

    if (poLayer == nullptr) {
        cerr << "Creation of " << name << " Layer failed." << endl;
        exit( 1 );
    }
    else if (poLayer->GetSpatialRef() == nullptr)
        cerr << "No SRS defined" << endl;

    CSLDestroy( papszOptions);

    // Create a field
    for (int i = 0; i < fieldNames.size(); i++) {
        OGRFieldType oFieldType = OFTString;
        if (fieldTypes[i] == INT)
            oFieldType = OFTInteger;
        else if (fieldTypes[i] == DOUBLE)
            oFieldType = OFTReal;
        else if (fieldTypes[i] != STR) {
            cerr << "Bad Field Type" << endl;
            exit(1);
        }

        OGRFieldDefn oField( fieldNames[i], oFieldType);
        oField.SetWidth(32);
        if (poLayer->CreateField( &oField ) != OGRERR_NONE) {
            printf("Creating %s field failed.\n", fieldNames[i]);
            exit(1);
        }
    }

//    // C code because bug in C++ OGRFieldDefn in version gdal-3.9.3
//    OGRFieldDefnH hFieldDefn;
//    for (int i = 0; i < fieldNames.size(); i++) {
//        if (fieldTypes[i] == INT)
//            hFieldDefn = OGR_Fld_Create(fieldNames[i], OFTInteger);
//        else if (fieldTypes[i] == DOUBLE)
//            hFieldDefn = OGR_Fld_Create(fieldNames[i], OFTReal);
//        else if (fieldTypes[i] == STR)
//            hFieldDefn = OGR_Fld_Create(fieldNames[i], OFTString);
//        else {
//            cerr << "Bad Field Type" << endl;
//            exit(1);
//        }
//        OGR_Fld_SetWidth(hFieldDefn, 32);
//
//        if (OGR_L_CreateField(poLayer, hFieldDefn, TRUE) != OGRERR_NONE) {
//            printf("Creating %s field failed.\n", fieldNames[i]);
//            exit(1);
//        }
//        OGR_Fld_Destroy(hFieldDefn);
//    }
    return poLayer;
}

OGRLayer* GeoDataset::create_Layer(const char* name,
                                   const vector<const char*>& fieldNames, const Type fieldType, const char* fid,
                                   const GeoCoord::CRS crs) const {
    if (!(poDS->TestCapability(ODsCCreateLayer))) {
        cerr << "Layer creation for " << filename << " is not supported" << endl;
        exit(1);
    }

    if (crs != GeoCoord::CRS::WGS84)
        cerr << "Not yet implemented. WGS84 is used." << endl;

    OGRSpatialReference oSRS;
    //Swiss CRS = EPSG:2056 not supported
    if (oSRS.SetWellKnownGeogCS( "CRS84" ) == OGRERR_FAILURE)
        cerr << "SRS failed" << endl;

    char **papszOptions = nullptr;
    papszOptions = CSLSetNameValue(papszOptions, "OVERWRITE", "YES");
    papszOptions = CSLSetNameValue(papszOptions, "FID", fid);

    OGRLayer* poLayer = poDS->CreateLayer(name, &oSRS, wkbUnknown, papszOptions);
    if (poLayer == nullptr) {
        cerr << "Creation of " << name << " Layer failed." << endl;
        exit( 1 );
    }
    else if (poLayer->GetSpatialRef() == nullptr)
        cerr << "No SRS defined" << endl;

    CSLDestroy( papszOptions);

    // Create a field
    OGRFieldType oFieldType = OFTString;
    if (fieldType == INT)
        oFieldType = OFTInteger;
    else if (fieldType == DOUBLE)
        oFieldType = OFTReal;
    else if (fieldType != STR) {
        cerr << "Bad Field Type" << endl;
        exit(1);
    }
    for (auto fieldName : fieldNames) {
        OGRFieldDefn oField( fieldName, oFieldType);
        oField.SetWidth(32);
        if (poLayer->CreateField( &oField ) != OGRERR_NONE) {
            printf("Creating %s field failed.\n", fieldName);
            exit(1);
        }
    }
    return poLayer;
}

OGRLayer* GeoDataset::create_Layer(const char* name,
                       const vector<string>& fieldNames, const vector<Type>& fieldTypes, const char* fid,
                       GeoCoord::CRS crs) const {
    vector<const char*> names;
    names.reserve(fieldNames.size());
    for (auto &fieldName : fieldNames)
        names.push_back(fieldName.c_str());
    return create_Layer(name, names, fieldTypes, fid, crs);
}

OGRLayer* GeoDataset::create_Layer(const char* name,
                       const vector<string>& fieldNames, Type fieldType, const char* fid,
                       GeoCoord::CRS crs) const {
    vector<const char*> names;
    names.reserve(fieldNames.size());
    for (auto &fieldName : fieldNames)
        names.push_back(fieldName.c_str());
    return create_Layer(name, names, fieldType, fid, crs);
}

vector<string> GeoDataset::get_Layers() {
    vector<string> layers;
    string layer;
    for (int i = 0; i < poDS->GetLayerCount(); i++) {
        layer = poDS->GetLayer(i)->GetName();
        layers.push_back(layer);
    }
    return layers;
}

OGRLayer* GeoDataset::get_Layer(const char* layername){
    return poDS->GetLayerByName(layername);
}

void GeoDataset::remove_Layer(const char* to_remove){
    for (int i = 0; i < poDS->GetLayerCount(); i++) {
        if ( strcmp(poDS->GetLayer(i)->GetName(), to_remove) == 0 ){
            poDS->DeleteLayer(i);
            break;
        }
    }
}

vector<string> GeoDataset::get_Fields_on_Layer(OGRLayer* poLayer){
    vector<string> fields;
    OGRFeatureDefn* ldefn = poLayer->GetLayerDefn();
    for (int i=0; i<ldefn->GetFieldCount(); i++) {
        OGRFieldDefn* fdefn = ldefn->GetFieldDefn(i);
        string name;
        name = fdefn->GetNameRef();
        fields.emplace_back(name);
    }
    return fields;
}

vector<string> GeoDataset::get_Fields_on_Layer(const char* layername){
    OGRLayer* poLayer;
    poLayer = poDS->GetLayerByName(layername);
    return get_Fields_on_Layer(poLayer);
}

void GeoDataset::remove_Field(const char* layername, const char* to_remove){
    OGRLayer* poLayer;
    poLayer = poDS->GetLayerByName(layername);
    OGRFeatureDefn* ldefn = poLayer->GetLayerDefn();
    int i = poLayer->FindFieldIndex(to_remove, true);
    poLayer->DeleteField(i);
}

vector<vector<string>> GeoDataset::get_layer_as_table(OGRLayer* poLayer) {
    vector<vector<string>> content;
    vector<string> line, titles;
    string value;
    OGRFeature* poFeature;

    // title line
    titles = get_Fields_on_Layer(poLayer);
    content.push_back(titles);
    // content
    while ((poFeature = poLayer->GetNextFeature()) != nullptr) {
        line.clear();
        for (auto &title : titles) {
            value = poFeature->GetFieldAsString(title.c_str());
            line.emplace_back(value);
        }
        content.push_back(line);
    }
    return content;
}

vector<vector<string>> GeoDataset::get_layer_as_table(const char* layername){
    OGRLayer* poLayer;
    poLayer = poDS->GetLayerByName(layername);
    return get_layer_as_table(poLayer);
}

void GeoDataset::stage_dataset() const {
    if (!(poDS->TestCapability(ODsCTransactions) || poDS->TestCapability(ODsCEmulatedTransactions))) {
        cerr << "Transaction for " << filename << " is not supported" << endl;
    }
    else if (poDS->StartTransaction() != OGRERR_NONE){
        cerr << "Transaction for " << filename << " cannot be started" << endl;
        exit(1);
    }
}

void GeoDataset::unstage_dataset() const {
    if(!(poDS->TestCapability(ODsCTransactions) || poDS->TestCapability(ODsCEmulatedTransactions))) {
        cerr << "Transaction for " << filename << " is not supported" << endl;
    }
    else if (poDS->RollbackTransaction() != OGRERR_NONE){
        cerr << "Transaction for " << filename << " cannot be rolled back" << endl;
        exit(1);
    }
}

void GeoDataset::commit_dataset() const {
    if(!(poDS->TestCapability(ODsCTransactions) || poDS->TestCapability(ODsCEmulatedTransactions))) {
        cerr << "Transaction for " << filename << " is not supported" << endl;
    }
    else if(poDS->CommitTransaction() != OGRERR_NONE){
        poDS->RollbackTransaction();
        cerr << "Transaction for " << filename << " cannot be committed" << endl;
        exit(1);
    }
}

void GeoDataset::vacuum_dataset() const {
    poDS->ExecuteSQL("VACUUM", nullptr, nullptr);
}

OGRLayer* GeoDataset::write_csv_as_layer(const char* inputfile, const char* layername, const char* fid,
                                         unsigned int skip, unsigned int Tln, char sep) const {

    vector<vector<string>> content = read_csv(inputfile, sep);
    vector<const char*> fieldNames;
    int fid_idx=-1;
    bool found=false;

    for (int i=0; i<content[Tln-1].size(); i++) {
        const char *field = content[Tln - 1][i].c_str();
        if (found || strcmp(field, fid) != 0)
            fieldNames.push_back(field);
        else {
            fid_idx = i;
            found = true;
        }
    }
    vector<Type> fieldTypes (fieldNames.size(), STR);

    // create layer
    OGRLayer* poLayer;
    poLayer = create_Layer(layername, fieldNames, fieldTypes, fid);
    bool transaction_is_active = (poLayer->StartTransaction() == OGRERR_NONE);    // true if transaction started

    // create features
    for (unsigned int i=skip; i<content.size(); i++) {      // loop over rows, skip header
        // Create a feature
        OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
        if (found)
            poFeature->SetFID(stoi(content[i][fid_idx]));

        for (int j=0; j<fieldNames.size(); j++) {
            int idx = poFeature->GetFieldIndex(fieldNames[j]);
            if (j<fid_idx)
                poFeature->SetField(idx, content[i][j].c_str());
            else
                poFeature->SetField(idx, content[i][j + (int)found].c_str());
        }
        if (poLayer->CreateFeature(poFeature ) != OGRERR_NONE ) {
            cerr << "Failed to create feature." << endl;
            exit( 1 );
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    if (transaction_is_active && poLayer->CommitTransaction() != OGRERR_NONE){
        poLayer->RollbackTransaction();
        cerr << "Transaction for " << inputfile << " cannot be committed" << endl;
        exit(1);
    }
    return poLayer;
}

//TODO read_csv: copy from QNetwork
vector<vector<string>> GeoDataset::read_csv(const char * inputfile, char sep) {
    vector<vector<string>> content = read_csv_if_exists(inputfile, sep);
    if (content.empty()) {
        cerr << "Could not open the file " << inputfile << "\n";
        exit(1);
    }
    return content;
}

vector<vector<string>> GeoDataset::read_csv_if_exists(const char * inputfile, char sep) {
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













/*** ------------------------- Child Class GeoReadDataset ------------------------------ ***/

// Constructor for input file using GDAL
GeoReadDataset::GeoReadDataset(const char* inputfile, const char* pointlayer, const char* tracklayer, bool verbose) :
                            GeoDataset(inputfile, READ, verbose) {
    trackLayerName = tracklayer;
    pointLayerName = pointlayer;
}

GeoReadDataset::GeoReadDataset(const char* inputfile, bool verbose, const char* pointlayer, const char* tracklayer) :
        GeoDataset(inputfile, READ, verbose) {
    trackLayerName = tracklayer;
    pointLayerName = pointlayer;
}


/**
 * Method to extract all tracks from the layer of the dataset
 *
 *It acts as a bridge between GDAL (Multi)LineString and GeoReadTrack
 *
 * The method ensures that all tracks are unique
 * @return vector<GeoReadTrack*>
 */

template<typename T>
std::vector<GeoReadTrack*> GeoReadDataset::extract_all_tracks(const char* filter) {
    std::vector<GeoReadTrack*> tracks;
    OGRFeature* poFeature;

    if (trackLayer == nullptr){
        // find trackLayer
        if ( (trackLayer = poDS->GetLayerByName(trackLayerName.c_str())) == nullptr) {
            cerr << "Track Layer is not available" << endl;
            exit(1);
        }
    }

    // set filter
    if (strcmp(filter, "") != 0)
        trackLayer->SetAttributeFilter(filter);

    while ((poFeature = trackLayer->GetNextFeature()) != nullptr) {
        OGRGeometry* poGeometry = poFeature->GetGeometryRef();
        auto* track = new T(poGeometry, poFeature);

        // check if the track is not a duplicate
        bool exists = false;
        for (const auto& existingTrack : tracks) {
            const double tol(0.01);
            if ((existingTrack->start->distance(*track->start) <= tol && existingTrack->end->distance(*track->end) <= tol) ||
                (existingTrack->start->distance(*track->end) <= tol && existingTrack->end->distance(*track->start) <= tol)) {
                exists = true;
                break;
            }
        }
        if (!exists)
            tracks.push_back(track);
    }
    trackLayer->ResetReading();

    // update the distance
    for (const auto& track : tracks) {
        track->compute_track();
        if (track->steps.empty())
            track->get_steps();
    }
    return tracks;
}
// explicit instantiations for template function:
template std::vector<GeoReadTrack*> GeoReadDataset::extract_all_tracks<GeoReadTrack>(const char* filter);
template std::vector<GeoReadTrack*> GeoReadDataset::extract_all_tracks<GeoRaccordTrack>(const char* filter);


/**
 * Method to extract all points from the layer of the dataset
 *
 * It acts as a bridge between GDAL (Multi)Point and GeoCoord
 *
 * @return vector<GeoCoord>
 */
std::vector<GeoCoord> GeoReadDataset::extract_all_points() {
    std::vector<GeoCoord> coords;
    OGRFeature* poFeature;

    if (pointLayer == nullptr){
        // find nodeLayer
        if ( (pointLayer = poDS->GetLayerByName(pointLayerName.c_str())) == nullptr) {
            cerr << "Node Layer is not available" << endl;
            exit(1);
        }
    }

    while ((poFeature = pointLayer->GetNextFeature()) != nullptr) {
        OGRGeometry* geometry = poFeature->GetGeometryRef();

        switch (wkbFlatten(geometry->getGeometryType())) {
            case wkbPoint: {
                OGRPoint* point = geometry->toPoint();
                if (point) {
                    coords.emplace_back(GeoCoord::CH1903plus, point->getX(), point->getY());
                }
                break;
            }
            case wkbMultiPoint: {
                OGRMultiPoint* multi_point = geometry->toMultiPoint();
                if (multi_point) {
                    for (int i = 0; i < multi_point->getNumGeometries(); i++) {
                        //OGRPoint* point = multi_point->getGeometryRef(i)->toPoint();          // deprecated
                        OGRPoint* point = multi_point->getGeometryRef(i);
                        if (point) {
                            coords.emplace_back(GeoCoord::CH1903plus, point->getX(), point->getY());
                        }
                    }
                }
                break;
            }
            default:
                std::cout << "Geometry type not supported. Expected 'Point' or 'MultiPoint'" << std::endl;
        }
        // Destroy feature regardless of geometry type
        OGRFeature::DestroyFeature(poFeature);
    }
    return coords;
}













/*** ------------------------- Child Class GeoWriteDataset ------------------------------ ***/

GeoWriteDataset::GeoWriteDataset(const char* outputfile, Accessmode rw,
                const char* nodelayer, const char* edgelayer, const char* sstlayer, bool verbose)
                : GeoDataset(outputfile, rw, verbose) {
    nodeLayerName = nodelayer;
    edgeLayerName = edgelayer;
    sstLayerName = sstlayer;
}

GeoWriteDataset::GeoWriteDataset(const char* outputfile, const char* nodelayer, const char* edgelayer,
                                 const char* sstlayer, Accessmode rw, bool verbose)
        : GeoDataset(outputfile, rw, verbose) {
    nodeLayerName = nodelayer;
    edgeLayerName = edgelayer;
    sstLayerName = sstlayer;
}

GeoWriteDataset::GeoWriteDataset(const char* outputfile, Accessmode rw, bool verbose,
                         const char* nodelayer, const char* edgelayer, const char* sstlayer)
                         : GeoDataset(outputfile, rw, verbose) {
    nodeLayerName = nodelayer;
    edgeLayerName = edgelayer;
    sstLayerName = sstlayer;
}

GeoWriteDataset::GeoWriteDataset(const char* outputfile, bool verbose, Accessmode rw,
                const char* nodelayer, const char* edgelayer, const char* sstlayer)
                : GeoDataset(outputfile, rw, verbose) {
    nodeLayerName = nodelayer;
    edgeLayerName = edgelayer;
    sstLayerName = sstlayer;
}


void GeoWriteDataset::set_Layers(){
    if (!nodeLayerName.empty() && (nodeLayer = poDS->GetLayerByName(nodeLayerName.c_str())) == nullptr)
        cerr << "Node Layer is not defined" << endl;
    if (!edgeLayerName.empty() && (edgeLayer = poDS->GetLayerByName(edgeLayerName.c_str())) == nullptr)
        cerr << "Edge Layer is not defined" << endl;
    if (!sstLayerName.empty() && (sstLayer  = poDS->GetLayerByName( sstLayerName.c_str())) == nullptr)
        cerr << "SST Layer is not defined" << endl;
}

void GeoWriteDataset::createNodeLayer(const vector<const char*>& fieldNames, const GeoCoord::CRS crs){
    vector<Type> fieldTypes (fieldNames.size());
    for (int i = 0; i < fieldNames.size(); i++) {
        if ( strcmp(fieldNames[i], "SST ID") == 0 )
            fieldTypes[i] = INT;
        else
            fieldTypes[i] = STR;
    }
    nodeLayer = create_Layer(nodeLayerName.c_str(), fieldNames, fieldTypes, "fid", crs);
}

void GeoWriteDataset::createEdgeLayer(const vector<const char*>& fieldNames, const GeoCoord::CRS crs) {
    vector<Type> fieldTypes (fieldNames.size());
    for (int i = 0; i < fieldNames.size(); i++) {
        if (strcmp(fieldNames[i], "DN") == 0 || strcmp(fieldNames[i], "Branch (l)") == 0 || strcmp(fieldNames[i], "Segment (k)") == 0)
            fieldTypes[i] = INT;
        else if (strcmp(fieldNames[i], "Length") == 0)
            fieldTypes[i] = DOUBLE;
        else
            fieldTypes[i] = STR;
    }
    edgeLayer = create_Layer(edgeLayerName.c_str(), fieldNames, fieldTypes, "fid", crs);
}

void GeoWriteDataset::createSSTLayer(const vector<const char*>& fieldNames, const GeoCoord::CRS crs){
    vector<Type> fieldTypes (fieldNames.size());
    for (int i = 0; i < fieldNames.size(); i++) {
        if (strcmp(fieldNames[i], "DN_raccord") == 0 )
            fieldTypes[i] = INT;
        else if (strcmp(fieldNames[i], "Altitude [m]") == 0 ||
                strcmp(fieldNames[i], "T_need (°C)") == 0 ||
                strcmp(fieldNames[i], "T_ref (°C)") == 0 ||
                strcmp(fieldNames[i], "L_raccord (m)") == 0)
            fieldTypes[i] = DOUBLE;
        else
            fieldTypes[i] = STR;        // TODO option "data"

        sstLayer = create_Layer(sstLayerName.c_str(), fieldNames, fieldTypes, "fid", crs);
    }
}

//void GeoWriteDataset::write_point(GeoCoord* coord, const char* field, const char* value){
//    OGRSpatialReference* oSRS;
//    oSRS->SetWellKnownGeogCS( "WGS84" );
//    pointsLayer->SetActiveSRS(0, oSRS);
//
//    // Create a feature
//    OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn() );
//    poFeature->SetField( field, value );
//
//    OGRPoint pt;
//    pt.setX( coord->lat );
//    pt.setY( coord->lon );
//
//    poFeature->SetGeometry( &pt );
//
//    if(poLayer->CreateFeature(poFeature ) != OGRERR_NONE ) {
//        printf( "Failed to create feature.\n" );
//        exit( 1 );
//    }
//    OGRFeature::DestroyFeature( poFeature );
//}
//
//void GeoWriteDataset::write_points(const vector<GeoCoord*>& coords, const vector<const char*>& values){
//    for (int i=0; i<coords.size(); i++) {
////        if (values.size() < coords.size())
////            write_point(coords[i]);     // use default values
////        else
////            write_point(coords[i], values[i]);
//    }
//}

void GeoWriteDataset::write_nodes(const vector<QNode*>& nodes) const{
    bool transaction_is_active = (nodeLayer->StartTransaction() == OGRERR_NONE);    // true if transaction started

    for (auto node : nodes) {
        // Create a feature
        OGRFeature *poFeature = OGRFeature::CreateFeature(nodeLayer->GetLayerDefn() );
        poFeature->SetFID(node->id);
        int idx;

        if (node->type == QNode::QNodeType::PROTON) {
            if ((idx=poFeature->GetFieldIndex("Type")) != -1)   // check if field exists
                poFeature->SetField(idx, "Proton");
            if ((idx=poFeature->GetFieldIndex("SST ID")) != -1)
                poFeature->SetField(idx, (int) node->sst_id);
            if ((idx=poFeature->GetFieldIndex("Branch (l)")) != -1)
                poFeature->SetField(idx, (int) node->l);
            if ((idx=poFeature->GetFieldIndex("Segment (k)")) != -1)
                poFeature->SetField(idx, (int) node->k);
        }
        else if (node->type == QNode::QNodeType::NEUTRON) {
            if ((idx=poFeature->GetFieldIndex("Type")) != -1)
                poFeature->SetField(idx, "Neutron");
            if ((idx=poFeature->GetFieldIndex("Neutron (n)")) != -1)
                poFeature->SetField(idx, (int)node->n);
            if ((idx=poFeature->GetFieldIndex("Branch (l)")) != -1)
                poFeature->SetField(idx, (int) node->l);
        }
        else if (node->type == QNode::QNodeType::CENTRAL)
            if ((idx=poFeature->GetFieldIndex("Type")) != -1)
                poFeature->SetField(idx, "Central Plant");

//        OGRSpatialReference oSRS;
//        oSRS.SetWellKnownGeogCS( "WGS84" );        //Swiss CRS = EPSG:2056
//        nodeLayer->SetActiveSRS(0, &oSRS);

        if (strcmp(nodeLayer->GetSpatialRef()->GetName(), "WGS 84 (CRS84)") != 0)
            cerr << "Not yet implemented. WGS84 is used." << endl;

        OGRPoint pt;
        pt.setX( node->coord->lon );
        pt.setY( node->coord->lat );

        poFeature->SetGeometry( &pt );

        if(nodeLayer->CreateFeature(poFeature ) != OGRERR_NONE ) {
            printf( "Failed to create feature.\n" );
            exit( 1 );
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    if (transaction_is_active && nodeLayer->CommitTransaction() != OGRERR_NONE){
        nodeLayer->RollbackTransaction();
        cerr << "Transaction for Nodes cannot be committed" << endl;
        exit(1);
    }
}

void GeoWriteDataset::write_node(QNode* node) const {
    vector<QNode*> nodes = {node};
    write_nodes(nodes);
}

void GeoWriteDataset::write_nodes(const QGraph& graph) const{
    vector<QNode*> nodes = {graph.source};
    nodes.insert(nodes.end(), graph.nodes.begin(), graph.nodes.end());
    write_nodes(nodes);
}

void GeoWriteDataset::write_edges(const vector<QEdge*>& edges) const {
    bool transaction_is_active = (edgeLayer->StartTransaction() == OGRERR_NONE);    // true if transaction started
    int idx;
    for (auto &edge : edges) {
        OGRFeature *poFeature = OGRFeature::CreateFeature(edgeLayer->GetLayerDefn());
        if ((idx=poFeature->GetFieldIndex("DN")) != -1)
            poFeature->SetField(idx, (int)edge->DN);
        if ((idx=poFeature->GetFieldIndex("Length")) != -1)
            poFeature->SetField(idx, edge->length);
        if ((idx=poFeature->GetFieldIndex("Branch (l)")) != -1)
            poFeature->SetField(idx, (int)edge->l);
        if ((idx=poFeature->GetFieldIndex("Segment (k)")) != -1)
            poFeature->SetField(idx, (int)edge->k);

        if (strcmp(nodeLayer->GetSpatialRef()->GetName(), "WGS 84 (CRS84)") != 0)
            cerr << "Not yet implemented. WGS84 is used." << endl;

        OGRLineString lstr;
        for (auto &step : edge->track->steps)
            lstr.addPoint(step->lon, step->lat);
        poFeature->SetGeometry(&lstr);

        if (edgeLayer->CreateFeature(poFeature) != OGRERR_NONE) {
            printf("Failed to create feature.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);
    }
    if (transaction_is_active && edgeLayer->CommitTransaction() != OGRERR_NONE){
        edgeLayer->RollbackTransaction();
        cerr << "Transaction for Edges cannot be committed" << endl;
        exit(1);
    }
}

void GeoWriteDataset::write_all_edges(const QGraph& graph) const{
    write_edges(graph.edges);
}

void GeoWriteDataset::write_mst_edges(const QGraph& graph) const{
    write_edges(graph.mst_edges);
}

void GeoWriteDataset::copy_file_to_sstLayer(const char* inputfile, int nTln, char sep) const {
    if (sstLayer == nullptr){
        cerr << "SST Layer is not defined" << endl;
        exit(1);
    }
    vector<vector<string>> content = read_csv(inputfile, sep);
    vector<string> fields = get_Fields_on_Layer(sstLayer);
    string field;
    field = "Substation nb";
    fields.emplace_back(field);
    map<const char*, short> columns;

    // get indices
    for (const auto& fieldName : fields){
        ptrdiff_t pos = find(content[0].begin(), content[0].end(), fieldName) - content[0].begin();
        if(pos >= content[0].size()) {
            cerr << "Column " << fieldName << " not found on SST file" << endl;
            //exit(1);
        }
        else
            columns[fieldName.c_str()] = (short) pos;
    }

    bool transaction_is_active = (sstLayer->StartTransaction() == OGRERR_NONE);    // true if transaction started
    for (int i=nTln; i<content.size(); i++) {      // loop over lines, skip header
        // Create a feature
        OGRFeature *poFeature = OGRFeature::CreateFeature(sstLayer->GetLayerDefn());
        poFeature->SetFID(stoi(content[i][columns["Substation nb"]]));

        for (const auto& column : columns) {
            int idx = poFeature->GetFieldIndex(column.first);
            field = column.first;
            if ( field == "Substation nb")
                continue;
            else if (field == "DN_raccord")
                poFeature->SetField(idx, stoi(content[i][column.second]));
            else if (field == "Altitude [m]" || field == "T_need (°C)" || field == "T_ref (°C)" || field == "L_raccord (m)")
                poFeature->SetField(idx, stod(content[i][column.second]));
            else
                poFeature->SetField(idx, content[i][column.second].c_str());     // TODO option "data": use AlterFieldDefn()
        }
        if (sstLayer->CreateFeature(poFeature ) != OGRERR_NONE ) {
            cerr << "Failed to create feature." << endl;
            exit( 1 );
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    if (transaction_is_active && sstLayer->CommitTransaction() != OGRERR_NONE){
        sstLayer->RollbackTransaction();
        cerr << "Transaction for SST cannot be committed" << endl;
        exit(1);
    }
}

void GeoWriteDataset::copy_layer_to_sstLayer(const char* infile, const char* inlayer) const {
    if (sstLayer == nullptr){
        cerr << "SST Layer is not defined" << endl;
        exit(1);
    }
    GeoDataset inDataset (infile, GeoDataset::Accessmode::READ);
    OGRLayer* inLayer = inDataset.get_Layer(inlayer);
    int fieldCount = sstLayer->GetLayerDefn()->GetFieldCount();
    string field;
    map<string, int> field_number;

    for (int i=0; i < fieldCount; i++){
        field = sstLayer->GetLayerDefn()->GetFieldDefn(i)->GetNameRef();
        int idx = inLayer->FindFieldIndex(field.c_str(), TRUE);
        if (idx >= 0)
            field_number[field] = idx;
        else
            cerr << "Field " << field << " not found on " << infile << endl;
    }

    // create and copy features
    bool transaction_is_active = (sstLayer->StartTransaction() == OGRERR_NONE);    // true if transaction started
    for (auto& inFeature: *inLayer) {      // loop over features, skip header
        // Create a feature
        OGRFeatureDefn* poLayerDefn = sstLayer->GetLayerDefn();
        OGRFeature *poFeature = OGRFeature::CreateFeature(poLayerDefn);
        poFeature->SetFID(inFeature->GetFID());

        for (int i=0; i < poLayerDefn->GetFieldCount(); i++) {  // loop over fields of sstLayer
            field = poLayerDefn->GetFieldDefn(i)->GetNameRef();

            if (field == "Substation nb")
                continue;
            else if (field == "DN_raccord")
                poFeature->SetField(i, inFeature->GetFieldAsInteger(field_number[field]));
            else if (field == "Altitude [m]" || field == "T_need (°C)" || field == "T_ref (°C)" || field == "L_raccord (m)")
                poFeature->SetField(i, inFeature->GetFieldAsDouble(field_number[field]));
            else
                poFeature->SetField(i, inFeature->GetFieldAsString(field_number[field]));     // TODO option "data": use AlterFieldDefn()
        }
        if (sstLayer->CreateFeature(poFeature) != OGRERR_NONE) {
            cerr << "Failed to create feature." << endl;
            exit( 1 );
        }
        OGRFeature::DestroyFeature( poFeature );
    }
    if (transaction_is_active && sstLayer->CommitTransaction() != OGRERR_NONE){
        sstLayer->RollbackTransaction();
        cerr << "Transaction for SST cannot be committed" << endl;
        exit(1);
    }
}

void GeoWriteDataset::copy_geometry_to_sstLayer(){
    int sst;
    bool transaction_is_active = (sstLayer->StartTransaction() == OGRERR_NONE);    // true if transaction started
    for(auto& nodeFeature: nodeLayer) {   // loop over all features in nodeLayer
        if (strcmp(nodeLayer->GetSpatialRef()->GetName(), "WGS 84 (CRS84)") != 0)
            cerr << "Not yet implemented. WGS84 is used." << endl;

        if ((sst=nodeFeature->GetFieldAsInteger("SST ID")) != 0){
            OGRFeature *sstFeature = sstLayer->GetFeature(sst);

            // copy position of corresponding proton
            OGRPoint *pt = nodeFeature->GetGeometryRef()->toPoint();
            sstFeature->SetGeometry(pt);

            if (sstLayer->UpsertFeature(sstFeature) != OGRERR_NONE) {
                printf("Failed to upsert feature.\n");
                exit(1);
            }
            OGRFeature::DestroyFeature(sstFeature);
        }
    }
    if (transaction_is_active && sstLayer->CommitTransaction() != OGRERR_NONE){
        sstLayer->RollbackTransaction();
        cerr << "Transaction for Edges cannot be committed" << endl;
        exit(1);
    }
}

void GeoWriteDataset::write_tracks_to_sstLayer(const vector<GeoReadTrack*>& raccord_track, double tol) {
    // create DN_raccord, L_raccord fields if not yet on Layer
    if (sstLayer->FindFieldIndex("DN_raccord", TRUE) == -1) {
        OGRFieldDefn oField( "DN_raccord", OFTInteger );
        oField.SetWidth(32);
        if( sstLayer->CreateField( &oField ) != OGRERR_NONE ) {
            printf( "Creating DN_raccord field failed.\n" );
            exit( 1 );
        }
    }
    if (sstLayer->FindFieldIndex("L_raccord (m)", TRUE) == -1) {
        OGRFieldDefn oField( "L_raccord", OFTReal );
        oField.SetWidth(32);
        if( sstLayer->CreateField( &oField ) != OGRERR_NONE ) {
            printf( "Creating L_raccord field failed.\n" );
            exit( 1 );
        }
    }

    bool transaction_is_active = (edgeLayer->StartTransaction() == OGRERR_NONE);    // true if transaction started
    // get DN_raccord, L_raccord from raccord_tracks
    for (auto &track : raccord_track){
        // identify feature
        OGRFeature *poFeature;
        auto* raccord = dynamic_cast<GeoRaccordTrack*>(track);
        if (raccord != nullptr && raccord->sst_id > 0) {
            // get or create feature
            if ((poFeature = sstLayer->GetFeature(raccord->sst_id)) == nullptr) {
                poFeature = OGRFeature::CreateFeature(sstLayer->GetLayerDefn());
                poFeature->SetFID(raccord->sst_id);
            }
        }
        else {  // no sst_id available, get from nodeLayer
            int sst;
            for (auto& nodeFeature: nodeLayer) {   // loop over all features in nodeLayer
                if (strcmp(nodeFeature->GetFieldAsString("Type"), "Proton") == 0){
                    OGRGeometry *geom = nodeFeature->GetGeometryRef();

                    if (strcmp(geom->getSpatialReference()->GetName(), "WGS 84 (CRS84)") != 0)
                        cerr << "Not yet implemented. WGS84 is used." << endl;
                    GeoCoord proton(GeoCoord::CRS::WGS84, geom->toPoint()->getX(), geom->toPoint()->getY());

                    if (distance_2d(proton, *track->start) < tol || distance_2d(proton, *track->end) < tol){
                        sst = nodeFeature->GetFieldAsInteger("SST ID");
                        poFeature = sstLayer->GetFeature(sst);
                        break;
                    }
                }
            }
        }

        // fill values
        int idx;
        if ((idx=poFeature->GetFieldIndex("DN_raccord")) != -1)
            poFeature->SetField(idx, (int)track->diameter);
        if ((idx=poFeature->GetFieldIndex("L_raccord (m)")) != -1)
            poFeature->SetField(idx, track->length);

        if (strcmp(nodeLayer->GetSpatialRef()->GetName(), "WGS 84 (CRS84)") != 0)
            cerr << "Not yet implemented. WGS84 is used." << endl;
        OGRLineString lstr;
        for (auto &step : track->steps)
            lstr.addPoint(step->lon, step->lat);
        poFeature->SetGeometry(&lstr);

        if (sstLayer->UpsertFeature(poFeature) != OGRERR_NONE) {
            printf("Failed to create feature.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);
    }
    if (transaction_is_active && sstLayer->CommitTransaction() != OGRERR_NONE){
        sstLayer->RollbackTransaction();
        cerr << "Transaction for Edges cannot be committed" << endl;
        exit(1);
    }
}

void GeoWriteDataset::extract_networktable(vector<vector<int>>& table) {
    int k, l, kmax=-1, lmax=-1;
    map<pair<int, int>, int> protons;   // (k,l) and SST-ID
    map<int, int> neutrons;             // l and n
    OGRFeature* poFeature;

    if (nodeLayer == nullptr){
        // find nodeLayer
        if ( (nodeLayer = poDS->GetLayerByName(nodeLayerName.c_str())) == nullptr) {
            cerr << "Node Layer is not available" << endl;
            exit(1);
        }
    }

    // check if all relevant fields exist
    vector<string> fields = get_Fields_on_Layer(nodeLayer);
    for (const char* field : {"Type", "Branch (l)", "Segment (k)", "Neutron (n)", "SST ID"})
        if (find(fields.begin(), fields.end(), field) == fields.end()){
            cerr << "Field " << field << " is missing in file " << filename;
            cerr << " on layer " << nodeLayer->GetName() << endl;
            exit(1);
        }

    // read features
    while ((poFeature = nodeLayer->GetNextFeature()) != nullptr) {
        // get protons
        if (strcmp(poFeature->GetFieldAsString("Type"), "Proton") == 0) {
            l = poFeature->GetFieldAsInteger("Branch (l)");
            k = poFeature->GetFieldAsInteger("Segment (k)");
            if (k > 0) {
                protons[make_pair(k, l)] = poFeature->GetFieldAsInteger("SST ID");
                if (abs(l) > lmax)
                    lmax = abs(l);
                if (k > kmax)
                    kmax = k;
            }
            else {
                cerr << "Something wrong: k not set for proton" << endl;
                exit(1);
            }
        }
            // get neutrons
        else if (strcmp(poFeature->GetFieldAsString("Type"), "Neutron") == 0) {
            l = poFeature->GetFieldAsInteger("Branch (l)");
            neutrons[l] = poFeature->GetFieldAsInteger("Neutron (n)");
            if (neutrons[l] <= 0) {
                cerr << "Something wrong: n not set for neutron" << endl;
                exit(1);
            }
        }
    }

    // fill table
    table.clear();
    table.resize(kmax+2, vector<int>(2*lmax+1));  // filled with default value 0
    for (int i=-lmax; i<=lmax; i++)
        table[0][i+lmax] = i;          // fill first line with branch numbers
    for (auto p : protons)
        table[p.first.first][p.first.second + lmax] = p.second;
    for (auto n : neutrons)
        table[kmax+1][n.first + lmax] = n.second;
}

void GeoWriteDataset::extract_diam_length_table(vector<vector<int>>& diamtable, vector<vector<double>>& lengthtable) {
    int k, l, kmax=-1, lmax=-1;
    map<pair<int, int>, pair<int, double>> data;
    OGRFeature* poFeature;

    if (edgeLayer == nullptr){
        // find edgeLayer
        if ( (edgeLayer = poDS->GetLayerByName(edgeLayerName.c_str())) == nullptr) {
            cerr << "Edge Layer is not available" << endl;
            exit(1);
        }
    }

    // check if all relevant fields exist
    vector<string> fields = get_Fields_on_Layer(edgeLayer);
    for (const char* field : {"Branch (l)", "Segment (k)", "DN", "Length"})
        if (find(fields.begin(), fields.end(), field) == fields.end()){
            cerr << "Field " << field << " is missing in file " << filename;
            cerr << " on layer " << edgeLayer->GetName() << endl;
            exit(1);
        }

    // read features
    while ((poFeature = edgeLayer->GetNextFeature()) != nullptr) {
        l = poFeature->GetFieldAsInteger("Branch (l)");
        k = poFeature->GetFieldAsInteger("Segment (k)");
        if (k>0) {
            data[make_pair(k, l)] = make_pair(poFeature->GetFieldAsInteger("DN"),
                                              poFeature->GetFieldAsDouble("Length"));
            if (abs(l) > lmax)
                lmax = abs(l);
            if (k > kmax)
                kmax = k;
        }
        else {
            cerr << "Something wrong: k not set for segment" << endl;
            exit(1);
        }
    }

    // fill tables
    diamtable.clear();
    diamtable.resize(kmax, vector<int>(2*lmax+1));  // default value 0
    lengthtable.clear();
    lengthtable.resize(kmax, vector<double>(2*lmax+1));
    for (auto d : data) {
        diamtable[d.first.first - 1][d.first.second + lmax] = d.second.first;
        lengthtable[d.first.first - 1][d.first.second + lmax] = d.second.second;
    }
}





/***  Outside of classes ***/
void extract_quantum_network(const char* inputfile, const char* outputfile, const char* sstfile,
                             const char* pointlayer, const char* tracklayer,
                             const char* nodelayer, const char* edgelayer, const char* sstlayer,
                             double tol, const char* distrib_filter, const char* raccord_filter) {
    GeoReadDataset readDataset(inputfile, pointlayer, tracklayer);
    vector<GeoReadTrack*> distrib_tracks (readDataset.extract_all_tracks<GeoReadTrack>(distrib_filter));
    vector<GeoReadTrack*> raccord_tracks (readDataset.extract_all_tracks<GeoRaccordTrack>(raccord_filter));

    // merge tracks, but do not mix distrib and raccord
    merge_all(distrib_tracks, tol);
    merge_all(raccord_tracks, tol);

    segment_tracks(distrib_tracks, tol);                           // Segment tracks with themselves --> get branches
    segment_tracks(distrib_tracks, raccord_tracks, tol);     // Segment tracks with raccordements --> get segments

    for (auto &track : raccord_tracks)  // TODO
        dynamic_cast<GeoRaccordTrack *>(track)->get_sst_id();

    QGraph graph(distrib_tracks, tol);
    graph.identify_source_and_sst(raccord_tracks, tol);
    graph.define_mst(graph.edges);
    graph.compute_quantum_numbers(true);


    GeoWriteDataset writeDataset(outputfile, nodelayer, edgelayer, sstlayer);
    writeDataset.createNodeLayer();
    writeDataset.write_nodes(graph);
    writeDataset.createEdgeLayer();
    writeDataset.write_mst_edges(graph);

    // get SST from gpkg or csv
    size_t length = strlen(sstfile);
    if (strncmp(sstfile + length - 4, ".csv", 4) == 0) {       // csv
        writeDataset.createSSTLayer();
        writeDataset.copy_file_to_sstLayer(sstfile);
    }
    else if (strncmp(sstfile + length - 5, ".gpkg", 5) == 0) {   // gpkg
        if (strcmp(outputfile, sstfile) != 0) {      // copy relevant info from sstfile to outputfile
            writeDataset.createSSTLayer();
            writeDataset.copy_layer_to_sstLayer(sstfile, sstlayer);
        }
    }
    else {
        cerr << "SST File must be either .csv or .gpkg" << endl;
        exit(1);
    }
    writeDataset.write_tracks_to_sstLayer(raccord_tracks);
}

void extract_quantum_network(const char* inputfile, const char* outputfile, const char* sstfile,
                             const char* pointlayer, const char* tracklayer, const char* sstlayer, const double tol,
                             const char* distrib_filter, const char* raccord_filter){
    extract_quantum_network(inputfile, outputfile, sstfile, pointlayer, tracklayer,
                            "Nodes", "Edges", sstlayer, tol, distrib_filter, raccord_filter);
}

void extract_quantum_network(const char* inputfile, const char* outputfile, const char* sstfile, const double tol,
                             const char* pointlayer, const char* tracklayer, const char* sstlayer,
                             const char* distrib_filter, const char* raccord_filter){
    extract_quantum_network(inputfile, outputfile, sstfile, pointlayer, tracklayer, sstlayer, tol, distrib_filter, raccord_filter);
}