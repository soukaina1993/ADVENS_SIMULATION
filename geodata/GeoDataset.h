//
// Created by jamil.maj on 19.09.2024.
//


#ifndef GEOREADDATASET_H
#define GEOREADDATASET_H

#include <ogrsf_frmts.h>
#include <vector>

#include "GeoReadTrack.h"
#include "QGraph.h"


/*** Base Class GeoDataset ***/
class GeoDataset {
public:
    enum Type { INT, DOUBLE, STR };
    enum Accessmode { READ, WRITE };

public:
    GeoDataset()=default;
    explicit GeoDataset(const char* ifilename, Accessmode rw, bool verbose = false);
    ~GeoDataset();

public:
    OGRLayer* create_Layer(const char* name,
                           const vector<const char*>& fieldNames, const vector<Type>& fieldTypes, const char* fid="fid",
                           GeoCoord::CRS crs=GeoCoord::CRS::WGS84) const;

    OGRLayer* create_Layer(const char* name,
                           const vector<const char*>& fieldNames, Type fieldType, const char* fid="fid",
                           GeoCoord::CRS crs=GeoCoord::CRS::WGS84) const;

    OGRLayer* create_Layer(const char* name,
                           const vector<string>& fieldNames, const vector<Type>& fieldTypes, const char* fid="fid",
                           GeoCoord::CRS crs=GeoCoord::CRS::WGS84) const;

    OGRLayer* create_Layer(const char* name,
                           const vector<string>& fieldNames, Type fieldType, const char* fid="fid",
                           GeoCoord::CRS crs=GeoCoord::CRS::WGS84) const;

    vector<string> get_Layers();
    OGRLayer* get_Layer(const char* layername);
    void remove_Layer(const char* to_remove);

    static vector<string> get_Fields_on_Layer(OGRLayer* poLayer);
    vector<string> get_Fields_on_Layer(const char* layername);
    void remove_Field(const char* layername, const char* to_remove);

    static vector<vector<string>> get_layer_as_table(OGRLayer* poLayer);
    vector<vector<string>> get_layer_as_table(const char* layername);

    // working with transactions
    void stage_dataset() const;
    void unstage_dataset() const;
    void commit_dataset() const;

    // remove unnecessary data
    void vacuum_dataset() const;


    /**
    * write a csv into a layer:
    * - the first 'skip' lines are ignored
    * - fid is found in column fid
    * - field names are taken from line Tln
    **/
    OGRLayer* write_csv_as_layer(const char* inputfile, const char* layername, const char* fid="fid",
                                 unsigned int skip=1, unsigned int Tln=1, char sep='\t') const;

    static vector<vector<string>> read_csv(const char* inputfile, char sep='\t');
    static vector<vector<string>> read_csv_if_exists(const char* inputfile, char sep='\t');

protected:
    string filename;
    GDALDataset* poDS=nullptr;
    GDALDriver* poDriver=nullptr;
};


/**
 *  GeoReadDataset intended to read a dataset from a .shp or .gpkg file using GDAL
 *  The class adds the ability to extract all lines from the layer of the dataset directly as a vector of OGRLineString
 *
 **/
class GeoReadDataset : public GeoDataset {
public: // Attributes
    // Constructor
    GeoReadDataset()=default;
    GeoReadDataset(const char* inputfile, const char* pointLayerName, const char* trackLayerName, bool verbose = false);
    explicit GeoReadDataset(const char* inputfile, bool verbose = false,
                   const char* pointLayerName="nodes", const char* trackLayerName="edges");

private:
    OGRLayer* trackLayer=nullptr;
    OGRLayer* pointLayer=nullptr;
    string trackLayerName, pointLayerName;

public:
    // T = GeoReadTrack or GeoRaccordTrack
    template<typename T> std::vector<GeoReadTrack*> extract_all_tracks(const char* filter="");
    vector<GeoCoord> extract_all_points();
};

/**
 *  GeoWriteDataset intended to write a dataset into a .gpkg file using GDAL
 *
 **/
class GeoWriteDataset : public GeoDataset {
public:
    // Constructor
    GeoWriteDataset()=default;
    GeoWriteDataset(const char* outputfile, Accessmode rw,
                    const char* nodelayer, const char* edgelayer, const char* sstlayer, bool verbose=false);
    GeoWriteDataset(const char* outputfile, const char* nodelayer, const char* edgelayer, const char* sstlayer,
                    Accessmode rw=Accessmode::WRITE, bool verbose=false);
    explicit GeoWriteDataset(const char* outputfile, Accessmode rw=Accessmode::WRITE, bool verbose=false,
                             const char* nodelayer="Nodes", const char* edgelayer="Edges", const char* sstlayer="SST");
    GeoWriteDataset(const char* outputfile, bool verbose, Accessmode rw=Accessmode::WRITE,
                    const char* nodelayer="Nodes", const char* edgelayer="Edges", const char* sstlayer="SST");


public:
    // set layer pointers from layer names
    void set_Layers();

    // create layers
    void createNodeLayer(const vector<const char*>& fieldNames=
            {"Type", "Neutron (n)", "Branch (l)", "Segment (k)", "SST ID"},
            GeoCoord::CRS crs=GeoCoord::CRS::WGS84);

    void createEdgeLayer(const vector<const char*>& fieldNames=
            {"DN", "Length", "Branch (l)", "Segment (k)"},
            GeoCoord::CRS crs=GeoCoord::CRS::WGS84);

    void createSSTLayer(const vector<const char*>& fieldNames=
            {"Energy [kWh/an]",	"Substation Type", "Substation Name", "Building Type", "Heat Type",
             "Altitude [m]", "T_need (C)", "T_ref (C)", "dT", "L_raccord (m)", "DN_raccord"},
             GeoCoord::CRS crs=GeoCoord::CRS::WGS84);

//    // not implemented
//    void write_point(GeoCoord* coord, const char* field, const char* value="-");
//    void write_points(const vector<GeoCoord*>& coords, const vector<const char*>& values={});
//    void write_point(GeoCoord* coord, int value);
//    void write_points(const vector<GeoCoord*>& coords, const vector<int>& values);
//    void write_point(GeoCoord* coord, double value);
//    void write_points(const vector<GeoCoord*>& coords, const vector<double>& values);

    // write functions
    void write_nodes(const vector<QNode*>& nodes) const;
    void write_node(QNode* node) const;
    void write_nodes(const QGraph& graph) const;      // writes source and all nodes of graph
    void write_edges(const vector<QEdge*>& edges) const;
    void write_all_edges(const QGraph& graph) const;    // writes graph.edges (potentially no (l,k) values)
    void write_mst_edges(const QGraph& graph) const;    // writes graph.mst_edges

    void copy_file_to_sstLayer(const char* inputfile, int nTln=1, char sep='\t') const;     // get all info from file
    void copy_layer_to_sstLayer(const char* infile, const char* inlayer) const; // get info from another layer, only works with same field names
    void copy_geometry_to_sstLayer();                               // copies the position of the corresponding proton
    void write_tracks_to_sstLayer(const vector<GeoReadTrack*>& raccord_track, double tol= 0.1);

    // read functions
    void extract_networktable(vector<vector<int>>& table);
    void extract_diam_length_table(vector<vector<int>>& diamtable, vector<vector<double>>& lengthtable);


public:
    OGRLayer* nodeLayer=nullptr;
    OGRLayer* edgeLayer=nullptr;
    OGRLayer* sstLayer=nullptr;
private:
    string nodeLayerName, edgeLayerName, sstLayerName;
};



/*** --------------------------------- Outside of classes -------------------------------------------------------------- ***/

/**
 * reads in a cleaned DH network and extracts its quantum network representation
 * - inputfile: the input gpkg file
 * - outputfile: the output gpkg file, default values are used for the layers
 * - sstfile: either a csv or gpkg file
 * - pointlayer: points in inputfile
 * - tracklayer: edges in inputfile
 * - sstlayer: if csv: not used / if gpkg: existing SST layer
 * - tol: tolerance for merging
 * - distrib_filter/ raccord_filter: filter to distinguish track type
**/
void extract_quantum_network(const char* inputfile, const char* outputfile, const char* sstfile,
                             const char* pointlayer, const char* tracklayer,
                             const char* nodelayer, const char* egelayer, const char* sstlayer,
                             double tol=0.1,
                             const char* distrib_filter="Genre_cond = 'Retour' AND Fonction = 'Conduite de distribution'",
                             const char* raccord_filter="Genre_cond = 'Retour' AND Fonction = 'Conduite de raccordement'");

void extract_quantum_network(const char* inputfile, const char* outputfile, const char* sstfile,
                             const char* pointlayer="nodes", const char* tracklayer="edges", const char* sstlayer="SST",
                             double tol=0.1,
                             const char* distrib_filter="Genre_cond = 'Retour' AND Fonction = 'Conduite de distribution'",
                             const char* raccord_filter="Genre_cond = 'Retour' AND Fonction = 'Conduite de raccordement'");

void extract_quantum_network(const char* inputfile, const char* outputfile, const char* sstfile, double tol,
                             const char* pointlayer="nodes", const char* tracklayer="edges", const char* sstlayer="SST",
                             const char* distrib_filter="Genre_cond = 'Retour' AND Fonction = 'Conduite de distribution'",
                             const char* raccord_filter="Genre_cond = 'Retour' AND Fonction = 'Conduite de raccordement'");


#endif //GEOREADDATASET_H
