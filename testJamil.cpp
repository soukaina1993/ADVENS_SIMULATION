#include "geodata/GeoReadTrack.h"
#include "geodata/QGraph.h"
#include <iostream>

#include "geodata/GeoDataset.h"
// #include "geodata/HexGrid.h"

#include "network/QNetwork.h"


int main() {

    // std::string NetworkFilePath("../network_data/Network/V_SHP_CAD_PIPE_CAD_ESTA_20240807.shp"); // test SHP format
    // std::string NetworkFilePath = "../network_data/Network/V_SHP_CAD_PIPE_CAD_ESTA_20240807.geojson"; // test GeoJSON format

    // in order to use local paths, specify the working directory in Clion's configuration settings
    const char* NetworkFilePath = "network_data/GPKG/PreprocessedNetwork.gpkg"; // test GPKG format

    // subNetwork
    // NetworkFilePath = "../network_data/SubNetwork/subNet.shp"; // test subNetwork SHP
    // NetworkFilePath = "../network_data/SubNetwork/subNet.geojson"; // test subNetwork geojson

    // Initialize a GeoReadDataset object
    GeoReadDataset readDataset(NetworkFilePath, "nodes", "edges", true);

    std::vector<GeoReadTrack*>
            distrib_tracks(readDataset.extract_all_tracks<GeoReadTrack>("Genre_cond = 'Retour' AND Fonction = 'Conduite de distribution'"));
    std::vector<GeoReadTrack*>
            raccord_tracks(readDataset.extract_all_tracks<GeoRaccordTrack>("Genre_cond = 'Retour' AND Fonction = 'Conduite de raccordement'"));

    vector<GeoCoord> points(readDataset.extract_all_points());
    cout << "Number of points: " << points.size() << endl;

//    // Initialize all the diameter fields             --> moved inside constructor
//    string diameter_field = "DIAM_NOMIN";
//
//    for (auto track: distrib_tracks)
//        track->init_diameter(diameter_field);
//    for (auto track: raccord_tracks)
//        track->init_diameter(diameter_field);

    cout << endl;
    cout << "Before merging, 'Distribution' -> tracks: " << distrib_tracks.size() << endl;
    cout << "Before merging, 'Raccordement' -> tracks: " << raccord_tracks.size() << endl;
    cout << endl;

    //total length of the network
    double distrib_total_length_before(0), slen_before=0;
    for (auto track: distrib_tracks) {
        distrib_total_length_before += track->length;
        for (auto & sl : track->step_lengths)
            slen_before += sl;
    }

    double raccord_total_length_before(0);
    for (auto track: raccord_tracks) {
        raccord_total_length_before += track->length;
    }

    // Check that for each track, the start is equal to the first step and the end is equal to the last step
    check_tracks(distrib_tracks);
    check_tracks(raccord_tracks);

    merge_all(distrib_tracks);
    merge_all(raccord_tracks);

//    for (int i = 0; i < distrib_tracks.size(); i++) {
//        string name = "../network_data/Created_graph/Track" + to_string(i) + ".geojson";
//        distrib_tracks[i]->write_geojson(name.c_str());
//    }


    cout << "After merging, 'Distribution'  -> tracks: " << distrib_tracks.size() << "  -> steps: ";
    for (auto track: distrib_tracks)
        cout << track->steps.size() << " ";
    cout << endl;

    cout << "After merging, 'Raccordement'  -> tracks: " << raccord_tracks.size() << "  -> steps: ";
    for (auto track: raccord_tracks)
        cout << track->steps.size() << " ";
    cout << endl << endl;

    // total length of the network
    double distrib_total_length_after = 0, slen_after = 0;
    for (auto track: distrib_tracks) {
        distrib_total_length_after += track->length;
        for (auto & sl : track->step_lengths)
            slen_after += sl;
    }

    double raccord_total_length_after = 0;
    for (auto track: raccord_tracks)
        raccord_total_length_after += track->length;

    if (abs(distrib_total_length_before - distrib_total_length_after) > 1e-6) {
        cerr << "Warning: The merging process might have lost data: error between both length = " << abs(distrib_total_length_before - distrib_total_length_after) << endl;
    }

    if (abs(raccord_total_length_before - raccord_total_length_after) > 1e-6) {
        cerr << "Warning: The merging process might have lost data: error between both length = " << abs(raccord_total_length_before - raccord_total_length_after)  << endl;
    }

    // extract all the unique FID in the dataset
    // std::set<double> FIDs_distrib;
    // for (auto track: distrib_tracks) {
    //     FIDs_distrib.insert(track->getFeature()->GetFieldAsDouble("FID"));
    // }
    // std::set<double> FIDs_raccord;
    // for (auto track: raccord_tracks) {
    //     FIDs_raccord.insert(track->getFeature()->GetFieldAsDouble("FID"));
    // }
    // cout << "Number of unique FID (distrib + raccord) : " << FIDs_distrib.size() << " + " << FIDs_raccord.size() <<
    //         endl;
    // cout << endl;

    // Test the methods for tee detection
    cout << "Distrib size: " << distrib_tracks.size() << endl;

    for (auto track: distrib_tracks) {
        cout << "Diameter before segmenting: " <<
             track->diameter;
        cout << endl;
        cout << "Length before segmenting: " <<
             track->length;
        cout << endl;
    }

    segment_tracks(distrib_tracks, 0.1);                        // Segment tracks with themselves --> get branches
    segment_tracks(distrib_tracks, raccord_tracks, 0.1);     // Segment tracks with raccordements --> get segments
    //segment_tracks(raccord_tracks, 0.1);

    for (auto track: distrib_tracks) {
        cout << "Diameter after segmenting: " <<
             track->diameter;
        cout << endl;
        cout << "Length after segmenting: " <<
             track->length;
        cout << endl;
    }

    for (auto track: distrib_tracks)
        track->calculate_relative_angles(5);

    for (auto track: raccord_tracks)
        track->calculate_relative_angles(5);

    cout << "After rearranging, 'Distribution'  -> tracks: " << distrib_tracks.size() << "  -> steps: ";
    for (auto track: distrib_tracks)
        cout << track->steps.size() << " ";
    cout << endl;

    cout << "After rearranging, 'Raccordement'  -> tracks: " << raccord_tracks.size() << "  -> steps: ";
    for (auto track: raccord_tracks)
        cout << track->steps.size() << " ";
    cout << endl << endl;

    for (auto &track : raccord_tracks)
        dynamic_cast<GeoRaccordTrack *>(track)->get_sst_id();

    QGraph graph1(distrib_tracks);
    //graph1.write_geojson(R"(network_data/GPKG/debug.geojson)", GeoCoord::CH1903plus);

    //graph1.write_geojson(R"(../network_data/Created_graph/distrib_graph_new.geojson)", GeoCoord::CH1903plus);
    //QGraph graph2(raccord_tracks);
    //graph2.write_geojson(R"(../network_data/Created_graph/raccord_graph_new.geojson)", GeoCoord::CH1903plus);



//    // append the raccord_tracks to the distrib_tracks in a new all_tracks vector
//    std::vector<GeoReadTrack*> all_tracks(distrib_tracks);
//    // std::vector<GeoReadTrack*> all_tracks(raccord_tracks);
//
//    all_tracks.insert(all_tracks.end(), raccord_tracks.begin(), raccord_tracks.end());
//
//    check_tracks(all_tracks);
//
//    // add the end of the GeoTrack to the vector of GeoCoord
//    GeoGraph graph(all_tracks);
//    graph.display();
//
//    const string graph_path = "../network_data/Created_graph/created_graph.geojson";
//    graph.write_geojson(graph_path.c_str(), GeoCoord::CH1903plus);

    // *************************************************
    // Extract subgraphs
    //
//     vector<GeoGraph> subgraphs = graph.extract_subgraphs();
//
//     int index = 1;
//     for (auto& subgraph : subgraphs) {
//         std::cout << "Subgraph " << index++ << ":\n";
//         for (const auto& node : subgraph.nodes) {
//             std::cout << "  Node ID: " << node->id << "\n";
//         }
//         for (const auto& edge : subgraph.edges) {
//             std::cout << "  Edge: " << edge->nStart->id << " -> " << edge->nEnd->id << "\n";
//         }
//         std::ostringstream filename;
//         filename << "../network_data/Created_graph/subgraph_" << index++ << ".geojson";
//
//         // Write the subgraph to a GeoJSON file
//         subgraph.write_geojson(filename.str().c_str(), GeoCoord::CH1903plus);
//     }


    // new code, cornelia.blanke
    graph1.identify_source_and_sst(raccord_tracks);
    graph1.define_mst(graph1.edges);
    graph1.compute_quantum_numbers(true);
    graph1.display_mst();
    graph1.write_QNet_network_file(cout);
    graph1.write_QNet_length_file(cout);
    graph1.write_QNet_diam_file(cout);
    //graph1.write_QNet_files("C:/Users/cornelia.blanke/Bureau", "Preparation");






//    // test GeoWriteDataset
//    const char* InputFilePath = "network_data/GPKG/InputData.gpkg";
//    GeoWriteDataset inputDataset(InputFilePath, true);
//    inputDataset.write_csv_as_layer("Cases/Esta_NEW/Esta_NEW_HP.csv", "HP Catalogue");
//    inputDataset.write_csv_as_layer("Cases/Esta_NEW/Esta_NEW_pipes.csv", "Pipe Catalogue", "ID", 2);
//    inputDataset.write_csv_as_layer("Cases/Esta_NEW/Esta_NEW_meteo.csv", "Meteo Data", "Timestep");
//    inputDataset.write_csv_as_layer("Cases/Esta_NEW/Esta_NEW_SST_demand.csv", "Demand Data", "Timestep");
//    inputDataset.write_csv_as_layer("Cases/Esta_NEW/Esta_NEW_SST_dT.csv", "dT Data", "Timestep");



    const char* QnetFilePath = "network_data/GPKG/QuantumizedNetwork.gpkg";
    GeoWriteDataset writeDataset(QnetFilePath, true);
    writeDataset.createNodeLayer();
    writeDataset.createEdgeLayer();
    writeDataset.createSSTLayer();
    writeDataset.write_nodes({graph1.source});
    writeDataset.write_nodes(graph1.nodes);
    writeDataset.write_edges(graph1.mst_edges);
    writeDataset.copy_file_to_sstLayer("network_data/GPKG/test2_SST.csv");
//    writeDataset.copy_geometry_to_sst();
    writeDataset.write_tracks_to_sstLayer(raccord_tracks);

    QNetwork net;
    net.read_tables(QnetFilePath);
    net.display();
    net.display_table(net.diamtable);

    return 0;

}