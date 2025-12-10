//
// Created by cornelia.blanke on 25.06.2024.
//

#include <random>
#include "geodata/Pixelation.h"
#include "geodata/GeoGraph.h"
#include "geodata/QGraph.h"


using namespace std;

int main() {
    // Constructors
    GeoCoord coord(GeoCoord::WGS84, 7.15800, 46.79402, 500);
    coord.display();

    GeoCoord coord1(GeoCoord::WGS84, 7.13840, 46.79448, 500);
    cout << distance_2d(coord, coord1) << "\n";
    cout << distance_3d(coord, coord1) << "\n";

    GeoTrack track(coord, coord1);
    track.display();
    //track.write_geojson(R"(C:\Users\cornelia.blanke\Bureau\mytest.geojson)");
    GeoFactorTrack ftrack(coord, coord1, 2);
    ftrack.display();
    GeoRoutingTrack rtrack(coord, coord1);
    rtrack.compute_track();
    rtrack.display();
//    rtrack.write_geojson(R"(C:\Users\cornelia.blanke\Bureau\mytest2.geojson)");

    cout << "\n\n";
    vector<GeoCoord> vect={coord, coord1};
    Pixelation pixelation(vect, 10);
    pixelation.display(GeoCoord::CH1903, true);
    GeoCoord coord2(GeoCoord::WGS84, 7.14217, 46.79213, 500);
    cout << pixelation.find_closest_to(coord2)[0]->as_string(GeoCoord::CH1903) << "\n";

    Node node1(coord1);
    QNode qnode1(node1);
    node1.display();
    Node node2(coord2);
    node2.display();
    Edge<Node> edge12(&node1, &node2);
    edge12.display();

    GeoCoord coord3(GeoCoord::WGS84, 7.142816, 46.797845, 500);
    GeoCoord coord4(GeoCoord::WGS84, 7.15269, 46.79955, 500);
    vector<GeoCoord> coords={coord1, coord2, coord3, coord4};

    GeoGraph mygraph(coords, Edge<Node>::BASE);
    mygraph.display();
    //mygraph.write_geojson(R"(C:\Users\cornelia.blanke\Bureau\mytest3.geojson)");

    QGraph myQgraph(coord, mygraph, QEdge::BASE);
    myQgraph.display();
    myQgraph.display_nodeinfo(4);
    myQgraph.find_mst_without_neutrons();
    myQgraph.display_mst();
    myQgraph.compute_quantum_numbers();
    myQgraph.write_QNet_network_file(cout);
    cout << "\n";
    myQgraph.write_QNet_length_file(cout);
    //myQgraph.write_QNet_files("C:/Users/cornelia.blanke/Bureau", "Test");






    // create and test a more complex random graph
    cout << "\nCheck testcase:\n";
    vector<GeoCoord> random_coords;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> x(7.12, 7.17);
    std::uniform_real_distribution<double> y(46.785, 46.815);

    GeoCoord source_coord = GeoCoord(GeoCoord::WGS84, x(mt), y(mt), 500);
    for (int i=0; i<20; ++i)
        random_coords.emplace_back(GeoCoord::WGS84, x(mt), y(mt), 500);

    QGraph test(source_coord, random_coords);
    test.find_mst_without_neutrons();
    test.display_mst();
    //test.write_geojson(R"(C:\Users\cornelia.blanke\Bureau\randomtest.geojson)");
    //test.write_mst_geojson(R"(C:\Users\cornelia.blanke\Bureau\randommst.geojson)");
    test.compute_quantum_numbers();
    test.write_QNet_network_file(cout);
    cout << endl;
    test.write_QNet_length_file(cout);

    cout << "Undirected: Test Adjacency and Floyd-warshall" << endl;
    Node n1(coord1);
    Node n2(coord2);
    Node n3(coord3);
    Node n4(coord4);
    GeoGraph floydgraph({&n1, &n2, &n3, &n4});
    floydgraph.remove(&n1, &n4);
    floydgraph.get_edge(&n1, &n3)->track->length = 2;
    floydgraph.get_edge(&n3, &n4)->track->length = 2;
    Edge<Node>* n4n2 = floydgraph.get_edge(&n4, &n2);
    n4n2->track->length = 1;
    n4n2->nStart = &n4;
    n4n2->nEnd = &n2;
    floydgraph.get_edge(&n2, &n3)->track->length = 3;
    Edge<Node>* n2n1 = floydgraph.get_edge(&n2, &n1);
    n2n1->track->length = 4;
    n2n1->nStart = &n2;
    n2n1->nEnd = &n1;



    floydgraph.compute_connected_edges();
    map<pair<Node*, Node*>, double> adj = floydgraph.get_adjacencyMatrix({&n1, &n2, &n3, &n4}, floydgraph.connected_edges);
    for (auto & it : adj)
        cout << it.first.first->id << "\t" << it.first.second->id << "\t" << it.second << "\n";
    cout << "\n";
    map<pair<Node*, Node*>, double> shortest;
    map<pair<Node*, Node*>, Node*> next_node;
    floydgraph.FloydWarshall_undirected(shortest, next_node, true);
    for (auto & it : shortest)
        cout << it.first.first->id << "\t" << it.first.second->id << "\t" << it.second << "\n";
    cout << "\n";
    for (auto & it : next_node)
        cout << it.first.first->id << "\t" << it.first.second->id << "\t" << it.second->id << "\n";

    cout << "Directed: Test Adjacency and Floyd-warshall" << endl;
    floydgraph.get_edge(&n1, &n3)->track->length = -2;
    n4n2->track->length = -1;
    adj = floydgraph.get_adjacencyMatrix_directed({&n1, &n2, &n3, &n4},{floydgraph.get_edge(&n1, &n3),
                                                                              floydgraph.get_edge(&n3, &n4),
                                                                              floydgraph.get_edge(&n4, &n2),
                                                                              floydgraph.get_edge(&n2, &n3),
                                                                              floydgraph.get_edge(&n2, &n1)});
    for (auto & it : adj)
        cout << it.first.first->id << "\t" << it.first.second->id << "\t" << it.second << "\n";
    cout << "\n";

    floydgraph.FloydWarshall_directed(shortest, next_node, true);
    for (auto & it : shortest)
        cout << it.first.first->id << "\t" << it.first.second->id << "\t" << it.second << "\n";
    cout << "\n";
    for (auto & it : next_node)
        cout << it.first.first->id << "\t" << it.first.second->id << "\t" << it.second->id << "\n";

    cout << "Reconstruct" << endl;
    vector<Node*> path = floydgraph.reconstruct_path(next_node, floydgraph.nodes[1], floydgraph.nodes[3]);
    for (auto &node : path)
        cout << node->id << "\t";
    cout << "\n\n";

    cout << "New test" << "\n";
    test.display();
    for (auto &node : test.nodes) {
        node->set_type();       // reset all
        if (node->id == 22 || node->id == 23){
            node->set_type(QNode::QNodeType::NEUTRON);
            test.remove(test.source, node);
        }
        //test.display_nodeinfo(node);
    }
    test.find_mst_with_neutrons();
    test.display_mst();
    test.compute_quantum_numbers();
    test.write_QNet_network_file(cout);
    cout << endl;
    test.write_mst_geojson(R"(C:\Users\cornelia.blanke\Bureau\newMST.geojson)");


    return 0;
}