//
// Created by cornelia.blanke on 25.06.2024.
//

#ifndef FLUIDS_GEOGRAPH_H
#define FLUIDS_GEOGRAPH_H

#include <unordered_map>
#include <map>
#include"GeoTrack.h"
#include"GeoFactorTrack.h"
#include"GeoRoutingTrack.h"
#include"GeoReadTrack.h"


// discussion about circular dependency in template classes:
// https://cplusplus.com/forum/general/115940/

/*** --------------------------------------- class Node ----------------------------------------- ***/
class Node {
public:
    static int nodecounter;     // automatic numbering of all nodes
    unsigned int id;
    GeoCoord* coord=nullptr;

public:
    Node() : id(nodecounter++) {};
    explicit Node(GeoCoord* coord) : id(nodecounter++), coord(coord) {};
    explicit Node(GeoCoord& coord) : id(nodecounter++), coord(&coord) {};
    virtual ~Node()=default;
    bool operator==(const Node& node) const {
        return (this->id == node.id);
    }
    bool operator<(const Node& node) const {
        return (this->id < node.id);
    }
    void display(GeoCoord::CRS crs=GeoCoord::CH1903plus) const;
    virtual void display_more() const {}    // overridden in inherited class
};


/*** --------------------------------------- class Edge ----------------------------------------- ***/
template <typename N>
class Edge {
public:
    N* nStart=nullptr, *nEnd=nullptr;
    GeoTrack* track=nullptr;
    enum TrackType {NONE=0, BASE=1, FACTOR=2, ROUTING=3, READ=4};

protected:
    bool track_created=false;

public:
    Edge()=default;
    Edge(N* A, N* B, Edge::TrackType type=Edge::BASE) : nStart(A), nEnd(B) {
        if (*A == *B){
            cerr << "Two identical Nodes are not allowed in an Edge" << endl;   // no circular edge
            exit(1);
        }
        if (type==Edge::BASE) {
            track = new GeoTrack(nStart->coord, nEnd->coord);
            track_created  =true;
        }
        else if (type==Edge::FACTOR) {
            track = new GeoFactorTrack(nStart->coord, nEnd->coord);
            track_created  =true;
        }
        else if (type==Edge::ROUTING) {
            track = new GeoRoutingTrack(nStart->coord, nEnd->coord);
            track_created  =true;
        }
    }
    virtual ~Edge() {if (track_created) delete track; };

    bool operator==(const Edge& edge) const {   // undirected edge
        return ((this->nStart == edge.nStart || this->nStart == edge.nEnd) &&
            (this->nEnd == edge.nStart || this->nEnd == edge.nEnd));
    }
    bool operator<(const Edge& edge) const {
        unsigned int min1 = min(this->nStart->id, this->nEnd->id);
        unsigned int max1 = max(this->nStart->id, this->nEnd->id);
        unsigned int min2 = min(edge.nStart->id, edge.nEnd->id);
        unsigned int max2 = max(edge.nStart->id, edge.nEnd->id);
        return (min1 < min2 || (min1 == min2 && max1 < max2));
    }
    // https://stackoverflow.com/questions/4660123/overloading-friend-operator-for-class-template/
    template <typename T>
    friend ostream& operator<<(ostream& out, const Edge<T>& edge);
    Edge::TrackType get_tracktype() const;
    bool contains(N* node);
    N* get_othernode(const N* node);
    void distance(const N* node, double& distance, double& angle_a, double& angle_b, bool from_track=false, double tol=1e-16);
    void display() const;
};

/*** --------------------------------- template class Graph -------------------------------------- ***/
template <typename N, typename E>
class Graph {
public:
    vector<N*> nodes;
    vector<E*> edges;
    unordered_map<N*, vector<E*>> connected_edges;

protected:
    // for destructor:
    vector<N*> nodes_to_delete;
    vector<E*> edges_to_delete;
    // very big dummy distance for unconnected nodes
    double NOT_CONNECTED=1.0e16;


public:
    Graph()=default;
    ~Graph();

    // from a vector of coords, create distinct nodes and edges
    explicit Graph(vector<GeoCoord>& coords, double tol=0.1, typename E::TrackType type=E::BASE);
    explicit Graph(const vector<GeoCoord*>& coords, double tol=0.1, typename E::TrackType type=E::BASE);
    Graph(vector<GeoCoord>& coords, typename E::TrackType type, double tol=0.1);
    Graph(const vector<GeoCoord*>& coords, typename E::TrackType type, double tol=0.1);

    // from nodes, create edges
    explicit Graph(const vector<N*>& nodevector, typename E::TrackType type=E::BASE);

    // from edges, extract nodes
    explicit Graph(const vector<E*>& edgevector);

    // from a vector of (unconnected) tracks
    explicit Graph(const vector<GeoTrack*>& tracks, double tol=0.1);
    explicit Graph(const vector<GeoReadTrack *> &tracks, double tol=0.1);

    void build_graph(const vector<GeoTrack*>& tracks, double tol=0.1);

    void add(vector<N*>& more_nodes);
    virtual void add(vector<E*>& more_edges);
    void add(E* edge);

    // if nodes removed -> connected edges are removed
    void remove(vector<N*>& less_nodes);
    void remove(vector<N>& less_nodes);
    void remove(vector<int>& less_nodeids);
    void remove(vector<unsigned int>& less_nodeids);
    void remove(N* node);
    void remove(N& node);

    void remove(vector<E*>& less_edges);
    void remove(vector<E>& less_edges);
    void remove(N* A, N* B);       // remove edge from A to B
    void remove(N& A, N& B);
    void remove(unsigned int A, unsigned int B);

    E* get_edge(N* A, N* B);
    E* get_edge(N& A, N& B);

    virtual vector<N*> get_allnodes();          // get nodes, overridden for QGraph (source+nodes)

    void display_nodeinfo(N* node);
    void display_nodeinfo(N& node);
    void display_nodeinfo(int id);
    virtual void display() const;

    void check_connected_edges();
    void compute_connected_edges();     // delete current connected_edges and re-compute
    bool check_connectivity();          // use connected_edges to determine connectivity of all edges

    int get_nb_of_clusters();           //get number of unconnected clusters


    static vector<GeoCoord*> suppress_duplicate_coords(vector<GeoCoord>& coords, double tol= 0.1);
    static vector<GeoCoord*> suppress_duplicate_coords(const vector<GeoCoord*>& coords, double tol= 0.1);

    void create_nodes_from_coords(const vector<GeoCoord*>& coords);
    void create_edges_from_nodes(typename E::TrackType type=E::BASE);

    void merge_nodes(double tol= 0.1);

    // create edges between all combinations of nodevector
    void connect_nodes(vector<N*>& nodevector, typename E::TrackType type=E::BASE);
    // create edges between all elements of nodevector1 and node
    void connect_nodes(vector<N*>& nodevector, N* node, typename E::TrackType type=E::BASE);
    // create edges between all combinations of nodevector1 and nodevector2
    void connect_nodes(vector<N*>& nodevector1, vector<N*>& nodevector2, typename E::TrackType type=E::BASE);


    // from curl "https://routing.openstreetmap.de..."
    // computes the full distance matrix, but does not output the intermediate steps for each track
    void connect_nodes_with_routing(const vector<N*>& nodevector);
    void connect_nodes_with_routing(vector<N*>& nodevector, N* node);
    void connect_nodes_with_routing(vector<N*>& nodevector1, vector<N*>& nodevector2);

    // Subgraph methods
    vector<vector<N*>> find_connected_subgraphs();
    vector<Graph<N, E>> extract_subgraphs();

    vector<E*> find_shortest_edges(N* node, unsigned int n_lowest=1);

    map<pair<N*, N*>, double> get_adjacencyMatrix(const vector<N*>& nodevector, unordered_map<N*, vector<E*>>& map_connected_edges);
    map<pair<N*, N*>, double> get_adjacencyMatrix();
    map<pair<N*, N*>, double> get_adjacencyMatrix_directed(const vector<N*>& nodevector, const vector<E*>& edgevector);  // directed (nStart, nEnd)

    // https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
    void FloydWarshall_directed(map<pair<N*, N*>, double>& distance_matrix, map<pair<N*, N*>, N*>& next_node, bool compute_next=true);
    void FloydWarshall_undirected(map<pair<N*, N*>, double>& distance_matrix, map<pair<N*, N*>, N*>& next_node, bool compute_next=true);
    map<pair<N*, N*>, double> FloydWarshall_dir();        // directed graph
    map<pair<N*, N*>, double> FloydWarshall_undir();      // undirected graph
    vector<N*> reconstruct_path(map<pair<N*, N*>, N*>& next_node, N* A, N* B);         // get shortest path from A to B
    vector<E*> get_edges_on_path(vector<N*> path, unordered_map<N*, vector<E*>>& map_connected_edges);
    vector<E*> get_edges_on_path(vector<N*> path);

    void write_geojson(const char *fileout, GeoCoord::CRS crs=GeoCoord::CH1903plus);

};

/*** -------------------------------------- GeoGraph --------------------------------------------- ***/
using GeoGraph = Graph<Node, Edge<Node>>;


#endif //FLUIDS_GEOGRAPH_H
