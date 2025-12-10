//
// Created by cornelia.blanke on 10.07.2024.
//

#ifndef FLUIDS_QGRAPH_H
#define FLUIDS_QGRAPH_H

#include<set>
#include<unordered_map>
#include "GeoGraph.h"


/*** --------------------------------------- class QNode ---------------------------------------- ***/
class QEdge;
class QNode : public Node {
public:
    enum QNodeType {CENTRAL=0, PROTON=1, NEUTRON=2};
    QNodeType type=QNodeType::PROTON;
    // CENTRAL does not use quantum numbers
    int l=0;            // for NEUTRON: quantum number of branch above
    unsigned int k=0;   // for NEUTRON: store n in k (use a reference)
    unsigned int& n = k;   // (alias only for readability)
    unsigned int sst_id = 0;     // only used for PROTON

public:
    QNode() : Node(), type(QNodeType::PROTON) {};
    explicit QNode(GeoCoord* coord, QNodeType type=QNodeType::PROTON) : Node(coord), type(type) {};
    explicit QNode(GeoCoord& coord, QNodeType type=QNodeType::PROTON) : Node(coord), type(type) {};
    explicit QNode(Node* node) : Node(*node) {};    // does not increase Node::nodecounter
    explicit QNode(Node& node) : Node(node) {};

    QNode(const QNode& node)  : Node(node) {
        type = node.type;
        l = node.l;
        k = node.k;
        n = k;
    }
    QNode& operator= (const QNode& node){
        id = node.id;
        coord = node.coord;
        type = node.type;
        l = node.l;
        k = node.k;
        n = k;
        return *this;
    }
    //~QNode() override =default;
    void set_type(QNodeType itype=QNodeType::PROTON, int il=0, unsigned int ik=0) { type=itype; l=il; k=ik; }
    void display_more() const override;
};

/*** --------------------------------------- class QEdge ---------------------------------------- ***/
class QEdge : public Edge<QNode> {
public:
    int l=0;
    unsigned int k=0;
    double length=-1;
    unsigned int DN=0;     // nominal diameter

public:
    QEdge() : Edge() {}
    QEdge(QNode* A, QNode* B, Edge::TrackType type=Edge::BASE) {
        nStart = A;
        nEnd = B;
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
        init_length();
    }
    ~QEdge() override = default;

    void init_length(double len=-1);
    void init_DN(unsigned int diam=0);

};

/*** --------------------------------------- class QGraph --------------------------------------- ***/
class QGraph : public Graph<QNode, QEdge> {
public:
    QNode* source=nullptr;
    vector<QEdge*> mst_edges;

private:
    bool source_created=false;


public:
    QGraph()=default;
    ~QGraph();

    // from a source and a vector of coords, compute tracks between all pairs
    explicit QGraph(GeoCoord& sourcecoord, vector<GeoCoord>& coords, double tol=0.1, QEdge::TrackType type=QEdge::BASE);
    explicit QGraph(GeoCoord* sourcecoord, const vector<GeoCoord*>& coords, double tol=0.1, QEdge::TrackType type=QEdge::BASE);
    QGraph(GeoCoord& sourcecoord, vector<GeoCoord>& coords, QEdge::TrackType type, double tol=0.1);
    QGraph(GeoCoord* sourcecoord, const vector<GeoCoord*>& coords, QEdge::TrackType type, double tol=0.1);

    // from a Graph, transform Node/Edge into QNode/QEdge, add a source and connect
    QGraph(GeoCoord& sourcecoord, GeoGraph& graph, QEdge::TrackType type=QEdge::BASE);

    // from a vector of (unconnected) tracks
    explicit QGraph(const vector<GeoTrack*>& tracks, double tol=0.1);

    // Overloaded constructor to create a GeoGraph from a vector of GeoReadTrack pointers
    explicit QGraph(const vector<GeoReadTrack *> &tracks, double tol=0.1);
    explicit QGraph(QNode* sourcenode, const vector<GeoReadTrack *> &tracks, double tol=0.1);


    void create_source(GeoCoord& coord);
    void create_source(Node* node);
    void create_source(QNode* node);

    void redefine_node_as_source(QNode* node);

    void identify_source_and_sst(const vector<GeoReadTrack*>& raccord_tracks, double tol= 0.1);     // identifies pipe with an open end --> source

    using Graph<QNode, QEdge>::add;
    void add(vector<QEdge*>& more_edges) override;

    vector<QNode*> get_allnodes() override;

    void display() const override;
    void display_mst() const;

    void define_mst(const vector<QEdge*>& edgevector);

    void find_mst_without_neutrons();   // all nodes are treated in the same way (protons/neutrons are computed afterwards)
    void find_mst_with_neutrons();      // distinguishes between protons and neutrons
    void find_mst_with_neutrons2();     // NOT WORKING: reverse-delete algorithm, distinguishes between protons and neutrons

    void write_mst_geojson(const char *fileout, GeoCoord::CRS crs=GeoCoord::CH1903plus);

    void compute_quantum_numbers(bool check_for_errors=true);
    void compute_quantum_numbers(QNode* node, set<QEdge*>& candidate_edges, unsigned int current_n, int current_l, unsigned int current_k);

    void write_QNet_files(const char* case_dir, const char* case_name);
    void write_QNet_network_file(ostream& out);
    void write_QNet_length_file(ostream& out);
    void write_QNet_diam_file(ostream& out);
};


#endif //FLUIDS_QGRAPH_H
