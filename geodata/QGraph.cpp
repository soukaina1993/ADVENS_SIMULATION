//
// Created by cornelia.blanke on 10.07.2024.
//

#include<algorithm>
#include <fstream>
#include<map>
#include "QGraph.h"

void QNode::display_more() const {
    cout << "Quantum Type " ;
    if (type == QNodeType::CENTRAL) {
        cout << "CENTRAL";
    }
    else if (type == QNodeType::PROTON) {
        cout << "PROTON";
        cout << " l=" << l << ", k=" << k;
    }
    else {
        cout << "NEUTRON";
        cout << " n=" << n;
    }
    cout << endl;
}

void QEdge::init_length(const double len) {
    if (len >= 0)
        length = len;
    else if (track != nullptr)
        length = track->length;
    else {
        cerr << "Please specify a length" << endl;
        exit(1);
    }
}

void QEdge::init_DN(const unsigned int diam) {
    if (diam > 0)
        DN = diam;
    else if (dynamic_cast<GeoReadTrack*>(track))
        DN = dynamic_cast<GeoReadTrack*>(track)->diameter;
    else {
        cerr << "Please specify a nominal diameter" << endl;
        exit(1);
    }
}



QGraph::QGraph(GeoCoord& sourcecoord, vector<GeoCoord>& coords, double tol, QEdge::TrackType type) {
    vector<GeoCoord*> singlecoords = suppress_duplicate_coords(coords, tol);
    create_source(sourcecoord);
    create_nodes_from_coords(singlecoords);

    // connect nodes by edges
    if (type == QEdge::BASE) {
        connect_nodes(nodes, source, type);
        connect_nodes(nodes, type);
    }
    else if (type == QEdge::FACTOR) {
        connect_nodes(nodes, source, type);
        connect_nodes(nodes, type);
    }
    else if (type == QEdge::TrackType::ROUTING) {
        connect_nodes_with_routing(nodes, source);
        connect_nodes_with_routing(nodes);
    }
}

QGraph::QGraph(GeoCoord* sourcecoord, const vector<GeoCoord*>& coords, double tol, QEdge::TrackType type){
    vector<GeoCoord*> singlecoords = suppress_duplicate_coords(coords, tol);
    create_source(*sourcecoord);
    create_nodes_from_coords(singlecoords);

    // connect nodes by edges
    if (type == QEdge::BASE) {
        connect_nodes(nodes, source, type);
        connect_nodes(nodes, type);
    }
    else if (type == QEdge::FACTOR) {
        connect_nodes(nodes, source, type);
        connect_nodes(nodes, type);
    }
    else if (type == QEdge::TrackType::ROUTING) {
        connect_nodes_with_routing(nodes, source);
        connect_nodes_with_routing(nodes);
    }
}

QGraph::QGraph(GeoCoord& sourcecoord, vector<GeoCoord>& coords, QEdge::TrackType type, double tol) :
        QGraph(sourcecoord, coords, tol, type) {}
QGraph::QGraph(GeoCoord* sourcecoord, const vector<GeoCoord*>& coords, QEdge::TrackType type, double tol) :
        QGraph(sourcecoord, coords, tol, type){}

QGraph::QGraph(GeoCoord& sourcecoord, GeoGraph& graph, QEdge::TrackType type){
    map<unsigned int, QNode*> get_QNode;
    map<unsigned int, QNode*>::iterator it;

    // convert edges
    for (auto &edge : graph.edges){
        auto* new_edge = new QEdge();
        new_edge->track = edge->track;
        edges.push_back(new_edge);
        edges_to_delete.push_back(new_edge);

        // convert nodes of edge
        for (auto &node : {edge->nStart, edge->nEnd}){
            it = get_QNode.find(node->id);
            if (it != get_QNode.end()) {         // node already in nodes, define edge
                (node == edge->nStart) ? (new_edge->nStart = it->second) : (new_edge->nEnd = it->second);
                connected_edges[it->second].push_back(new_edge);
            }
            else {                                          // new QNode, define edge
                auto* new_node = new QNode(node);
                get_QNode[node->id] = new_node;
                nodes.push_back(new_node);
                nodes_to_delete.push_back(new_node);
                (node == edge->nStart) ? (new_edge->nStart = new_node) : (new_edge->nEnd = new_node);
                connected_edges[new_node].push_back(new_edge);
            }
        }
    }

    // convert standalone nodes
    for (auto &node : graph.nodes){
        if (get_QNode.count(node->id) == 0){    // node not yet converted (i.e. not connected to any edge)
            auto* new_node = new QNode(node);   // new QNode
            get_QNode[node->id] = new_node;
            nodes.push_back(new_node);
            nodes_to_delete.push_back(new_node);
            connected_edges[new_node] = {};
        }
    }

    // add source
    create_source(sourcecoord);
    // connect source with rest of graph
    if (type == QEdge::BASE)
        connect_nodes(nodes, source, type);
    else if (type == QEdge::FACTOR)
        connect_nodes(nodes, source, type);     //TODO track->compute_penalty()
    else if (type == QEdge::ROUTING)
        connect_nodes_with_routing(nodes, source);
}

QGraph::QGraph(const vector<GeoTrack*>& tracks, const double tol) : Graph<QNode, QEdge>(tracks, tol) {
    for (auto &edge : edges)
        edge->init_length();
}

QGraph::QGraph(const vector<GeoReadTrack *> &tracks, const double tol) : Graph<QNode, QEdge>(tracks, tol) {
    for (auto &edge : edges) {
        edge->init_length();
        edge->init_DN();
    }
}

QGraph::QGraph(QNode* sourcenode, const vector<GeoReadTrack *> &tracks, const double tol) : QGraph(tracks, tol){
    source = sourcenode;
    sourcenode->type = QNode::QNodeType::CENTRAL;
    for (auto &edge : edges) {
        edge->init_length();
        edge->init_DN();
    }
}

QGraph::~QGraph(){
    if (source_created)
        delete source;
}





void QGraph::create_source(GeoCoord& coord) {
    source_created=true;
    source = new QNode(coord, QNode::QNodeType::CENTRAL);
    connected_edges[source] = {};
}
void QGraph::create_source(Node* node){
    source_created=true;
    source = new QNode(node);
    connected_edges[source] = {};
}

void QGraph::create_source(QNode* node){
    source = node;
    node->type = QNode::QNodeType::CENTRAL;
    connected_edges[source] = {};
}

void QGraph::redefine_node_as_source(QNode* node){
    vector<QNode*>::iterator it;
    it = find(nodes.begin(), nodes.end(), node);

    // remove node from nodes
    if (it == nodes.end()){
        cerr << "Node " << node->id << " is not found on nodes vector" << endl;
        exit(1);
    }
    else
        nodes.erase(it);

    // define node as source
    source = node;
    source->type = QNode::QNodeType::CENTRAL;
}

void QGraph::identify_source_and_sst(const vector<GeoReadTrack*>& raccord_tracks, const double tol) {
    double tol2 = pow(tol, 2);
    set<GeoReadTrack*> candidate_tracks (raccord_tracks.begin(), raccord_tracks.end());
    bool found;
    const vector <QNode*> all_nodes = nodes;     // do not loop over nodes as it will be modified
    for (auto &node: all_nodes) {
        found = false;
        for (auto raccord: candidate_tracks) {
            if (sqr_distance_2d(*raccord->start, *node->coord) < tol2 || sqr_distance_2d(*raccord->end, *node->coord) < tol2) {
                node->type = QNode::PROTON;
                candidate_tracks.erase(raccord);     // raccordement belongs to node
                if (auto r=dynamic_cast<GeoRaccordTrack*>(raccord))
                    node->sst_id = r->sst_id;
                found = true;
                break;
            }
        }
        if (!found) {
            if (connected_edges[node].size() == 1)
                redefine_node_as_source(node);
            else if (connected_edges[node].size() >= 3)
                node->type = QNode::NEUTRON;
            else {
                cerr << "Something wrong: Node " << node->id << " is neither Central Plant, Proton, nor Neutron" << endl;
                exit(1);
            }
        }
    }
    if (!candidate_tracks.empty()){
        cerr << candidate_tracks.size() << " raccord_tracks are not connected to network" << endl;
    }
}


void QGraph::add(vector<QEdge*>& more_edges){
    // add nodes if not yet present
    for (auto &edge : more_edges){
        for (auto &node : {edge->nStart, edge->nEnd}) {
            //create node if does not exist
            if (node->type == QNode::QNodeType::CENTRAL) {
                if (source == nullptr || source == node)
                    source = node;
                else {
                    cerr << "Do not override source" << endl;
                    exit(1);
                }
            } else    // PROTON or NEUTRON
                if (find(nodes.begin(), nodes.end(), node) == nodes.end())    //not yet in nodes
                    nodes.push_back(node);

            connected_edges[node].push_back(edge);
        }
    }
    edges.insert( edges.end(), more_edges.begin(), more_edges.end());
}

vector<QNode*> QGraph::get_allnodes(){
    vector<QNode*> allnodes (1, source);
    allnodes.insert(allnodes.end(), nodes.begin(), nodes.end());
    return allnodes;
}


void QGraph::display() const {
    if (source != nullptr)
        source->display();
    Graph::display();
}

void QGraph::display_mst() const {
    cout << "\n" << mst_edges.size() << " MST Edges: ";
    for (auto &edge : mst_edges)
        cout << *edge << " ";
    cout << endl;
}

void QGraph::define_mst(const vector<QEdge*>& edgevector){
    mst_edges = edgevector;
}

void QGraph::find_mst_without_neutrons() {
    if (edges.empty()){
        cerr << "No edges defined" << endl;
        exit(1);
    }
    set<QNode*> search_nodes;
    unordered_map<QNode*, unsigned short> edge_count;
    set<QEdge*> candidate_edges (edges.begin(), edges.end());    // all edges
    QEdge* min_edge, *node_min_edge;

    // start with source, find the shortest edge
    mst_edges = find_shortest_edges(source, 1);

    for (auto &n_edge : connected_edges[source])
        candidate_edges.erase(n_edge);

    QNode* other_node = mst_edges[0]->get_othernode(source);
    //(mst_edges[0]->nStart == source) ? (other_node = mst_edges[0]->nEnd) : (other_node = mst_edges[0]->nStart);
    search_nodes.insert(other_node);
    edge_count[other_node] = 1;

    // maximum three connected edges
    while (!search_nodes.empty()) {
        min_edge = nullptr;
        for (auto &node: search_nodes) {
            node_min_edge = nullptr;
            for (auto &edge: connected_edges[node]) {
                if ((candidate_edges.count(edge) > 0) &&         // edge is a candidate
                    (node_min_edge == nullptr || (edge->length < node_min_edge->length)))
                        node_min_edge = edge;
            }
            if (node_min_edge != nullptr &&
                (min_edge == nullptr || (node_min_edge->length < min_edge->length)))
                    min_edge = node_min_edge;
        }
        if (min_edge != nullptr) {
            mst_edges.push_back(min_edge);
            candidate_edges.erase(min_edge);
            for (auto &node : {min_edge->nStart, min_edge->nEnd}) {
                if (edge_count.find(node) == edge_count.end())
                    edge_count[node] = 0;       // add to map

                if (edge_count[node] > 1) {                   // is used for the third time
                    edge_count[node]++;
                    search_nodes.erase(node);
                    for (auto &n_edge: connected_edges[node])
                        candidate_edges.erase(n_edge);
                }
                else {
                    search_nodes.insert(node);
                    edge_count[node]++;
                    for (auto &n_edge: connected_edges[node])
                        if (edge_count.find(n_edge->nStart) != edge_count.end() &&
                            edge_count.find(n_edge->nEnd) != edge_count.end())
                            candidate_edges.erase(n_edge);
                    }
                }
            }
        else
            break;      // no edge found, stop while-loop
    }
    if (edge_count.size() < nodes.size()){
        cerr << "Unconnected node(s)" << endl;
        exit(1);
    }
//    else
//        sort(mst_edges.begin(), mst_edges.end(), [](const Edge* a, const Edge* b){return (*a < *b);});
}

void QGraph::find_mst_with_neutrons() {
    if (edges.empty()){
        cerr << "No edges defined: MST cannot be computed" << endl;
        exit(1);
    }
    mst_edges.clear();
    vector<QNode*> protons, neutrons;

    for (auto &node : nodes)
        (node->type == QNode::PROTON) ? protons.push_back(node) : neutrons.push_back(node);

    map<pair<QNode*, QNode*>, double> distance_matrix;
    map<pair<QNode*, QNode*>, QNode*> next_node;

    FloydWarshall_undirected(distance_matrix, next_node, true);  // contains source+protons+neutrons

    QGraph subgraph;
    subgraph.create_source(source);
    subgraph.add(neutrons);

    //create edges from source to neutrons with length of shortest path
    for (auto neutron : neutrons){
        auto* edge = new QEdge(source, neutron);
        edge->length = distance_matrix[make_pair(source, neutron)];
        subgraph.add(edge);
        subgraph.edges_to_delete.push_back(edge);
    }
    //create edges between all pairs of neutrons
    if (!neutrons.empty()) {
        for (auto &n1: neutrons) {
            for (auto &n2: neutrons) {
                if (n1 != n2) {
                    auto *edge = new QEdge(n1, n2);
                    edge->length = distance_matrix[make_pair(n1, n2)];
                    subgraph.add(edge);
                    subgraph.edges_to_delete.push_back(edge);
                }
            }
        }
        //find MST on this subgraph
        subgraph.find_mst_without_neutrons();
    }

    vector<QNode*> pathnodes;
    vector<QEdge*> pathedges;
    map<QNode*, unsigned int> free_neutrons;    // free neutrons with their number of connected edges
    set<QNode*> free_protons(protons.begin(), protons.end());

    if (neutrons.empty()) {     // connect nearest proton
        double nearest_dist=NOT_CONNECTED;
        QNode* nearest_p;
        for (auto &proton  : protons){
            if (distance_matrix[make_pair(source, proton)] < nearest_dist){
                nearest_dist = distance_matrix[make_pair(source, proton)];
                nearest_p = proton;
            }
        }
        // add shortest edge from source
        mst_edges.push_back(get_edge(source, nearest_p));
        free_protons.erase(nearest_p);
    }

    for(auto &edge : subgraph.mst_edges){
        if (edge->nStart != source)
            free_neutrons[edge->nStart]++;      // creates entry if not yet present
        if (edge->nEnd != source)
            free_neutrons[edge->nEnd]++;
        pathnodes = reconstruct_path(next_node, edge->nStart, edge->nEnd);
        for (auto &pathnode : pathnodes)
            free_protons.erase(pathnode);
        pathedges = get_edges_on_path(pathnodes);
        mst_edges.insert(mst_edges.end(), pathedges.begin(), pathedges.end());
    }
    pathnodes.clear();
    pathedges.clear();
    //TODO clear subgraph?!


    for (auto &entry : free_neutrons){
        if (entry.second > 3){
            cerr << "More than three edges: Error in graph calculation" << endl;
            exit(1);
        }
        else if (entry.second == 3)      // not free
            free_neutrons.erase(entry.first);
    }

    // connect free neutrons to nearest free protons (nothing done if empty)
    map<double, pair<QNode*, QNode*>> swapped;      // first: distance, second: free nodes (a map is sorted by its keys)
    double dist;
    for(auto &free_n : free_neutrons)
        for(auto &free_p : free_protons) {
            dist = distance_matrix[make_pair(free_n.first, free_p)];
            if (swapped.find(dist) != swapped.end())     // key already there, increase slightly
                dist *= 1.00001;
            swapped[dist] = make_pair(free_n.first, free_p);
        }

    for (auto &swap :swapped){
        for (auto &free_n : free_neutrons){
            if (swap.second.first == free_n.first) {
                //store connection from proton to neutron in MST
                QEdge* edge = get_edge(free_n.first, swap.second.second);
                if (edge == nullptr){
                    cerr << "No edge found" << endl;
                    exit(1);
                }
                else
                    mst_edges.push_back(edge);
                free_n.second++;
                free_protons.erase(swap.second.second);
                if (free_n.second == 3) {      // completely connected
                    free_neutrons.erase(free_n.first);
                    if (free_neutrons.empty())      // last neutron removed
                        break;
                }
            }
        }
    }
    if (!free_neutrons.empty()){
        cerr << "Neutron " << (*free_neutrons.begin()).first->n << " is still not connected" << endl;
        exit(1);
    }
    else {
        //connect free_protons
        set<QNode*> candidates(neutrons.begin(), neutrons.end());
        for (auto &proton : protons){
            if (free_protons.find(proton) == free_protons.end())    // proton not free
                candidates.insert(proton);
        }

        while (!free_protons.empty()){
            QNode *free_proton, *nearest_cand;
            double nearest_dist=NOT_CONNECTED;

            for(auto &fp : free_protons) {
                // find nearest proton or neutron
                for (auto &candidate: candidates) {
                    if ((dist=distance_matrix[make_pair(fp, candidate)]) < nearest_dist) {
                        nearest_dist = dist;
                        free_proton = fp;
                        nearest_cand = candidate;
                    }
                }
            }

            // TODO
//            double angle_a, angle_b;
//            QEdge* nearest_edge;
//            // find nearest free-proton/mst-edge pair
//            for(auto &fp : free_protons) {
//                // find nearest proton or neutron
//                for (auto &edge: mst_edges) {
//                    edge->distance(fp, dist, angle_a, angle_b);
//                    if (dist < nearest_dist){
//                        free_proton = fp;
//                        nearest_dist = dist;
//                        nearest_edge = edge;
//                    }
//                }
//            }
//            cout << "Nearest Edge: " << *nearest_edge << endl;







            // get neighbours of nearest_cand
            vector<QNode *> other;
            for (auto &edge : mst_edges){
                if (edge->contains(nearest_cand))
                    other.push_back(edge->get_othernode(nearest_cand));
            }
            if (other.size() == 1) {      // proton on end of branch
                QEdge *edge = get_edge(nearest_cand, free_proton);
                if (edge != nullptr)
                    mst_edges.push_back(get_edge(nearest_cand, free_proton));
            }
            else if (other.size() == 2) {   // proton in the middle of a branch
                if (distance_matrix[make_pair(free_proton, other[0])] >= NOT_CONNECTED &&
                    distance_matrix[make_pair(free_proton, other[1])] >= NOT_CONNECTED) {
                    cerr << "No edge to neighbours found! Case not implemented" << endl;
                    exit(1);
                } else {
                    short idx;
//                    (distance_matrix[make_pair(free_proton, other[0])] <
//                     distance_matrix[make_pair(free_proton, other[1])]) ? (idx = 0) : (idx = 1);
                    (distance_matrix[make_pair(free_proton, other[0])] - distance_matrix[make_pair(nearest_cand, other[0])] <
                            distance_matrix[make_pair(free_proton, other[1])] - distance_matrix[make_pair(nearest_cand, other[1])]) ? (idx = 0) : (idx = 1);

                    auto it = find(mst_edges.begin(), mst_edges.end(), get_edge(nearest_cand, other[idx]));
                    if (it != mst_edges.end())
                        *it = get_edge(nearest_cand, free_proton);  // replace edge
                    mst_edges.push_back(get_edge(free_proton, other[idx]));    // add second edge
                }
            }
            else if (other.size() == 3) {   // neutron
                if (distance_matrix[make_pair(free_proton, other[0])] >= NOT_CONNECTED &&
                   distance_matrix[make_pair(free_proton, other[1])] >= NOT_CONNECTED &&
                   distance_matrix[make_pair(free_proton, other[2])] >= NOT_CONNECTED) {
                    cerr << "No edge to neighbours found! Case not implemented" << endl;
                    exit(1);
                } else {
                    short idx;
                    (distance_matrix[make_pair(free_proton, other[0])] - distance_matrix[make_pair(nearest_cand, other[0])] <
                     distance_matrix[make_pair(free_proton, other[1])] - distance_matrix[make_pair(nearest_cand, other[1])]) ? (idx = 0) : (idx = 1);
                    if (distance_matrix[make_pair(free_proton, other[2]) ] - distance_matrix[make_pair(nearest_cand, other[2])] <
                                    distance_matrix[make_pair(free_proton, other[idx])] - distance_matrix[make_pair(nearest_cand, other[idx])])
                        idx = 2;
                    auto it = find(mst_edges.begin(), mst_edges.end(), get_edge(nearest_cand, other[idx]));
                    if (it != mst_edges.end())
                        *it = get_edge(nearest_cand, free_proton);  // replace edge
                    mst_edges.push_back(get_edge(free_proton, other[idx]));    // add second edge
                }
            }
            else {
                cerr << "Something wrong!" << endl;
                exit(1);
            }
            // successfully connected
            candidates.insert(free_proton);
            free_protons.erase(free_proton);
        }
    }
}

void QGraph::find_mst_with_neutrons2() {
    if (edges.empty()) {
        cerr << "No edges defined" << endl;
        exit(1);
    }
    mst_edges.clear();
    set<QEdge *> candidate_edges(edges.begin(), edges.end());    // all edges
    unordered_map<QNode *, unsigned short> edge_count;
    // initialize edge_count and candidate_edges
    for (auto &node: nodes) {
        edge_count[node] = connected_edges[node].size();
        if (node->type == QNode::QNodeType::CENTRAL) {
            if (edge_count[node] == 0) {
                cerr << "Central Plant not connected" << endl;
                exit(1);
            } else if (edge_count[node] == 1) {
                candidate_edges.erase(connected_edges[node][0]);    // keep one and only edge to source
                mst_edges.push_back(connected_edges[node][0]);
            }
        } else if (node->type == QNode::QNodeType::NEUTRON) {
            if (edge_count[node] < 3) {
                cerr << "Neutron does not have enough edges" << endl;
                exit(1);
            } else if (edge_count[node] == 3)                             // keep three edges to neutron
                for (auto &edge: connected_edges[node]) {
                    candidate_edges.erase(edge);
                    mst_edges.push_back(edge);
                }
        } else {          // PROTON
            if (edge_count[node] == 0) {
                cerr << "Proton not connected" << endl;
                exit(1);
            } else if (edge_count[node] == 1) {
                candidate_edges.erase(connected_edges[node][0]);    // keep one and only edge to proton
                mst_edges.push_back(connected_edges[node][0]);
            }
        }
    }
    // sort candidates by length (descending)
    vector<QEdge *> sorted_candidates(candidate_edges.begin(), candidate_edges.end());
    sort(sorted_candidates.begin(), sorted_candidates.end(),
         [](const QEdge *a, const QEdge *b) { return a->track->length > b->track->length; });

    // take longest edge and analyze
    for (auto &edge: sorted_candidates) {
        if (candidate_edges.find(edge) != candidate_edges.end()) {      // edge is still in candidate_edges
            // delete edge from candidates
            candidate_edges.erase(edge);

            // check if removal of edge leads to disconnectivity
//            if (edge->nStart->type == QNode::QNodeType::PROTON && edge->nEnd->type == QNode::QNodeType::PROTON
//                && edge_count[edge->nStart] == 2 && edge_count[edge->nEnd] == 2) {
            if (edge->nStart->type == QNode::QNodeType::PROTON && edge->nEnd->type == QNode::QNodeType::PROTON){
                vector<QEdge *> current_edges(mst_edges.begin(), mst_edges.end());
                current_edges.insert(current_edges.end(), candidate_edges.begin(), candidate_edges.end());

                QGraph subgraph;
                subgraph.add(current_edges);
                if (subgraph.get_nb_of_clusters() > 1) {
                    mst_edges.push_back(edge);  // avoid disconnectivity, keep edge
                    continue;                   // go to next edge
                }
            }

            // edge is definitively removed, check consequences
            for (auto &node: {edge->nStart, edge->nEnd}) {
                edge_count[node]--;
                //edge_count[edge->get_othernode(node)]--;
                if ((node->type == QNode::QNodeType::CENTRAL ) && edge_count[node] == 1) {
                    // keep last edge
                    for (auto &conn_edge: connected_edges[node])
                        if (candidate_edges.find(conn_edge) != candidate_edges.end()) {
                            candidate_edges.erase(conn_edge);
                            mst_edges.push_back(conn_edge);
                            break;
                        }
                }
                else if ((node->type == QNode::QNodeType::PROTON) && edge_count[node] == 2) {
                    // keep last edge
                    for (auto &conn_edge: connected_edges[node])
                        if (candidate_edges.find(conn_edge) != candidate_edges.end()) {
                            candidate_edges.erase(conn_edge);
                            mst_edges.push_back(conn_edge);
                            break;
                        }
                }

                else if (node->type == QNode::QNodeType::NEUTRON && edge_count[node] == 3) {
                    // keep remaining edges
                    for (auto &conn_edge: connected_edges[node])
                        if (candidate_edges.find(conn_edge) != candidate_edges.end()) {
                            candidate_edges.erase(conn_edge);
                            mst_edges.push_back(conn_edge);
                        }
                }
            }
        }
    }
    // add remaining edges
    for (auto &edge : candidate_edges)
        mst_edges.push_back(edge);
}





void QGraph::write_mst_geojson(const char *fileout, GeoCoord::CRS crs){
    ofstream out (fileout, ofstream::out);
    if (out)
        out << GeoTrack::write_header();
    for (auto &edge : mst_edges)
        out << edge->track->write_feature(crs) << ",";
    out.seekp(-1, std::ios_base::end);      // remove last character
    out << GeoTrack::write_closing();
}

void QGraph::compute_quantum_numbers(QNode* node, set<QEdge*>& candidate_edges,
                                     const unsigned int current_n, const int current_l, const unsigned int current_k){
    vector<QEdge*> node_edges;
    static unsigned int max_n = 0;

    // find intersection of node->connected_edges and candidate_edges
    for (auto &edge : connected_edges[node]){
        if (candidate_edges.count(edge) > 0)
            node_edges.push_back(edge);
    }

    if (node_edges.empty()){    // end of branch
        if (node != source) {   // write information into QNode
            node->type = QNode::PROTON;
            node->l = current_l;
            node->k = current_k;
        }
    }
    else if (node_edges.size() == 1){   // source or proton
        if (node != source) {
            node->type = QNode::PROTON;
            node->l = current_l;
            node->k = current_k;
        }
        QEdge* edge = node_edges[0];
        edge->l = current_l;
        edge->k = current_k + 1;
        candidate_edges.erase(node_edges[0]);
        // call compute_quantum_numbers() for the other end point of edge
        compute_quantum_numbers(edge->get_othernode(node), candidate_edges, current_n, edge->l, edge->k);
    }
    else if (node_edges.size() == 2){   // neutron
        node->type = QNode::NEUTRON;
        max_n++;
        node->n = max_n;
        node->l = current_l;      // the l of the branch before the neutron (needed for network table)

        for (int i=0; i<2; i++){
            node_edges[i]->l = (int)pow(-1, i+1) * (int)(node->n);  // l=-n if i=0, l=+n if i=1
            node_edges[i]->k = 1;
            candidate_edges.erase(node_edges[i]);
        }
        for (auto &edge : node_edges)
            compute_quantum_numbers(edge->get_othernode(node), candidate_edges, node->n, edge->l, 1);
    }
    else {
        cerr << "Computation of quantum numbers failed" << endl;
        exit(1);
    }
}

void QGraph::compute_quantum_numbers(const bool check_for_errors){
    set<QEdge*> candidate_edges (mst_edges.begin(), mst_edges.end());
    compute_quantum_numbers(source, candidate_edges, 0, 0, 0);

    if (check_for_errors){
        for (auto &edge : mst_edges){
            if (edge->k == 0)
                cerr << "Quantum number k not set for edge " << *edge << endl;
        }
        for (auto &node : nodes){
            if (node->k == 0)
                cerr << "Quantum number k not set for node " << node->id << endl;
            if (node->type == QNode::PROTON && node->sst_id == 0)
                cerr << "SST number not set for proton " << node->id << endl;
        }
    }
}

void QGraph::write_QNet_files(const char* case_dir, const char* case_name){
    //string file = string(case_dir) + "/" + string(case_name);
    //string file ( string(case_dir) + "/" + case_name );     // faster!
    string file;
    file = case_dir;
    file += "/" ;
    file += case_name;
    ofstream networkfile (file + ".csv", ofstream::out);
    if (networkfile) {
        write_QNet_network_file(networkfile);
        networkfile.close();
    }
    else
        cerr << file << ".csv cannot be written" << endl;

    ofstream lengthfile (file + "_length.csv", ofstream::out);
    if (lengthfile) {
        lengthfile << fixed << setprecision(1);
        write_QNet_length_file(lengthfile);
        lengthfile.close();
    }
    else
        cerr << file << "_length.csv cannot be written" << endl;

    ofstream diamfile (file + "_diam.csv", ofstream::out);
    if (diamfile) {
        diamfile << fixed << setprecision(1);
        write_QNet_diam_file(diamfile);
        diamfile.close();
    }
    else
        cerr << file << "_diam.csv cannot be written" << endl;
}

const char orbit[] = {' ', 's', 'p', 'd', 'f'};

void QGraph::write_QNet_network_file(ostream& out){
    unsigned int l_max=0, k_max=0;
    vector<vector<unsigned int>> pos_table;     // columns from 0 to l_max
    vector<vector<unsigned int>> neg_table;     // columns from -1 to -l_max
    vector<unsigned int> pos_last_line;
    vector<unsigned int> neg_last_line;

    // get relevant data from nodes
    for (auto &node : nodes){
        auto * q_node = dynamic_cast<QNode*>(node);
        if (q_node->type == QNode::NEUTRON){
            if (q_node->l >= 0){
                if (pos_last_line.size() < q_node->l + 1)
                    pos_last_line.resize(q_node->l + 1);
                pos_last_line[q_node->l] = q_node->n;
            }
            else {
                if (neg_last_line.size() < -(q_node->l))
                    neg_last_line.resize(-(q_node->l));
                neg_last_line[-(q_node->l) - 1] = q_node->n;
            }
        }
        else if (q_node->type == QNode::PROTON){
            if (q_node->l >= 0){
                if (pos_table.size() < q_node->k) {
                    pos_table.resize(q_node->k);
                    if (q_node->k > k_max)
                        k_max = q_node->k;
                }
                if (pos_table[q_node->k-1].size() < q_node->l + 1) {
                    pos_table[q_node->k - 1].resize(q_node->l + 1);
                    if (q_node->l > l_max)
                        l_max = q_node->l;
                }
                pos_table[q_node->k-1][q_node->l] = q_node->sst_id;
            }
            else {
                if (neg_table.size() < q_node->k) {
                    neg_table.resize(q_node->k);
                    if (q_node->k > k_max)
                        k_max = q_node->k;
                }
                if (neg_table[q_node->k-1].size() < -(q_node->l)) {
                    neg_table[q_node->k - 1].resize(-(q_node->l));
                    if (-(q_node->l) > l_max)
                        l_max = -(q_node->l);
                }
                neg_table[q_node->k-1][-(q_node->l)-1] = node->sst_id;
            }
        }
    }

    // print out
    out << "l\t0";
    for (int l=-(int)l_max; l<=(int)l_max; l++)
        out << "\t" << l;
    out << "\n";

    bool big_enough;   // switch if the vector has the desired size, else print out "0"
    for (int k=1; k<=k_max; k++) {
        if (k_max < 14)     // print s p d f g h i j k l m n o
            out << "k\t" << (char)(orbit[min(k, 4)] + max(0, k-4));
        else                // print integers
            out << "k\t" << k;
        big_enough = (neg_table.size() >= k);
        for (int i = (int)l_max-1; i >= 0; i--){
            if (big_enough && neg_table[k-1].size() > i)
                out << "\t" << neg_table[k - 1][i];
            else
                out << "\t0";
        }
        big_enough = (pos_table.size() >= k);
        for (int i = 0; i <= l_max; i++){
            if (big_enough && pos_table[k-1].size() > i)
                out << "\t" << pos_table[k - 1][i];
            else
                out << "\t0";
        }
        out << "\n";
    }

    out << "n\t0";
    for (int i = (int)l_max-1; i >= 0; i--){
        if (neg_last_line.size() > i)
            out << "\t" << neg_last_line[i];
        else
            out << "\t0";
    }
    for (int i = 0; i <= l_max; i++){
        if (pos_last_line.size() > i)
            out << "\t" << pos_last_line[i];
        else
            out << "\t0";
    }
    out << "\n";
}

void QGraph::write_QNet_length_file(ostream& out){
    unsigned int l_max=0, k_max=0;
    vector<vector<double>> pos_table;     // columns from 0 to l_max
    vector<vector<double>> neg_table;     // columns from -1 to -l_max

    for (auto &edge : mst_edges){
        if (edge->l >= 0){
            if (pos_table.size() < edge->k) {
                pos_table.resize(edge->k);
                if (edge->k > k_max)
                    k_max = edge->k;
            }
            if (pos_table[edge->k - 1].size() <= edge->l) {
                pos_table[edge->k - 1].resize(edge->l + 1);
                if (edge->l > l_max)
                    l_max = edge->l;
            }
            pos_table[edge->k - 1][edge->l] = edge->length;
        }
        else {
            if (neg_table.size() < edge->k) {
                neg_table.resize(edge->k);
                if (edge->k > k_max)
                    k_max = edge->k;
            }
            if (neg_table[edge->k - 1].size() < -(edge->l)) {
                neg_table[edge->k - 1].resize(-(edge->l));
                if (-(edge->l) > l_max)
                    l_max = -(edge->l);
            }
            neg_table[edge->k - 1][-(edge->l) - 1] = edge->length;
        }
    }

    // print out
    out << "-";
    for (int l=-(int)l_max; l<=(int)l_max; l++)
        out << "\t" << l;
    out << "\n";

    bool big_enough;
    for (int k=1; k<=k_max; k++) {
        if (k_max < 14)     // print s p d f g h i j k l m n o
            out << (char) (orbit[min(k, 4)] + max(0, k - 4));
        else                // print integers
            out << k;
        big_enough = (neg_table.size() >= k);
        for (int i = (int)l_max-1; i >= 0; i--){
            if (big_enough && neg_table[k-1].size() > i)
                out << "\t" << neg_table[k - 1][i];
            else
                out << "\t0";
        }
        big_enough = (pos_table.size() >= k);
        for (int i = 0; i <= l_max; i++){
            if (big_enough && pos_table[k-1].size() > i)
                out << "\t" << pos_table[k - 1][i];
            else
                out << "\t0";
        }
        out << "\n";
    }
}

void QGraph::write_QNet_diam_file(ostream& out){
    unsigned int l_max=0, k_max=0;
    vector<vector<unsigned int>> pos_table;     // columns from 0 to l_max
    vector<vector<unsigned int>> neg_table;     // columns from -1 to -l_max

    for (auto &edge : mst_edges){
        if (edge->l >= 0){
            if (pos_table.size() < edge->k) {
                pos_table.resize(edge->k);
                if (edge->k > k_max)
                    k_max = edge->k;
            }
            if (pos_table[edge->k - 1].size() <= edge->l) {
                pos_table[edge->k - 1].resize(edge->l + 1);
                if (edge->l > l_max)
                    l_max = edge->l;
            }
            pos_table[edge->k - 1][edge->l] = edge->DN;
        }
        else {
            if (neg_table.size() < edge->k) {
                neg_table.resize(edge->k);
                if (edge->k > k_max)
                    k_max = edge->k;
            }
            if (neg_table[edge->k - 1].size() < -(edge->l)) {
                neg_table[edge->k - 1].resize(-(edge->l));
                if (-(edge->l) > l_max)
                    l_max = -(edge->l);
            }
            neg_table[edge->k - 1][-(edge->l) - 1] = edge->DN;
        }
    }

    // print out
    out << "-";
    for (int l=-(int)l_max; l<=(int)l_max; l++)
        out << "\t" << l;
    out << "\n";

    bool big_enough;
    for (int k=1; k<=k_max; k++) {
        if (k_max < 14)     // print s p d f g h i j k l m n o
            out << (char) (orbit[min(k, 4)] + max(0, k - 4));
        else                // print integers
            out << k;
        big_enough = (neg_table.size() >= k);
        for (int i = (int)l_max-1; i >= 0; i--){
            if (big_enough && neg_table[k-1].size() > i)
                out << "\t" << neg_table[k - 1][i];
            else
                out << "\t0";
        }
        big_enough = (pos_table.size() >= k);
        for (int i = 0; i <= l_max; i++){
            if (big_enough && pos_table[k-1].size() > i)
                out << "\t" << pos_table[k - 1][i];
            else
                out << "\t0";
        }
        out << "\n";
    }
}
