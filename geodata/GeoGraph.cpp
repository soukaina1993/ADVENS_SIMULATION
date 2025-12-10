//
// Created by cornelia.blanke on 25.06.2024.
//

#include <algorithm>
#include<numeric>
#include <fstream>
#include<set>
#include "GeoGraph.h"
#include "Pixelation.h"



void Node::display(GeoCoord::CRS crs) const {
    cout << "Node ID " << id << "\n";
    if (coord != nullptr) {
        cout << "Coordinate: ";
        cout << coord->as_string(crs);
    }
    cout << endl;
    display_more();
}

template <typename N>
ostream& operator<<(ostream& out, const Edge<N>& edge){
    if (edge.nStart != nullptr && edge.nEnd != nullptr)
        out << "[" << edge.nStart->id << ", " << edge.nEnd->id << "]";
    else
        out << "Nodes of Edge not defined";
    return out;
}

template <typename N>
typename Edge<N>::TrackType Edge<N>::get_tracktype() const{
    if (track == nullptr)
        return Edge::NONE;
    else if (dynamic_cast<GeoFactorTrack*>(track))
        return Edge::FACTOR;
    else if (dynamic_cast<GeoRoutingTrack*>(track))
        return Edge::ROUTING;
    else if (dynamic_cast<GeoReadTrack*>(track))
        return Edge::READ;
    else
        return Edge::BASE;
}

template <typename N>
bool Edge<N>::contains(N* node){
    return (nStart == node || nEnd == node);
}

template <typename N>
N* Edge<N>::get_othernode(const N* node){
    if (nStart == node)
        return nEnd;
    else
        return nStart;
}

template <typename N>
void Edge<N>::distance(const N* node, double& distance, double& angle_a, double& angle_b, const bool from_track, const double tol){
    distance=0.0;
    if (from_track) //TODO
        cerr << "Not yet implemented! Direct connection between nStart and nEnd is used" << endl;

    struct myvector {
        myvector(double ix, double iy) : x(ix), y(iy) {}
        double x, y;
    };
    myvector ab(nEnd->coord->east-nStart->coord->east, nEnd->coord->north-nStart->coord->north);
    myvector ac(node->coord->east-nStart->coord->east, node->coord->north-nStart->coord->north);
    myvector bc(node->coord->east-nEnd->coord->east, node->coord->north-nEnd->coord->north);

    if (track->dist == 0){
        cerr << "Edge has length 0" << endl;
        exit(1);
    }
    else if ((abs(ac.x) > tol || abs(ac.y) > tol) && (abs(bc.x) > tol || abs(bc.y) > tol)) {    // c is neither a nor b
        angle_a = acos((ab.x * ac.x + ab.y * ac.y) / (track->dist * sqrt(pow(ac.x, 2) + pow(ac.y, 2)))) * 180 / M_PI;
        angle_b = acos((- ab.x * bc.x - ab.y * bc.y) / (track->dist * sqrt(pow(bc.x, 2) + pow(bc.y, 2)))) * 180 / M_PI;

        if (angle_a >= 90)            // distance from a
            distance = distance_2d(*nStart->coord, *node->coord);
        else if (angle_b >= 90)       // distance from b
            distance = distance_2d(*nEnd->coord, *node->coord);
        else
            distance = abs(ab.x * ac.y - ab.y * ac.x) / track->dist;      // perpendicular distance from edge
    }
}

template <typename N>
void Edge<N>::display() const {
    cout << *this << endl;
}


/*** --------------------------------- template class Graph -------------------------------- ***/

int Node::nodecounter = 1;

template <typename N, typename E>
Graph<N, E>::Graph(vector<GeoCoord>& coords, const double tol, const typename E::TrackType type){
    vector<GeoCoord*> singlecoords = suppress_duplicate_coords(coords, tol);
    create_nodes_from_coords(singlecoords);
    create_edges_from_nodes(type);
}

template <typename N, typename E>
Graph<N, E>::Graph(const vector<GeoCoord*>& coords, const double tol, const typename E::TrackType type){
    vector<GeoCoord*> singlecoords = suppress_duplicate_coords(coords, tol);
    create_nodes_from_coords(singlecoords);
    create_edges_from_nodes(type);
}

template <typename N, typename E>
Graph<N, E>::Graph(vector<GeoCoord>& coords, const typename E::TrackType type, const double tol) : Graph(coords, tol, type) {}

template <typename N, typename E>
Graph<N, E>::Graph(const vector<GeoCoord*>& coords, const typename E::TrackType type, const double tol) : Graph(coords, tol, type) {}

template <typename N, typename E>
Graph<N, E>::Graph(const vector<N*>& nodevector, const typename E::TrackType type){
    nodes = nodevector;
    create_edges_from_nodes(type);
}

template <typename N, typename E>
Graph<N, E>::Graph(const vector<E*>& edgevector){
    edges = edgevector;
    set<N*> distinct_nodes;
    for (auto &edge : edges){
        for (auto &node : {edge->nStart, edge->nEnd})
            distinct_nodes.insert(node);
    }
    nodes.insert(nodes.end(), distinct_nodes.begin(), distinct_nodes.end());
}

template <typename N, typename E>
Graph<N, E>::Graph(const vector<GeoTrack*>& tracks, const double tol) {
    build_graph(tracks, tol);
}

// Overloaded constructor to create a GeoGraph from a vector of GeoReadTrack pointers.
template <typename N, typename E>
Graph<N, E>::Graph(const vector<GeoReadTrack*>& tracks, const double tol) {
    vector<GeoTrack*> geotracks;
    geotracks.reserve(tracks.size());
    for (auto & track : tracks)
        geotracks.push_back(track);
    build_graph(geotracks, tol);
}



template <typename N, typename E>
Graph<N, E>::~Graph(){
    for (auto &node: nodes_to_delete)
        delete node;
    for (auto &edge: edges_to_delete)
        delete edge;
}

template <typename N, typename E>
void Graph<N, E>::build_graph(const vector<GeoTrack*>& tracks, const double tol){
    Pixelation pixelation;
    unordered_map<GeoCoord*, N*> identify_node;
    vector<GeoCoord*> closest_coord;
    N* current_node;
    E* current_edge;
    short count;

    edges.reserve(tracks.size());
    edges_to_delete.reserve(tracks.size());
    for (auto &track : tracks){
        count = 0;
        // create Edge
        current_edge = new E;
        current_edge->track = track;
        edges.push_back(current_edge);
        edges_to_delete.push_back(current_edge);

        for (auto &coord : {track->start, track->end}){
            closest_coord.clear();
            closest_coord = pixelation.find_closest_to(*coord, 1, tol);
            if (closest_coord.empty()) {
                // create node
                pixelation.add_coord(*coord);
                current_node = new N(coord);
                nodes.push_back(current_node);
                nodes_to_delete.push_back(current_node);
                identify_node[coord] = current_node;
            }
            else {
                // use existing node
                current_node = identify_node.at(closest_coord[0]);  // throws 'std::out_of_range' if not present
                count++;
            }
            if (coord  == track->start)    // i.e. we are in start-loop
                current_edge->nStart = current_node;
            else                                    // i.e. we are in end-loop
                current_edge->nEnd = current_node;

            if (count == 2){    // both nodes were already present, only possible in end-loop
                for(auto &edge : connected_edges[current_node]){
                    if ((edge->nStart == current_edge->nStart) || (edge->nEnd == current_edge->nStart)){
                        // there is already an edge from current_edge->nStart to current_edge->nEnd
                        cerr << "Two tracks from "<< current_edge->nStart->coord->as_string(GeoCoord::CH1903plus);
                        cerr << " to " << current_edge->nEnd->coord->as_string(GeoCoord::CH1903plus) << endl;
                        exit(1);
                    }
                }
            }
            connected_edges[current_node].push_back(current_edge);
        }
    }
}

template <typename N, typename E>
void Graph<N, E>::add(vector<N*>& more_nodes){
    nodes.insert(nodes.end(), more_nodes.begin(), more_nodes.end());
    for (auto &node : more_nodes)
        connected_edges[node] = {};
}

template <typename N, typename E>
void Graph<N, E>::add(vector<E*>& more_edges){
    // add nodes if not yet present
    for (auto &edge : more_edges){
        for (auto &node : {edge->nStart, edge->nEnd}) {
            if ((find(nodes.begin(), nodes.end(), node) == nodes.end()))        // node not found
                nodes.push_back(edge->nStart);
            connected_edges[node].push_back(edge);
        }
    }
    edges.insert( edges.end(), more_edges.begin(), more_edges.end());
}

template <typename N, typename E>
void Graph<N, E>::add(E* edge){
    vector<E*> edgevector = {edge};
    add(edgevector);
}

template <typename N, typename E>
void Graph<N, E>::remove(vector<N*>& less_nodes){
    // complexity O(NlogN)
    unsigned int current_index=0, count=0;
    std::sort(less_nodes.begin(), less_nodes.end());    // pointers are sorted by their addresses
    for (auto &node : nodes) {
        if (!std::binary_search(less_nodes.begin(), less_nodes.end(), node)) {
            // not found, i.e. keep node
            nodes[current_index] = node;
            current_index++;
        } else {
            // remove connected edges
            remove(connected_edges[node]);
            connected_edges.erase(node);
            if ((++count) == less_nodes.size())   // all less_nodes found
                break;
        }
    }
    // erase last nodes (now unused)
    nodes.erase(nodes.begin() + current_index, nodes.end());
}

template <typename N, typename E>
void Graph<N, E>::remove(vector<N>& less_nodes){
    unsigned int current_index=0, count=0;
    std::sort(less_nodes.begin(), less_nodes.end());    // sorted by overloaded operator<
    for (auto &node : nodes) {
        if (!std::binary_search(less_nodes.begin(), less_nodes.end(), *node)) {
            nodes[current_index] = node;
            current_index++;
        } else {
            remove(connected_edges[node]);
            connected_edges.erase(node);
            if ((++count) == less_nodes.size())
                break;
        }
    }
    nodes.erase(nodes.begin() + current_index, nodes.end());
}
template <typename N, typename E>
void Graph<N, E>::remove(vector<int>& less_nodeids){
    unsigned int current_index=0, count=0;
    std::sort(less_nodeids.begin(), less_nodeids.end());    // sort integers
    for (auto &node : nodes){
        if (!std::binary_search (less_nodeids.begin(), less_nodeids.end(), node->id)) {
            nodes[current_index] = node;
            current_index++;
        }
        else {
            remove(connected_edges[node]);
            connected_edges.erase(node);
            if ((++count) == less_nodeids.size())
                break;
        }
        nodes.erase(nodes.begin() + current_index, nodes.end());
    }
}
template <typename N, typename E>
void Graph<N, E>::remove(vector<unsigned int>& less_nodeids){
    unsigned int current_index=0, count=0;
    std::sort(less_nodeids.begin(), less_nodeids.end());
    for (auto &node : nodes){
        if (!std::binary_search (less_nodeids.begin(), less_nodeids.end(), node->id)) {
            nodes[current_index] = node;
            current_index++;
        }
        else {
            remove(connected_edges[node]);
            connected_edges.erase(node);
            if ((++count) == less_nodeids.size())
                break;
        }
        nodes.erase(nodes.begin() + current_index, nodes.end());
    }
}
template <typename N, typename E>
void Graph<N, E>::remove(N* node){
    vector<N*> vect = {node};
    remove(vect);
}
template <typename N, typename E>
void Graph<N, E>::remove(N& node){
    vector<N> vect = {node};
    remove(vect);
}
template <typename N, typename E>
void Graph<N, E>::remove(vector<E*>& less_edges){
    unsigned int current_index=0, count=0;
    typename vector<E*>::iterator it;
    std::sort(less_edges.begin(), less_edges.end());
    for (auto &edge : edges) {
        if (!std::binary_search(less_edges.begin(), less_edges.end(), edge)) {
            // not found, i.e. keep edge
            edges[current_index] = edge;
            current_index++;
        }
        else {
            for (auto &node : {edge->nStart, edge->nEnd}) {
                it = find(connected_edges[node].begin(), connected_edges[node].end(), edge);
                connected_edges[node].erase(it);
            }
            if ((++count) == less_edges.size())
                break;      // ell less_edges found
        }
        // erase last edges
        edges.erase(edges.begin() + current_index, edges.end());
    }
}
template <typename N, typename E>
void Graph<N, E>::remove(vector<E>& less_edges){
    unsigned int current_index=0, count=0;
    typename vector<E*>::iterator it;
    std::sort(less_edges.begin(), less_edges.end());
    for (auto &edge : edges) {
        if (!std::binary_search(less_edges.begin(), less_edges.end(), *edge)) {
            // not found, i.e. keep edge
            edges[current_index] = edge;
            current_index++;
        }
        else {
            for (auto &node : {edge->nStart, edge->nEnd}) {
                it = find(connected_edges[node].begin(), connected_edges[node].end(), edge);
                connected_edges[node].erase(it);
            }
            if ((++count) == less_edges.size())
                break;
        }
        // erase last edges
        edges.erase(edges.begin() + current_index, edges.end());
    }
}
template <typename N, typename E>
void Graph<N, E>::remove(N* A, N* B){
    typename vector<E*>::iterator it, it2;
    for (it = edges.begin(); it != edges.end(); it++) {
        if (((*it)->nStart == A || (*it)->nEnd == A) && ((*it)->nStart == B || (*it)->nEnd == B)) { // found
            for (auto &node : {(*it)->nStart, (*it)->nEnd}) {
                it2 = find(connected_edges[node].begin(), connected_edges[node].end(), (*it));
                connected_edges[node].erase(it2);
            }
            edges.erase(it);
            break;      // no duplicate edges
        }
    }
}
template <typename N, typename E>
void Graph<N, E>::remove(N& A, N& B){
    typename vector<E*>::iterator it, it2;
    for (it = edges.begin(); it != edges.end(); it++) {
        if ((*(*it)->nStart == A || *(*it)->nEnd == A) && (*(*it)->nStart == B || *(*it)->nEnd == B)) { // found
            for (auto &node : {(*it)->nStart, (*it)->nEnd}) {
                it2 = find(connected_edges[node].begin(), connected_edges[node].end(), (*it));
                connected_edges[node].erase(it2);
            }
            edges.erase(it);
            break;      // no duplicate edges
        }
    }
}
template <typename N, typename E>
void Graph<N, E>::remove(unsigned int A, unsigned int B){
    typename vector<E*>::iterator it, it2;
    for (it = edges.begin(); it != edges.end(); it++) {
        if (((*it)->nStart->id == A || (*it)->nEnd->id == A) && ((*it)->nStart->id == B || (*it)->nEnd->id == B)) { // found
            for (auto &node : {(*it)->nStart, (*it)->nEnd}) {
                it2 = find(connected_edges[node].begin(), connected_edges[node].end(), (*it));
                connected_edges[node].erase(it2);
            }
            edges.erase(it);
            break;      // no duplicate edges
        }
    }
}

template <typename N, typename E>
E* Graph<N, E>::get_edge(N* A, N* B){
    typename vector<E*>::iterator it;

    if (connected_edges.count(A) > 0) {
        for (it = connected_edges[A].begin(); it != connected_edges[A].end(); it++) {
            if (((*it)->nStart == A || (*it)->nEnd == A) && ((*it)->nStart == B || (*it)->nEnd == B))  // found
                return *it;
        }
    }
    // if not yet found, search in all edges
    for (it = edges.begin(); it != edges.end(); it++) {
        if (((*it)->nStart == A || (*it)->nEnd == A) && ((*it)->nStart == B || (*it)->nEnd == B))  // found
            return *it;
    }
    return nullptr;     //not found
}

template <typename N, typename E>
E* Graph<N, E>::get_edge(N& A, N& B){
    typename vector<E*>::iterator it;

    if (connected_edges.count(&A) > 0) {
        for (it = connected_edges[&A].begin(); it != connected_edges[&A].end(); it++) {
            if ((*((*it)->nStart) == A || *((*it)->nEnd) == A) && (*((*it)->nStart) == B || *((*it)->nEnd) == B))  // found
                return *it;
        }
    }
    // if not yet found, search in all edges
    for (it = edges.begin(); it != edges.end(); it++) {
        if ((*((*it)->nStart) == A || *((*it)->nEnd) == A) && (*((*it)->nStart) == B || *((*it)->nEnd) == B))  // found
            return *it;
    }
    return nullptr;     //not found
}

template <typename N, typename E>
vector<N*> Graph<N, E>::get_allnodes(){
    return nodes;
}

template <typename N, typename E>
void Graph<N, E>::display_nodeinfo(N* node) {
    node->display();
    cout << "Connected Edges: ";
    for (auto &edge : connected_edges[node])
        cout << *edge << " ";
    cout << "\n";
}

template <typename N, typename E>
void Graph<N, E>::display_nodeinfo(N& node) {
    display_nodeinfo(&node);
}

template <typename N, typename E>
void Graph<N, E>::display_nodeinfo(const int id) {
    for (auto& node : nodes){
        if (node->id == id)
            display_nodeinfo(node);
    }
}

template <typename N, typename E>
void Graph<N, E>::display() const {
    cout << nodes.size() << " Nodes: ";
    for (auto &node : nodes)
        cout << node->id << " ";
    cout << "\n" << edges.size() << " Edges: ";
    for (auto &edge : edges)
        cout << *edge << " ";
    cout << "\n";
}

template <typename N, typename E>
void Graph<N, E>::check_connected_edges(){
    unordered_map<N*, unsigned int> edge_count;
    if (connected_edges.size() > nodes.size())
        cerr << "Too many entries in connected_edges" << endl;
    else {
        for (auto& node : nodes)
            if (!connected_edges.count(node))
                cerr << "Node ID " << node->id << " not found in connected_edges" << endl;
        for (auto& edge : edges) {
            for (auto& node : {edge->nStart, edge->nEnd}) {
                if (find(connected_edges[node].begin(), connected_edges[node].end(), edge) == connected_edges[node].end())
                    cerr << "Edge " << *edge << " not found in connected_edges of Node ID " << node->id << endl;
                else
                    edge_count[node]++;
            }
        }
        for (const auto& it : connected_edges)
            if (it.second.size() > edge_count[it.first])
                cerr << "Node ID " << it.first->id << " has too many edges in connected_edges" << endl;
    }
}

template <typename N, typename E>
void Graph<N, E>::compute_connected_edges(){
    connected_edges.clear();
    for (auto& node : nodes)
        connected_edges[node] = {};
    for (auto& edge : edges)
        for (auto& node : {edge->nStart, edge->nEnd})
            connected_edges[node].push_back(edge);
}

template <typename N, typename E>
bool Graph<N, E>::check_connectivity(){
    set<N*> visited_nodes, current_nodes, next_nodes;
    N* other;
    vector<N*> allnodes = get_allnodes();
    current_nodes.insert(allnodes[0]);
    visited_nodes = current_nodes;

    while (!current_nodes.empty()){
        for (auto &node: current_nodes) {
            for (auto &edge: connected_edges[node]) {
                other = edge->get_othernode(node);
                if (visited_nodes.find(other) == visited_nodes.end()) {     // not yet visited
                    visited_nodes.insert(other);
                    next_nodes.insert(other);
                }
            }
        }
        current_nodes = next_nodes;
        next_nodes.clear();
    }
    return (visited_nodes.size() == allnodes.size());   // false or true
}

template <typename N, typename E>
int Graph<N, E>::get_nb_of_clusters(){
    int smaller, greater, maxidx=0, count=0;
    map<N*, unsigned int> idx;
    for (auto &node : nodes)    // initialize with 0
        idx[node] = 0;

    for (auto &edge : edges){
        if (idx[edge->nStart] == idx[edge->nEnd]) {
            if (idx[edge->nStart] == 0) {                             // new cluster
                idx[edge->nStart] = idx[edge->nEnd] = ++maxidx;
                count++;
            }
            // else edge is inside same cluster, nothing to do
        }
        else if (idx[edge->nStart] == 0)
            idx[edge->nStart] = idx[edge->nEnd];
        else if (idx[edge->nEnd] == 0)
            idx[edge->nEnd] = idx[edge->nStart];
        else {                          // connect two different clusters
            if (idx[edge->nStart] < idx[edge->nEnd]){
                smaller = idx[edge->nStart];
                greater = idx[edge->nEnd];
            }
            else {
                greater = idx[edge->nStart];
                smaller = idx[edge->nEnd];
            }
            idx[edge->nStart] = idx[edge->nEnd] = smaller;
            for (auto it = idx.begin(); it != idx.end(); ++it)
                if (it->second == greater)
                    it->second = smaller;
            count--;
        }
    }
    return count;
}


template <typename N, typename E>
vector<GeoCoord*> Graph<N, E>::suppress_duplicate_coords(vector<GeoCoord>& coords, const double tol){
    if (!coords.empty()) {
        Pixelation pixelation(coords[0]);
        for (auto &coord: coords) {
            if (pixelation.find_closest_to(coord, 1, tol).empty())
                pixelation.add_coord(coord);
        }
        return pixelation.get_coords();
    }
    else {
        vector<GeoCoord*> vect;
        return vect;
    }
}

template <typename N, typename E>
vector<GeoCoord*> Graph<N, E>::suppress_duplicate_coords(const vector<GeoCoord*>& coords, const double tol){
    if (!coords.empty()) {
        Pixelation pixelation(*coords[0]);
        for (auto &coord : coords) {
            if(pixelation.find_closest_to(*coord, 1, tol).empty())
                pixelation.add_coord(*coord);
        }
        return pixelation.get_coords();
    }
    else {
        vector<GeoCoord*> vect;
        return vect;
    }
}

template <typename N, typename E>
void Graph<N, E>::create_nodes_from_coords(const vector<GeoCoord*>& coords){
    nodes.reserve(nodes.size() + coords.size());
    nodes_to_delete.reserve(nodes_to_delete.size() + coords.size());
    for (auto &coord : coords) {
        N* node = new N(coord);
        nodes.push_back(node);
        nodes_to_delete.push_back(node);
        connected_edges[node] = {};
    }
}

template <typename N, typename E>
void Graph<N, E>::create_edges_from_nodes(const typename E::TrackType type){
    if (nodes.empty()){
        cerr << "Init nodes first" << endl;
        exit(1);
    }
    if (type == E::BASE)
        connect_nodes(nodes, type);
    else if (type == E::FACTOR)
        connect_nodes(nodes, type);     //TODO track->compute_penalty()
    else if (type == E::ROUTING)
        connect_nodes_with_routing(nodes);    // more efficient than tracking one by one
}

template <typename N, typename E>
void Graph<N, E>::merge_nodes(const double tol) {
    Pixelation pixelation(*nodes[0]->coord);
    unordered_map<GeoCoord*, N*> identify_node;
    vector<GeoCoord*> vect;
    N* node_to_keep;

    for (auto &node : nodes) {
        if(!(vect = pixelation.find_closest_to(*node->coord, 1, tol)).empty()) {
            // coord already there --> merge nodes
            node_to_keep = identify_node[vect[0]];
            for (auto &edge : connected_edges[node]) {
                // exchange node with node_to_keep
                if(edge->nStart == node)
                    edge->nStart = node_to_keep;
                else
                    edge->nEnd = node_to_keep;
                connected_edges[node_to_keep].push_back(edge);
            }
            connected_edges[node].clear();
            remove(node);
        }
        else {
            pixelation.add_coord(*node->coord);
            identify_node[node->coord] = node;
        }
    }
}

template <typename N, typename E>
void Graph<N, E>::connect_nodes(vector<N*>& nodevector, const typename E::TrackType type){
    int size=(int)nodevector.size();
    unsigned long long newsize = edges.size() + (size * (size - 1)) / 2;
    edges.reserve(newsize);
    edges_to_delete.reserve(newsize);
    for (int i=0; i<size-1; i++){
        for (int j=i+1; j<size; j++){
            E* edge = new E(nodevector[i], nodevector[j], type);
            edges.push_back(edge);
            edges_to_delete.push_back(edge);
            connected_edges[nodevector[i]].push_back(edge);
            connected_edges[nodevector[j]].push_back(edge);
        }
    }
}
template <typename N, typename E>
void Graph<N, E>::connect_nodes(vector<N*>& nodevector, N* node, const typename E::TrackType type){
    int size=(int)nodevector.size();
    unsigned long long newsize = edges.size() + size;
    edges.reserve(newsize);
    edges_to_delete.reserve(newsize);
    for (int i=0; i<size; i++){
        E* edge = new E(nodevector[i], node, type);
        edges.push_back(edge);
        edges_to_delete.push_back(edge);
        connected_edges[nodevector[i]].push_back(edge);
        connected_edges[node].push_back(edge);
    }
}

template <typename N, typename E>
void Graph<N, E>::connect_nodes(vector<N*>& nodevector1, vector<N*>& nodevector2, const typename E::TrackType type){
    for (auto &node : nodevector2)
        connect_nodes(nodevector1, node, type);
}


template <typename N, typename E>
void Graph<N, E>::connect_nodes_with_routing(const vector<N*>& nodevector){
    // prepare curl call
    ostringstream curl;
    streamsize ss = curl.precision();
    curl << fixed << setprecision(10);
    //curl << "curl \"https://routing.openstreetmap.de/routed-foot/table/v1/driving/";    // by foot
    curl << "curl \"http://localhost:5000/table/v1/driving/";
    for (auto &node : nodevector)
        curl << node->coord->lon << "," << node->coord->lat << ";";
    curl.seekp(-1, std::ios_base::end);      // remove last character
    curl << "?annotations=distance";
    curl << "&generate_hints=false";
    curl << "&skip_waypoints=true";
    curl << "\"2>&0";                   // 2>&0 don't show error message
    curl << defaultfloat << setprecision(int(ss));
    const string strcurl = std::move(curl).str();
    //cout << strcurl << endl;

    string result = GeoRoutingTrack::execute_curl(strcurl);
    //cout << result << endl;

    // process result (depends on output format)
    if (result.size() > 11 && result.substr(1, 11)==R"("code":"Ok")"){
        stringstream strstr(result.substr(25, result.size() - 26));
        //cout << strstr.str() << endl;
        string row, length;
        unsigned int i=0, j;
        while (i < nodevector.size() && getline(strstr, row, ']')) {
            j = 0;
            stringstream str(row.substr(2));
            while (getline(str, length, ',')) {
                if (j > i) {
                    E* edge = new E(nodevector[i], nodevector[j], E::ROUTING);
                    edge->track->length = stod(length);
                    edges.push_back(edge);
                    edges_to_delete.push_back(edge);
                    connected_edges[nodevector[i]].push_back(edge);
                    connected_edges[nodevector[j]].push_back(edge);
                }
                j++;
            }
            i++;
        }
    }
    else {
        cerr << "Routing not successful." << endl;
    }
}

template <typename N, typename E>
void Graph<N, E>::connect_nodes_with_routing(vector<N*>& nodevector, N* node){
    vector<N*> nodevector2 = {node};
    connect_nodes_with_routing(nodevector, nodevector2);
    //connect_nodes_with_routing(nodevector2, nodevector);
}

template <typename N, typename E>
void Graph<N, E>::connect_nodes_with_routing(vector<N*>& nodevector1, vector<N*>& nodevector2){
    // prepare curl call
    ostringstream curl;
    streamsize ss = curl.precision();
    curl << fixed << setprecision(10);
    // curl << "curl \"https://routing.openstreetmap.de/routed-foot/table/v1/driving/";    // by foot
    curl << "curl \"http://localhost:5000/table/v1/driving/";
    for (auto &pnode : nodevector1)
        curl << pnode->coord->lon << "," << pnode->coord->lat << ";";
    for (auto &pnode : nodevector2)
        curl << pnode->coord->lon << "," << pnode->coord->lat << ";";
    curl.seekp(-1, std::ios_base::end);      // remove last character

    unsigned int idx=0, size=nodevector1.size();
    curl << "?sources=";
    while (idx < size)
        curl << idx++ << ";";
    curl.seekp(-1, std::ios_base::end);      // remove last character

    size += nodevector2.size();
    curl << "&destinations=";
    while (idx < size)
        curl << idx++ <<";";
    curl.seekp(-1, std::ios_base::end);      // remove last character

    curl << "&annotations=distance";
    curl << "&generate_hints=false";
    curl << "&skip_waypoints=true";
    curl << "\"2>&0";                   // 2>&0 don't show error message
    curl << defaultfloat << setprecision(int(ss));
    const string strcurl = std::move(curl).str();

    string result = GeoRoutingTrack::execute_curl(strcurl);

    // process result (depends on output format)
    if (result.size() > 11 && result.substr(1, 11)==R"("code":"Ok")"){
        stringstream strstr(result.substr(25, result.size() - 26));
        cout << strstr.str() << endl;
        string row, length;
        unsigned int i=0, j;
        size = nodevector1.size();
        while (i < size && getline(strstr, row, ']')) {
            j = 0;
            stringstream str(row.substr(2));
            while (getline(str, length, ',')) {
                E* edge = new E(nodevector1[i], nodevector2[j], E::ROUTING);
                edge->track->length = stod(length);
                edges.push_back(edge);
                edges_to_delete.push_back(edge);
                connected_edges[nodevector1[i]].push_back(edge);
                connected_edges[nodevector2[i]].push_back(edge);
                j++;
            }
            i++;
        }
    }
    else {
        cerr << "Routing not successful." << endl;
    }
}

// Note: implemented using AI and https://www.geeksforgeeks.org/ to understand the algorithm
template <typename N, typename E>
vector<vector<N*>> Graph<N, E>::find_connected_subgraphs() {
    vector<vector<N*>> connectedSubgraphs;
    set<N*> visited;

    // Helper function for DFS
    auto dfs = [&](N* node, vector<N*>& component, auto&& dfs_ref) -> void {
        visited.insert(node);
        component.push_back(node);
        for (auto &edge : connected_edges[node]) {
            N* neighbor = (edge->nStart == node) ? edge->nEnd : edge->nStart;
            if (visited.find(neighbor) == visited.end())
                dfs_ref(neighbor, component, dfs_ref);
        }
    };
    // Iterate through all nodes to find connected components
    for (auto &node : nodes) {
        if (visited.find(node) == visited.end()) {
            std::vector<N*> component;
            dfs(node, component, dfs); // Start a DFS for this component
            connectedSubgraphs.push_back(component);
        }
    }
    return connectedSubgraphs;
}

template <typename N, typename E>
vector<Graph<N, E>> Graph<N, E>::extract_subgraphs() {
    vector<Graph<N, E>> subgraphs;
    vector<vector<N*>> connectedSubgraphs = find_connected_subgraphs();

    for (const auto& component : connectedSubgraphs) {
        Graph<N, E> subgraph;

        // Filter nodes for this component
        for (N* node : component) {
            // Add the node to the new subgraph's nodes vector
            if (std::find(subgraph.nodes.begin(), subgraph.nodes.end(), node) == subgraph.nodes.end())
                subgraph.nodes.push_back(node);
        }

        // Filter edges for this component
        for (auto &edge : edges) {
            bool startInComponent = std::find(component.begin(), component.end(), edge->nStart) != component.end();
            bool endInComponent = std::find(component.begin(), component.end(), edge->nEnd) != component.end();

            // Add edge if both nodes belong to this component
            if (startInComponent && endInComponent) {
                if (std::find(subgraph.edges.begin(), subgraph.edges.end(), edge) == subgraph.edges.end())
                    subgraph.edges.push_back(edge);
            }
        }
        subgraphs.push_back(subgraph);
    }
    return subgraphs;
}














template <typename N, typename E>
vector<E*> Graph<N ,E>::find_shortest_edges(N* node, unsigned int n_lowest) {
    vector<double> lengths;
    vector<E*> result;
    n_lowest = min(n_lowest, (unsigned int)(connected_edges[node].size()));
    lengths.reserve(connected_edges[node].size());
    result.reserve(n_lowest);

    for (auto &edge : connected_edges[node])
        lengths.push_back(edge->track->length);
    vector<unsigned int> indices (lengths.size());
    iota(indices.begin(), indices.end(), 0);
    partial_sort(indices.begin(), indices.begin() + n_lowest, indices.end(),
                 [&lengths](const unsigned int &a, const unsigned int &b){return lengths[a] < lengths[b];});
    for (int i=0; i<n_lowest; i++)
        result.push_back(connected_edges[node][indices[i]]);
    return result;
}

template <typename N, typename E>
map<pair<N*, N*>, double> Graph<N ,E>::get_adjacencyMatrix(const vector<N*>& nodevector, unordered_map<N*, vector<E*>>& map_connected_edges){
    map<pair<N*, N*>, double> adjacency;
    N* other;
    for (auto &node : nodevector) {
        for (auto &node2 : nodevector)
            if (node != node2)      // diagonal not needed
                adjacency[make_pair(node, node2)] = NOT_CONNECTED;
        if (map_connected_edges.count(node)) {
            for (auto &edge: map_connected_edges[node]) {
                other = edge->get_othernode(node);
                //(edge->nStart == node) ? (other = edge->nEnd) : (other = edge->nStart);
                if (edge->track == nullptr) {
                    cerr << "Edge " << *edge << " has no track" << endl;
                    exit(1);
                } else
                    adjacency[make_pair(node, other)] = edge->track->length;
            }
        }
    }
    return adjacency;
}

template <typename N, typename E>
map<pair<N*, N*>, double> Graph<N ,E>::get_adjacencyMatrix(){
    vector<N*> allnodes = get_allnodes();
    return get_adjacencyMatrix(allnodes, connected_edges);
}

template <typename N, typename E>
map<pair<N*, N*>, double> Graph<N ,E>::get_adjacencyMatrix_directed(const vector<N*>& nodevector, const vector<E*>& edgevector){
    map<pair<N*, N*>, double> adjacency;
    // initialize adjacency
    for (auto &node : nodevector) {
        for (auto &node2: nodevector)
            if (node != node2)      // diagonal not needed
                adjacency[make_pair(node, node2)] = NOT_CONNECTED;
    }
    // fill values
    for (auto &edge : edgevector) {
        if (edge->track == nullptr){
            cerr << "Edge " << edge << " has no track" << endl;
            exit(1);
        }
        else if (find(nodevector.begin(), nodevector.end(), edge->nStart) == nodevector.end()){
            cerr << "Node " << edge->nStart->id << " not found" << endl;
            exit(1);
        }
        else if (find(nodevector.begin(), nodevector.end(), edge->nEnd) == nodevector.end()){
            cerr << "Node " << edge->nEnd->id << " not found" << endl;
            exit(1);
        }
        else
            adjacency[make_pair(edge->nStart, edge->nEnd)] = edge->track->length;
    }
    return adjacency;
}

template <typename N, typename E>
map<pair<N*, N*>, double> Graph<N ,E>::FloydWarshall_dir() {
    map<pair<N*, N*>, double> shortest;
    map<pair<N*, N*>, N*> next_node;
    FloydWarshall_directed(shortest, next_node, false);
    return shortest;
}

template <typename N, typename E>
void Graph<N ,E>::FloydWarshall_directed(map<pair<N*, N*>, double>& distance_matrix, map<pair<N*, N*>, N*>& next_node,
                                                  const bool compute_next) {
    distance_matrix = get_adjacencyMatrix_directed(nodes, edges);
    typename map<pair<N*, N*>, double>::iterator it;
    set<N*> snode;
    for(it = distance_matrix.begin(); it != distance_matrix.end(); ++it)
        snode.insert(it->first.first);
    vector<N*> vnode(snode.begin(), snode.end());   // collect all nodes
    snode.clear();

    int i, j;
    pair<N*, N*> ij, ik, kj;
    double dist;

    // initialize next_node
    if (compute_next) {
        next_node.clear();
        for(it = distance_matrix.begin(); it != distance_matrix.end(); ++it)
            next_node[it->first] = nullptr;
        for (auto &edge: edges) {
            try {
                next_node.at(make_pair(edge->nStart, edge->nEnd)) = edge->nEnd;
            }
            catch (const std::out_of_range& ex) {
                cerr << "Node of edge " << edge << " not found" << endl;
                exit(1);
            }
        }
    }
    // Floyd-Warshall algorithm
    for (auto &node_k : vnode)
        // for all pairs (i,j): check if path via node_k is shorter
        for (i = 0; i < vnode.size(); i++)
            for (j = 0; j < vnode.size(); j++) {
                if (i != j) {
                    ij = make_pair(vnode[i], vnode[j]);
                    ik = make_pair(vnode[i], node_k);
                    kj = make_pair(node_k, vnode[j]);
                    if ((vnode[i] != node_k) && (vnode[j] != node_k) &&
                        (dist = distance_matrix[ik] + distance_matrix[kj]) < distance_matrix[ij]) {
                        distance_matrix[ij] = dist;
                        if (compute_next)
                            next_node[ij] = next_node[ik];
                    }
                }
            }
}

template <typename N, typename E>
map<pair<N*, N*>, double> Graph<N ,E>::FloydWarshall_undir() {
    map<pair<N*, N*>, double> shortest;
    map<pair<N*, N*>, N*> next_node;
    FloydWarshall_undirected(shortest, next_node, false);
    return shortest;
}

template <typename N, typename E>
void Graph<N ,E>::FloydWarshall_undirected(map<pair<N*, N*>, double>& distance_matrix, map<pair<N*, N*>, N*>& next_node,
                                const bool compute_next) {
    typename map<pair<N*, N*>, double>::iterator it;
    set<N*> snode;
    distance_matrix = get_adjacencyMatrix();
    for(it = distance_matrix.begin(); it != distance_matrix.end(); ++it)
        snode.insert(it->first.first);
    vector<N*> vnode(snode.begin(), snode.end());   // collect all nodes
    snode.clear();

    int i, j;
    pair<N*, N*> ij, ji, ik, jk, kj;
    double dist;

    // initialize next_node
    if (compute_next) {
        next_node.clear();
        for(it = distance_matrix.begin(); it != distance_matrix.end(); ++it)
            next_node[it->first] = nullptr;

        for (auto &edge: edges) {
            try {
                next_node.at(make_pair(edge->nStart, edge->nEnd)) = edge->nEnd;
                next_node.at(make_pair(edge->nEnd, edge->nStart)) = edge->nStart;
            }
            catch (const std::out_of_range& ex) {
                cerr << "Node of edge " << edge << " not found" << endl;
                exit(1);
            }
        }
    }
    for (auto &node_k : vnode)
        // for all pairs (i,j): check if path via node_k is shorter
        for (i = 0; i < vnode.size(); i++)
            for (j = i+1; j < vnode.size(); j++) {
                ij = make_pair(vnode[i], vnode[j]);
                ji = make_pair(vnode[j], vnode[i]);
                ik = make_pair(vnode[i], node_k);
                kj = make_pair(node_k, vnode[j]);
                jk = make_pair(vnode[j], node_k);
                if ((vnode[i]!=node_k) && (vnode[j]!=node_k) && (dist = distance_matrix[ik] + distance_matrix[kj]) < distance_matrix[ij]) {
                    distance_matrix[ij] = dist;
                    distance_matrix[ji] = dist;
                    if (compute_next) {
                        next_node[ij] = next_node[ik];
                        next_node[ji] = next_node[jk];
                    }
                }
            }
}

template <typename N, typename E>
vector<N*> Graph<N, E>::reconstruct_path(map<pair<N*, N*>, N*>& next_node, N* A, N* B){
    vector<N*> path = {A};
    int count=0;
    try {
        N* current = A;
        while (path.back() != B && count < next_node.size()) {
            current = next_node.at(make_pair(current, B));
            path.push_back(current);
            count++;
        }
    }
    catch (const std::out_of_range& ex) {
        cerr << "Reconstruction of path failed" << endl;
        exit(1);
    }
    return path;
}

template <typename N, typename E>
vector<E*> Graph<N, E>::get_edges_on_path(const vector<N*> path, unordered_map<N*, vector<E*>>& map_connected_edges){
    vector<E*> edgevector;
    typename vector<E*>::iterator it;
    for (int i=1; i<path.size(); i++){
        map_connected_edges[path[i-1]];
        if (map_connected_edges.count(path[i-1]) > 0) {
            for (it = map_connected_edges[path[i-1]].begin(); it != map_connected_edges[path[i-1]].end(); it++) {
                if (((*it)->nStart == path[i-1] || (*it)->nEnd == path[i-1]) &&
                                    ((*it)->nStart == path[i] || (*it)->nEnd == path[i])) {  // found
                    edgevector.push_back(*it);
                    break;
                }
                else if (it == map_connected_edges[path[i-1]].end()){
                    cerr << "Edge not found" << endl;
                    exit(1);
                }
            }
        }
    }
    return edgevector;
}

template <typename N, typename E>
vector<E*> Graph<N, E>::get_edges_on_path(const vector<N*> path){
    return get_edges_on_path(path, connected_edges);
}


template <typename N, typename E>
void Graph<N, E>::write_geojson(const char *fileout, GeoCoord::CRS crs){
    ofstream out (fileout, ofstream::out);
    if (out) {
        out << GeoTrack::write_header();
        for (auto &edge: edges)
            out << edge->track->write_feature(crs) << ",";
        out.seekp(-1, std::ios_base::end);      // remove last character
        out << GeoTrack::write_closing();
    }
}


#include "QGraph.h"
class QNode;
class QEdge;
template class Edge<Node>;
template class Edge<QNode>;
template class Graph<Node, Edge<Node>>;
template class Graph<QNode, QEdge>;