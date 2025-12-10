//
// Created by cornelia.blanke on 01.07.2024.
//

#include<algorithm>
#include<numeric>
#include "Pixelation.h"

Pixelation::Pixelation(GeoCoord& myorigin, const double precision) :
    origin(myorigin), bottom_left(myorigin), top_right(myorigin), precision(precision){
}

Pixelation::Pixelation(vector<GeoCoord>& coords, const double precision) :
    origin(coords[0]), bottom_left(coords[0]), top_right(coords[0]), precision(precision) {
        compute(coords, precision);
}

Pixelation::Pixelation(vector<GeoCoord>& coords, GeoCoord& myorigin, const double precision) :
    origin(myorigin), bottom_left(myorigin), top_right(myorigin), precision(precision) {
        compute(coords, precision);
}


void Pixelation::compute(vector<GeoCoord>& coords, const double prec){
    Index index;
    pixel.clear();
    precision = prec;
    for (auto& coord : coords){
        index = get_index(coord);
        add_to_pixelation(index, &coord);
    }
}

void Pixelation::compute(vector<GeoCoord*>& ptr_coords, const double prec){
    Index index;
    pixel.clear();
    precision = prec;
    for (auto& pcoord : ptr_coords){
        index = get_index(*pcoord);
        add_to_pixelation(index, pcoord);
    }
}

vector<GeoCoord*> Pixelation::get_coords(){
    vector<GeoCoord*> ptr_coords;
    for(auto& pix : pixel){
        ptr_coords.insert(ptr_coords.end(), pix.second.begin(), pix.second.end());
    }
    return ptr_coords;
}

void Pixelation::modify_precision(const int prec){
    vector<GeoCoord*> ptr_coords = get_coords();
    compute(ptr_coords, prec);
}

void Pixelation::modify_origin(GeoCoord& new_origin){
    origin=new_origin;
    vector<GeoCoord*> ptr_coords = get_coords();
    compute(ptr_coords, precision);
}

void Pixelation::add_coord(GeoCoord& coord){
    Index index = get_index(coord);
    add_to_pixelation(index, &coord);
}

void Pixelation::add_to_pixelation(const Index& key, GeoCoord* value){
    if (pixel.count(key)==0)
        pixel[key] = {value};
    else
        pixel[key].push_back(value);

    bool bl_changed=false, tr_changed=false;
    if (value->east < bottom_left.east){
        bottom_left.east = value->east;
        bl_changed=true;
    }
    else if (value->east > top_right.east) {
        top_right.east = value->east;
        tr_changed=true;
    }
    if (value->north < bottom_left.north) {
        bottom_left.north = value->north;
        bl_changed=true;
    }
    else if (value->north > top_right.north) {
        top_right.north = value->north;
        tr_changed=true;
    }
    if (bl_changed)
        bottom_left.compute_wgs84();
    if (tr_changed)
        top_right.compute_wgs84();
}

Index Pixelation::get_index(const GeoCoord& coord) const {
    Index index(int((coord.east - origin.east) / precision), int((coord.north - origin.north) / precision));
    return index;
}

vector<GeoCoord*> Pixelation::area_around(const Index& index, int radius) const{
    vector<GeoCoord*> coords_in_area;
//    const unsigned int max_count = (int)pow(2 * radius + 1,2);
//    unsigned int count=0;
//    for(const auto& pix : pixel) {
//        if ((abs(pix.first.i - index.i) <= radius) && (abs(pix.first.j - index.j) <= radius)) {
//            coords_in_area.insert(coords_in_area.end(), pix.second.begin(), pix.second.end());
//            count++;
//            if (count == max_count)     // found all pixels of interest
//                break;
//        }
//    }
    for (int i=index.i-radius; i<=index.i+radius; i++){
        for (int j=index.j-radius; j<=index.j+radius; j++){
            auto it = pixel.find(Index(i,j));
            if (it != pixel.end())  // index (i,j) found
                coords_in_area.insert(coords_in_area.end(), it->second.begin(), it->second.end());
        }
    }
    return coords_in_area;
}

vector<GeoCoord*> Pixelation::area_around(const GeoCoord& coord, int radius) const{
    Index index = get_index(coord);
    vector<GeoCoord*> coords_in_area = area_around(index, radius);
    return coords_in_area;
}

vector<GeoCoord*> Pixelation::ring_around(const Index& index, int radius) const{
    vector<GeoCoord*> coords_in_area;
    for (int i=index.i-radius; i<=index.i+radius; i++){
        for (int j=index.j-radius; j<=index.j+radius; j++){
            if (abs(i-index.i)==radius || abs(j-index.j)==radius) {
                auto it = pixel.find(Index(i, j));
                if (it != pixel.end())  // index (i,j) found
                    coords_in_area.insert(coords_in_area.end(), it->second.begin(), it->second.end());
            }
        }
    }
    return coords_in_area;
}

vector<GeoCoord*> Pixelation::ring_around(const GeoCoord& coord, int radius) const{
    Index index = get_index(coord);
    vector<GeoCoord*> coords_in_area = ring_around(index, radius);
    return coords_in_area;
}

vector<GeoCoord*> Pixelation::find_closest_to(const GeoCoord& coord, unsigned int count, const double maxdist) const{
    int radius=1, maxradius=ceil(maxdist/precision);
    double sqrmaxdist = pow(maxdist, 2);
    vector<double> sqrdist;
    vector<GeoCoord*> coords;
    count = min((int)count, (int)pixel.size());     // if count too big
    if (count == 0)     // nothing to do
        return coords;
    vector<GeoCoord*> candidates = area_around(coord, radius);

    while (candidates.size() < count && radius <= maxradius) {
        coords = ring_around(coord, ++radius);
        if (!coords.empty()) {
            candidates.insert(candidates.end(), coords.begin(), coords.end());
            coords.clear();
        }
    }
    if (candidates.empty())                     // no candidates found
        return coords;
    else if (candidates.size() < count)         // not enough candidates
        count = candidates.size();
    else if (radius < maxradius) {
        // now we have a square with count candidates, but they could be in the corners
        // --> create a bigger square
        candidates = area_around(coord, min((int)ceil(1.415*(radius+0.5)-0.5), maxradius));
    }
    for (auto cand : candidates)
        sqrdist.push_back(sqr_distance_2d(coord, *cand));

    vector<unsigned int> indices (sqrdist.size());
    std::iota(indices.begin(), indices.end(), 0);   // creates a vector 0, 1, 2...
    partial_sort(indices.begin(), indices.begin() + count, indices.end(),
                  [&sqrdist](const unsigned int &a, const unsigned int &b){return sqrdist[a] < sqrdist[b];});
    for (int i=0; i<count; i++) {
        if (sqrdist[indices[i]] > sqrmaxdist)
            break;
        else
            coords.push_back(candidates[indices[i]]);
    }
    return coords;
}

void Pixelation::display(GeoCoord::CRS crs, bool display_coords){
    cout << "Origin: " << origin.as_string(crs) << "\n";
    cout << "Bounding Box: ";
    Index bl, tr;
    cout << bottom_left.as_string(crs) << " " << top_right.as_string(crs);
    bl = get_index(bottom_left);
    tr = get_index(top_right);
    cout << "\n" << "Total number of pixels: ";
    cout << (tr.i - bl.i + 1) * (tr.j - bl.j + 1);
    cout << "\n" << "Number of filled pixels: " << pixel.size() << "\n";

    if(display_coords){
        cout << "Coordinates: " << "\n";
        for (const auto& pix : pixel){
            cout << "[" << pix.first.i << ", " << pix.first.j << "]: ";
            for (auto coord : pix.second)
                cout << coord->as_string(crs) << " ";
            cout <<"\n";
        }
    }
}
