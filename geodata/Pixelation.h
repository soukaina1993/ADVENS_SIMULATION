//
// Created by cornelia.blanke on 01.07.2024.
//

#ifndef FLUIDS_PIXELATION_H
#define FLUIDS_PIXELATION_H

#include<unordered_map>
#include "GeoCoord.h"



// https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
class Index{
public:
    int i=0; int j=0;

public:
    Index()=default;
    Index(int i, int j) : i(i), j(j) {}
    Index(const Index& copy){
        this->i = copy.i;
        this->j = copy.j;
    }
    void set(int a, int b) {i=a; j=b;}

    Index& operator= (const Index& index)= default;

    bool operator==(const Index& index) const {
        if (this->i == index.i && this->j == index.j) return true;
        else return false;
    }
    struct Hash {
        size_t operator()(const Index& index) const {
            size_t rowHash = std::hash<int>()(index.i);
            size_t colHash = std::hash<int>()(index.j) << 1;
            return rowHash ^ colHash;
        }
    };
};

/*** ------------------------------------------------------------------------------------------------- ***/

class Pixelation {
public:
    GeoCoord origin;     // origin of Pixelation's coordinate system
    GeoCoord bottom_left, top_right;    // Bounding box
    double precision=1.0;
    unordered_map<Index, vector<GeoCoord*>, Index::Hash> pixel;     // points to all GeoCoords in pixel i,j


public:
    Pixelation()=default;
    explicit Pixelation(GeoCoord& myorigin, double precision=1.0);
    explicit Pixelation(vector<GeoCoord>& coords, double precision=1.0);   // uses first entry as origin
    Pixelation(vector<GeoCoord>& coords, GeoCoord& myorigin, double precision=1.0);
    ~Pixelation()=default;

    void compute(vector<GeoCoord>& coords, double prec=1.0);
    void compute(vector<GeoCoord*>& ptr_coords, double prec=1.0);
    vector<GeoCoord*> get_coords();
    void modify_precision(int prec);
    void modify_origin(GeoCoord& origin);


    void add_coord(GeoCoord& coord);
    void add_to_pixelation(const Index& key, GeoCoord* value);
    Index get_index(const GeoCoord& coord) const;

    vector<GeoCoord*> area_around(const Index& index, int radius= 1) const;
    vector<GeoCoord*> area_around(const GeoCoord& coord, int radius= 1) const;
    vector<GeoCoord*> ring_around(const Index& index, int radius= 1) const;
    vector<GeoCoord*> ring_around(const GeoCoord& coord, int radius= 1) const;

    vector<GeoCoord*> find_closest_to(const GeoCoord& coord, unsigned int count=1, double maxdist=1e8) const;

    void display(GeoCoord::CRS crs=GeoCoord::CH1903plus, bool display_coords=false);


};


#endif //FLUIDS_PIXELATION_H
