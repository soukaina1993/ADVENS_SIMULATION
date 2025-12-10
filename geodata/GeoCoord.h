//
// Created by cornelia.blanke on 25.06.2024.
//

#ifndef FLUIDS_GEOCOORD_H
#define FLUIDS_GEOCOORD_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include<vector>



using namespace std;


class GeoCoord {
public:
    enum CRS {CH1903plus=0, CH1903=1, WGS84=2};
    //  WGS84 in longitude/latitude decimal degrees

    double north=2.6e6, east=1.2e6, altitude_CH=0.0;      // Zimmerwald geostation near Bern
    double lon=7.43864, lat=46.95108, altitude_wgs84=0.0;


public:
    // Constructors
    GeoCoord();
    GeoCoord(double north, double east, double altitude=0.0);   // default CH1903+
    GeoCoord(GeoCoord::CRS crs, double value1, double value2, double altitude=0.0);

    GeoCoord(const GeoCoord& coord);          // copy

    bool operator==(const GeoCoord& coord) const {
        return ((this->north == coord.north) && (this->east == coord.east) && (this->altitude_CH == coord.altitude_CH));
    }

    // member functions
    string as_string(GeoCoord::CRS crs, char bracket='(', int precision= -1) const;
    void display(GeoCoord::CRS crs) const;
    void display() const;

    void compute_wgs84();
    void compute_CH1903plus();

    double sqr_distance(const GeoCoord& b_coord, bool threedim=false) const;
    double distance(const GeoCoord& b_coord, bool threedim=false) const;
};


/*** outside class ***/
double sqr_distance_2d(const GeoCoord& a, const GeoCoord& b);
double sqr_distance_3d(const GeoCoord& a, const GeoCoord& b);
double distance_2d(const GeoCoord& a, const GeoCoord& b);
double distance_3d(const GeoCoord& a, const GeoCoord& b);




#endif //FLUIDS_GEOCOORD_H
