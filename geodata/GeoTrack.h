//
// Created by cornelia.blanke on 26.06.2024.
//

#ifndef FLUIDS_GEOTRACK_H
#define FLUIDS_GEOTRACK_H

#include"GeoCoord.h"

/// Base class for Tracking is a straight Line ///
class GeoTrack {
public:
    GeoCoord* start= nullptr;
    GeoCoord* end= nullptr;
    double dist=-1.0, length=-1.0;
    vector<GeoCoord*> steps;         // GeoCoords on track

public:
    GeoTrack();
    GeoTrack(GeoCoord* pstart, GeoCoord* pend);
    GeoTrack(GeoCoord& start, GeoCoord& end);
    ~GeoTrack();

    virtual void get_steps();
    void check();       // check that the first and last step of a track are equal to the start and end of the track
    virtual void compute_track();
    string write_steps(GeoCoord::CRS crs=GeoCoord::CH1903plus);
    string write_feature(GeoCoord::CRS crs=GeoCoord::CH1903plus);
    virtual string write_properties();
    static string write_header();
    static string write_closing();
    void write_geojson(const char *fileout, GeoCoord::CRS crs=GeoCoord::CH1903plus);
    virtual void display() const;
};

#endif //FLUIDS_GEOTRACK_H
