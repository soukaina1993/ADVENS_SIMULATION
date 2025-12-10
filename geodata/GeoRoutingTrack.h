//
// Created by cornelia.blanke on 26.06.2024.
//

#ifndef FLUIDS_GEOROUTINGTRACK_H
#define FLUIDS_GEOROUTINGTRACK_H

#include "GeoTrack.h"

/// Use Routing on openstreetmap ///
class GeoRoutingTrack : public GeoTrack {
//public:

private:
    bool steps_created=false;

public:
    GeoRoutingTrack();
    GeoRoutingTrack(GeoCoord* pstart, GeoCoord* pend);
    GeoRoutingTrack(GeoCoord& start, GeoCoord& end);
    ~GeoRoutingTrack();


    static string execute_curl(const string& strcurl, short buffersize=256);

    void get_steps() override;
    void compute_track() override;
    void get_steps_from_string(const string& instring);
};


#endif //FLUIDS_GEOROUTINGTRACK_H
