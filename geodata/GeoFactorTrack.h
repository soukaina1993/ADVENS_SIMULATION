//
// Created by cornelia.blanke on 26.06.2024.
//

#ifndef FLUIDS_GEOFACTORTRACK_H
#define FLUIDS_GEOFACTORTRACK_H

#include "GeoTrack.h"

/// The length of a line multiplied by a penalty factor ///
class GeoFactorTrack : public GeoTrack {
public:
    double penalty=1.0;
public:
    GeoFactorTrack();
    GeoFactorTrack(GeoCoord* pstart, GeoCoord* pend, double factor=1.0);
    GeoFactorTrack(GeoCoord& start, GeoCoord& end, double factor=1.0);
    ~GeoFactorTrack();

    void setPenalty(double factor);

    void compute_penalty();  //TODO, if needed

    void compute_track() override;
};


#endif //FLUIDS_GEOFACTORTRACK_H
