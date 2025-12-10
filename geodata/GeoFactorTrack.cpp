//
// Created by cornelia.blanke on 26.06.2024.
//

#include "GeoFactorTrack.h"

GeoFactorTrack::GeoFactorTrack()=default;

GeoFactorTrack::GeoFactorTrack(GeoCoord* pstart, GeoCoord* pend, double factor) : GeoTrack(pstart, pend), penalty(factor) {
    length = penalty * dist;
}

GeoFactorTrack::GeoFactorTrack(GeoCoord& start, GeoCoord& end, double factor) : GeoTrack(start, end), penalty(factor) {
    length = penalty * dist;
}

GeoFactorTrack::~GeoFactorTrack()=default;


void GeoFactorTrack::setPenalty(const double factor) {
    penalty = factor;
}

void GeoFactorTrack::compute_penalty() {
    penalty = 1.0;
}

void GeoFactorTrack::compute_track() {
    length = dist * penalty;
}
