//
// Created by cornelia.blanke on 26.06.2024.
//

#include <fstream>
#include "GeoTrack.h"

GeoTrack::GeoTrack()=default;

GeoTrack::GeoTrack(GeoCoord* pstart, GeoCoord* pend) : start(pstart), end(pend) {
    dist = distance_2d(*start, *end);
    length = dist;
    if (dist == 0){
        cerr << "Invalid Track" << endl;
        exit(1);
    }
}

GeoTrack::GeoTrack(GeoCoord& cstart, GeoCoord& cend) {
    start = &cstart;
    end = &cend;
    dist = distance_2d(*start, *end);
    length = dist;
    if (dist == 0){
        cerr << "Invalid Track" << endl;
        exit(1);
    }
}

GeoTrack::~GeoTrack()=default;

void GeoTrack::get_steps() {
    steps.resize(2);
    steps[0] = start;
    steps[1] = end;
}

void GeoTrack::check() {
    if (steps.empty()) {
        cerr << "Track has no steps." << endl;
        exit(1);
    }
    if (start->distance(*steps.front()) > 1e-6) {
        cerr << "Track start is not equal to the first step." << endl;
        // track->steps.insert(track->steps.begin(), track->start);
        exit(1);
    }
    if (end->distance(*steps.back()) > 1e-6) {
        cerr << "Track end is not equal to the last step." << endl;
        // track->steps.push_back(track->end);
        exit(1);
    }
}

void GeoTrack::compute_track() {
    if (dist < 0){
        cerr << "Distance between start and end not defined" << endl;
        exit(1);
    }
    else if (dist == 0){
        cerr << "Invalid Track" << endl;
        exit(1);
    }
    else
        length = dist;
}

string GeoTrack::write_steps(const GeoCoord::CRS crs){
    if (steps.empty())
        get_steps();
    stringstream ss;
    for (const auto& step : steps)
        ss << step->as_string(crs, '[') << ",";
    string s = std::move(ss).str();
    s.pop_back();      // remove last character
    return s;
}

string GeoTrack::write_feature(const GeoCoord::CRS crs){
    stringstream ss;
    ss << "\n\t  {\n\t\t\"type\": \"Feature\",\n";
    ss << "\t\t\"geometry\": {\n\t\t  \"type\": \"LineString\",\n\t\t  \"coordinates\": [\n\t\t\t";
    ss << write_steps(crs);
    ss << "\n\t\t  ]\n\t\t},\n\t\t\"properties\": {" << write_properties() << "}\n\t  }";
    string s = std::move(ss).str();
    //s.pop_back();      // remove last character
    return s;
}

string GeoTrack::write_properties(){
    stringstream ss;
    ss << "\"distance\":" << dist << ",";
    ss << "\"length\":" << length;
    string s = std::move(ss).str();
    return s;
}

string GeoTrack::write_header() {
    string s;
    s = "{\n  \"type\": \"FeatureCollection\",\n  \"features\": [";
    return s;

}
string GeoTrack::write_closing(){
    string s;
    s = "\n  ]\n}\n";
    return s;
}

void GeoTrack::write_geojson(const char *fileout, const GeoCoord::CRS crs) {
    ofstream out (fileout, ofstream::out);
    out << write_header();
    out << write_feature(crs);
    out << write_closing();
}

void GeoTrack::display() const {
    cout << "Track start:\n";
    start->display();
    cout << "Track end:\n";
    end->display();
    cout << "Distance: " << dist << "\n";
    cout << "Length: " << length << "\n";
}